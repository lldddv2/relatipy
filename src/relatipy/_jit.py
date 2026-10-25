"""Compile C source into a cached shared library and load it with ctypes.

This private module turns generated C99 (for instance the source of an
``rp_external_acceleration_fn``) into a shared library with the system C
compiler and returns the address of a function in it.  Libraries are cached
in ``$RELATIPY_CACHE_DIR`` or, when it is unset, ``~/.cache/relatipy/jit``,
under a SHA-256 key of the source, the compiler identity, the flags, the
RelatiPy version and the platform.  Files are written to a temporary name and
renamed into place, so concurrent processes never load a partial library.
"""

from __future__ import annotations

import ctypes
import functools
import hashlib
import os
import platform
import shlex
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

__all__ = [
    "CompilationError",
    "CompilerNotFoundError",
    "JitFunction",
    "FLAGS",
    "cache_directory",
    "compile_function",
    "compiler_command",
]

FLAGS: tuple[str, ...] = ("-std=c11", "-O3", "-fPIC", "-shared")
"""Compiler flags of every JIT library; ``-lm`` is appended when linking."""

_LIBRARIES = ("-lm",)


class CompilerNotFoundError(RuntimeError):
    """No C compiler is available to build a run-time force."""


class CompilationError(RuntimeError):
    """The C compiler rejected generated source.

    Attributes
    ----------
    command : tuple of str
        Compiler command that failed.
    stderr : str
        Compiler diagnostics.
    """

    def __init__(self, message: str, command: tuple[str, ...], stderr: str) -> None:
        super().__init__(message)
        self.command = command
        self.stderr = stderr


@dataclass(frozen=True)
class JitFunction:
    """A function in a loaded JIT library.

    The object keeps the ``ctypes`` library handle alive, so ``address``
    stays valid while it is referenced.

    Attributes
    ----------
    name : str
        Symbol name of the function.
    address : int
        Address of the function in the loaded library.
    path : pathlib.Path
        Path of the cached shared library.
    library : ctypes.CDLL
        Loaded library handle.
    """

    name: str
    address: int
    path: Path
    library: ctypes.CDLL = field(repr=False, compare=False)

    def as_ctypes(self, restype, *argtypes):
        """The function as a ``ctypes`` callable with the given prototype."""
        prototype = ctypes.CFUNCTYPE(restype, *argtypes)
        return prototype(self.address)


def cache_directory() -> Path:
    """Directory of cached JIT libraries (created on first use)."""
    configured = os.environ.get("RELATIPY_CACHE_DIR")
    if configured:
        return Path(configured).expanduser()
    return Path.home() / ".cache" / "relatipy" / "jit"


def compiler_command() -> tuple[str, ...]:
    """Compiler command from ``$CC`` (split like a shell word list) or ``cc``.

    Raises
    ------
    CompilerNotFoundError
        If the compiler executable cannot be found.
    """
    words = tuple(shlex.split(os.environ.get("CC", "") or "cc"))
    if not words:
        words = ("cc",)
    if shutil.which(words[0]) is None:
        raise CompilerNotFoundError(
            f"no C compiler found: {words[0]!r} is not on PATH; install a C11 "
            "compiler or set the CC environment variable to compile run-time "
            "forces"
        )
    return words


@functools.lru_cache(maxsize=16)
def _compiler_identity(command: tuple[str, ...], resolved: str) -> str:
    try:
        completed = subprocess.run(
            [*command, "--version"], check=False, capture_output=True, text=True, timeout=60,
        )
        version = completed.stdout.strip() or completed.stderr.strip()
    except (OSError, subprocess.SubprocessError):
        version = ""
    return f"{resolved}\n{version}"


def _relatipy_version() -> str:
    try:
        from importlib.metadata import PackageNotFoundError, version

        return version("relatipy")
    except PackageNotFoundError:
        return "0+unknown"


def cache_key(source: str, command: tuple[str, ...], extra_flags: tuple[str, ...] = ()) -> str:
    """SHA-256 key of source, compiler, flags, version and platform."""
    resolved = shutil.which(command[0]) or command[0]
    digest = hashlib.sha256()
    for part in (
        source,
        _compiler_identity(command, resolved),
        " ".join(command[1:]),
        " ".join(FLAGS + extra_flags + _LIBRARIES),
        _relatipy_version(),
        sys.platform,
        platform.machine(),
        platform.python_implementation(),
    ):
        encoded = part.encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "little"))
        digest.update(encoded)
    return digest.hexdigest()


def _atomic_write(path: Path, data: bytes) -> None:
    descriptor, temporary = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(data)
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def _build(source: str, command: tuple[str, ...], extra_flags: tuple[str, ...], library: Path) -> None:
    directory = library.parent
    directory.mkdir(parents=True, exist_ok=True)
    source_path = library.with_suffix(".c")
    if not source_path.exists():
        _atomic_write(source_path, source.encode("utf-8"))
    descriptor, temporary = tempfile.mkstemp(dir=directory, prefix=f".{library.name}.", suffix=".tmp")
    os.close(descriptor)
    full = (*command, *FLAGS, *extra_flags, str(source_path), "-o", temporary, *_LIBRARIES)
    try:
        completed = subprocess.run(full, check=False, capture_output=True, text=True)
        if completed.returncode != 0:
            raise CompilationError(
                f"C compilation failed with status {completed.returncode}:\n{completed.stderr}",
                full,
                completed.stderr,
            )
        os.replace(temporary, library)
    finally:
        Path(temporary).unlink(missing_ok=True)


def compile_function(
    source: str,
    name: str,
    *,
    cache_dir: str | os.PathLike | None = None,
    include_dirs: tuple[str | os.PathLike, ...] = (),
) -> JitFunction:
    """Compile ``source`` (or reuse its cached build) and load ``name``.

    Parameters
    ----------
    source
        Complete C translation unit.
    name
        Name of the exported function to load.
    cache_dir
        Cache directory; defaults to :func:`cache_directory`.
    include_dirs
        Extra ``-I`` directories; their absolute paths enter the cache key
        (the contents of included files do not).

    Returns
    -------
    JitFunction
        Loaded function with its address.

    Raises
    ------
    CompilerNotFoundError
        If no compiler is available.
    CompilationError
        If the compiler rejects the source.
    AttributeError
        If the library does not export ``name``.
    """
    command = compiler_command()
    directory = Path(cache_dir) if cache_dir is not None else cache_directory()
    extra_flags = tuple(f"-I{Path(item).resolve()}" for item in include_dirs)
    key = cache_key(source, command, extra_flags)
    library = directory / f"rp_{key[:32]}.so"
    if not library.exists():
        _build(source, command, extra_flags, library)
    handle = ctypes.CDLL(str(library))
    function = getattr(handle, name)
    address = ctypes.cast(function, ctypes.c_void_p).value
    return JitFunction(name=name, address=int(address), path=library, library=handle)
