"""The private JIT compiles C into a cached shared library."""

from __future__ import annotations

import ctypes
from pathlib import Path

import pytest

from relatipy import _jit

SOURCE = """
#include <math.h>
double rp_test_hypot(double x, double y) { return sqrt(x * x + y * y); }
"""


def _has_compiler() -> bool:
    try:
        _jit.compiler_command()
    except _jit.CompilerNotFoundError:
        return False
    return True


pytestmark = pytest.mark.skipif(not _has_compiler(), reason="no C compiler")


def _call(function: _jit.JitFunction, x: float, y: float) -> float:
    return function.as_ctypes(ctypes.c_double, ctypes.c_double, ctypes.c_double)(x, y)


def test_compile_loads_function_address_and_name(tmp_path: Path) -> None:
    function = _jit.compile_function(SOURCE, "rp_test_hypot", cache_dir=tmp_path)
    assert function.name == "rp_test_hypot"
    assert isinstance(function.address, int) and function.address != 0
    assert function.path.parent == tmp_path and function.path.suffix == ".so"
    assert _call(function, 3.0, 4.0) == 5.0


def test_cache_is_reused_and_keyed_by_source(tmp_path: Path) -> None:
    first = _jit.compile_function(SOURCE, "rp_test_hypot", cache_dir=tmp_path)
    stamp = first.path.stat().st_mtime_ns
    second = _jit.compile_function(SOURCE, "rp_test_hypot", cache_dir=tmp_path)
    assert second.path == first.path
    assert second.path.stat().st_mtime_ns == stamp
    other = _jit.compile_function(SOURCE + "\n/* changed */\n", "rp_test_hypot", cache_dir=tmp_path)
    assert other.path != first.path
    assert not [p for p in tmp_path.iterdir() if p.name.endswith(".tmp")]


def test_cache_key_covers_compiler_flags_version_and_platform(monkeypatch) -> None:
    command = _jit.compiler_command()
    key = _jit.cache_key(SOURCE, command)
    assert len(key) == 64
    assert _jit.cache_key(SOURCE + " ", command) != key
    assert _jit.cache_key(SOURCE, command, ("-I/somewhere",)) != key
    assert _jit.cache_key(SOURCE, (*command, "-m64")) != key
    monkeypatch.setattr(_jit, "_relatipy_version", lambda: "999.0")
    assert _jit.cache_key(SOURCE, command) != key
    monkeypatch.undo()
    monkeypatch.setattr(_jit.sys, "platform", "other-platform")
    assert _jit.cache_key(SOURCE, command) != key


def test_cache_directory_from_environment(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setenv("RELATIPY_CACHE_DIR", str(tmp_path / "jit-cache"))
    assert _jit.cache_directory() == tmp_path / "jit-cache"
    function = _jit.compile_function(SOURCE, "rp_test_hypot")
    assert function.path.parent == tmp_path / "jit-cache"
    monkeypatch.delenv("RELATIPY_CACHE_DIR")
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    assert _jit.cache_directory() == tmp_path / "home" / ".cache" / "relatipy" / "jit"


def test_missing_compiler_raises_a_clear_error(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setenv("CC", "relatipy-no-such-compiler")
    with pytest.raises(_jit.CompilerNotFoundError, match="CC environment variable"):
        _jit.compile_function(SOURCE, "rp_test_hypot", cache_dir=tmp_path)


def test_compiler_from_environment_is_used(monkeypatch, tmp_path: Path) -> None:
    command = _jit.compiler_command()
    monkeypatch.setenv("CC", " ".join(command))
    assert _jit.compiler_command() == command


def test_compilation_error_reports_diagnostics(tmp_path: Path) -> None:
    with pytest.raises(_jit.CompilationError) as error:
        _jit.compile_function("int broken(void) { return }", "broken", cache_dir=tmp_path)
    assert error.value.stderr
    assert "-std=c11" in error.value.command and "-O3" in error.value.command
    assert not list(tmp_path.glob("*.so"))
    assert not [p for p in tmp_path.iterdir() if p.name.endswith(".tmp")]


def test_missing_symbol_raises_attribute_error(tmp_path: Path) -> None:
    with pytest.raises(AttributeError):
        _jit.compile_function(SOURCE, "rp_test_absent", cache_dir=tmp_path)
