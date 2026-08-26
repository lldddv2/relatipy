"""Build hook: compile ctypes shared library radau_core for the Radau Kerr integrator."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

from setuptools import setup
from setuptools.command.bdist_wheel import bdist_wheel as _bdist_wheel
from setuptools.command.build_py import build_py as _build_py

ROOT = Path(__file__).parent.resolve()
RADAU_SRC = (
    ROOT
    / "src"
    / "relatipy"
    / "numeric"
    / "geodesic"
    / "integrators"
    / "kerr"
    / "radau_core.c"
)
INTEGRATORS_REL = Path("relatipy") / "numeric" / "geodesic" / "integrators" / "kerr"


def _lib_name() -> str:
    return "radau_core.dll" if sys.platform == "win32" else "radau_core.so"


def _compile_radau_core(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / _lib_name()
    if not RADAU_SRC.is_file():
        raise RuntimeError(f"radau_core.c not found at {RADAU_SRC}")
    if sys.platform == "win32":
        cc = os.environ.get("CC", "gcc")
        cmd = [cc, "-O3", "-shared", "-o", str(out), str(RADAU_SRC)]
    else:
        cc = os.environ.get("CC", "cc")
        cmd = [cc, "-std=c11", "-O3", "-fPIC", "-shared", "-o", str(out), str(RADAU_SRC), "-lm"]
    subprocess.check_call(cmd)


class BuildPy(_build_py):
    def run(self) -> None:
        super().run()
        self._build_radau()

    def _build_radau(self) -> None:
        dests: list[Path] = []
        if self.build_lib:
            dests.append(Path(self.build_lib) / INTEGRATORS_REL)
        dests.append(ROOT / "src" / INTEGRATORS_REL)
        seen: set[Path] = set()
        for d in dests:
            key = d.resolve()
            if key in seen:
                continue
            seen.add(key)
            _compile_radau_core(d)


class RelatipyBdistWheel(_bdist_wheel):
    """radau_core.so/.dll is native; wheel tag must be platform-specific (not py3-none-any)."""

    def finalize_options(self) -> None:
        super().finalize_options()
        self.root_is_pure = False


setup(
    cmdclass={
        "build_py": BuildPy,
        "bdist_wheel": RelatipyBdistWheel,
    },
)
