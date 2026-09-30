"""Build the private native evaluator and Python package."""

from pathlib import Path

import numpy
from Cython.Build import cythonize
from setuptools import Extension, find_packages, setup


ROOT = Path(__file__).parent
NATIVE = ROOT / "native" / "src"
MCMC_SOURCES = [
    ROOT / "bindings" / "cython" / "_mcmc_core.pyx",
    NATIVE / "utils" / "numeric.c",
    NATIVE / "utils" / "tensor.c",
    NATIVE / "metric" / "physic" / "kerr.c",
    NATIVE / "metric" / "kerr.c",
    NATIVE / "geodesic" / "physic" / "kerr.c",
    NATIVE / "geodesic" / "physic" / "kerr_observables.c",
    NATIVE / "geodesic" / "initial" / "convert.c",
    NATIVE / "geodesic" / "initial" / "bound.c",
    NATIVE / "geodesic" / "kerr.c",
    NATIVE / "geodesic" / "kerr_mcmc.c",
    NATIVE / "geodesic" / "integrators" / "integrator.c",
    NATIVE / "geodesic" / "integrators" / "dop853.c",
    NATIVE / "geodesic" / "integrators" / "dp45.c",
    NATIVE / "geodesic" / "integrators" / "radau.c",
]

SOLUTION_SOURCES = [
    ROOT / "bindings" / "cython" / "_core.pyx",
    NATIVE / "utils" / "numeric.c",
    NATIVE / "utils" / "tensor.c",
    NATIVE / "metric" / "physic" / "kerr.c",
    NATIVE / "metric" / "kerr.c",
    NATIVE / "metric" / "kerr_properties.c",
    NATIVE / "geodesic" / "solution" / "reconstruct.c",
    NATIVE / "geodesic" / "solution" / "solve.c",
    NATIVE / "geodesic" / "solution" / "preview.c",
    NATIVE / "geodesic" / "initial" / "convert.c",
    NATIVE / "geodesic" / "initial" / "bound.c",
    NATIVE / "geodesic" / "integrators" / "integrator.c",
    NATIVE / "geodesic" / "integrators" / "dop853.c",
    NATIVE / "geodesic" / "integrators" / "dp45.c",
    NATIVE / "geodesic" / "integrators" / "radau.c",
    NATIVE / "geodesic" / "integrators" / "kerr.c",
    NATIVE / "geodesic" / "physic" / "kerr.c",
    NATIVE / "geodesic" / "kerr.c",
]

extensions = [
    Extension(
        name,
        sources=[str(path.relative_to(ROOT)) for path in sources],
        include_dirs=[
            numpy.get_include(),
            str(ROOT / "native" / "include"),
            str(NATIVE),
        ],
        # -O3 keeps IEEE semantics; fast-math style flags stay excluded (architecture §10.2).
        extra_compile_args=["-std=c11", "-O3"],
        libraries=["m"],
    )
    for name, sources in (
        ("relatipy._mcmc_core", MCMC_SOURCES),
        ("relatipy._core", SOLUTION_SOURCES),
    )
]

setup(
    package_dir={"": "src"},
    packages=find_packages("src"),
    ext_modules=cythonize(extensions, compiler_directives={"language_level": "3"}),
)
