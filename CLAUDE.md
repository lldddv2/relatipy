# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Active Work (remove on merge)

- Branch: `feature/paths-bank`
- Tasks: `.project/tasks/paths_bank/tasks.md`
- Goal: `PathsBank` caching for `Geodesic.get_path` under `src/relatipy/numeric/utils/banks/paths_bank/`.
- Resuming a fresh chat: read `tasks.md`, then open the first `stage_N_*.md` whose box is unchecked and continue from there. Each stage file contains goal, files, acceptance criteria, and commit template.

## Language Rule

**Always respond in English**, regardless of the language used in the user's input. This is a hard rule for this project.

## Project Overview

**RelatiPy** is a Python library for symbolic, numerical, and visual exploration of relativistic geometry, focused on Kerr black holes. It works in **geometric units** (`G = c = 1`) with solar mass as the default reference mass.

## Commands

### Install & Build
```bash
pip install -e .           # Install in editable mode + compile C extension (radau_core)
python -m build            # Build wheel for distribution
```

### Tests
```bash
pytest                     # Run all tests
pytest tests/numeric/geodesic/test_yoshida6.py   # Run a single test file
pytest -q                  # Quiet output
```

### Linting & Formatting
```bash
black src/ tests/
isort src/ tests/
flake8 src/ tests/
```

## Architecture

The library has three main namespaces exposed from `src/relatipy/__init__.py`:

```
relatipy/
├── numeric/       # Numerical computations (core)
│   ├── metrics/       # Schwarzschild, Kerr metric tensors + Christoffel symbols
│   ├── coordinates/   # Boyer-Lindquist, orbital elements, Cartesian conversions
│   ├── geodesic/      # Geodesic integrators (test particle trajectories)
│   │   └── integrators/  # Yoshida6 (symplectic), Radau IIA (stiff, C-backed)
│   ├── constants.py   # Physical constants in geometric units
│   └── utils/         # Unit conversion helpers
├── symbolic/      # SymPy/EinsteinPy symbolic calculations
│   ├── metrics/       # Symbolic Kerr metric
│   └── coordinates/   # Symbolic coordinate systems
└── visualization/ # Plotting
    ├── _2D/           # Matplotlib (SciSubplot — serif, inward ticks)
    └── _3D/           # Plotly interactive (PlotKerr, PlotSchwarzschild)
```

**Typical usage flow:**
1. Instantiate a metric (`Schwarzschild` or `Kerr`) from `numeric.metrics`
2. Set initial conditions via `numeric.coordinates` (orbital elements or Boyer-Lindquist)
3. Integrate with `numeric.geodesic.integrators` (Yoshida6 or Radau)
4. Visualize with `visualization._2D.SciSubplot` or `visualization._3D.PlotKerr`

## C Extension

`radau_core.c` implements the Radau IIA order-5 integrator and is compiled during `pip install -e .` via the custom `BuildPy` class in `setup.py`. The compiled shared library (`radau_core.so` on Unix, `.dll` on Windows) is loaded at runtime via `ctypes` with graceful fallback if missing. When modifying the C file, reinstall to recompile.

The Kerr-specific integrators live in `src/relatipy/numeric/geodesic/integrators/kerr/` (currently under active development).

## Commit & Branch Conventions

**Commit prefix tags** (from `docs/development/`):
- `[FEAT]`, `[FIX]`, `[DOC]`, `[REFACTOR]`, `[TEST]`, `[STYLE]`, `[RELEASE]`

**Branch naming**:
- `feature/NUM-ShortDesc-###` (numeric), `feature/SYM-*`, `feature/VIS-*`, `feature/DOC-*`
- `bugfix/*`, `refactor/*`, `docupdate/*`

## Key Physics Concepts

- **Geometric units**: `G = c = 1`; mass `M` sets the length/time scale
- **Boyer-Lindquist coordinates** `(t, r, θ, φ)`: native coordinate system for Kerr
- **Conserved quantities**: energy `E`, z-angular momentum `L_z`, Carter constant `Q`
- **ISCO**: Innermost Stable Circular Orbit — a key observable computed from the metric
- **Constraint projection**: Radau integrator applies Newton iteration to maintain the geodesic constraint `g_μν u^μ u^ν = -1`
