# Python tests

Each directory answers one question about the code under test.

```text
tests/
├── api/                 Public contracts, grouped like src/relatipy/
│   ├── coordinates/     containers, units, immutability, identity
│   ├── metrics/         Kerr value object
│   ├── geodesic/        State, Orbit, Solution, integration options and errors
│   ├── observables/     KerrMcmcModel and the Kerr observable evaluator
│   └── plotting/        static, interactive and temporal plots
├── internal/            Private helpers and direct calls to relatipy._core
├── reference/           Scientific validation against independent references
├── fixtures/            Frozen JSON references and validation reports
└── support/             Shared test-only builders and paths
```

## Where a new test goes

1. It checks a native computation through a peer library, invariants or an
   exact solution: `reference/`.
2. It imports a private name (`relatipy._core`, `relatipy._validation`, a
   `_`-prefixed module) as the subject under test: `internal/`.
3. Otherwise: `api/<subpackage>/`, following the public module it exercises.

A test that uses a private module only to observe or patch a public call (for
example `monkeypatch` on `_core`) stays in `api/`.

## Rules

- Shared builders live in `support/`; never import one test module from
  another. Import them as `from support.builders import make_state`.
- Test file names are unique across the tree; `api/` subdirectories have no
  `__init__.py`.
- `reference/` and `fixtures/` keep their paths: frozen reports record module
  names (`reference.*`), relative paths and SHA-256 checksums of these files.
  Moving them invalidates that provenance.
- Native C tests live in `../native/tests/`; see `../native/README.md`.

## Running

`pyproject.toml` sets `testpaths` and puts `tests/` on `sys.path`, so run from
the repository root:

```console
uv run --with pytest pytest                   # whole suite
uv run --with pytest pytest tests/api         # public contracts only
uv run --with pytest pytest tests/api/geodesic/test_state.py
```

Peer-library tests in `reference/` skip when KerrGeoPy or PyGRO is missing.
