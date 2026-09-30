# Provenance of the DOP853 coefficients

- Upstream project: SciPy.
- Upstream version: `1.17.1`.
- Source file: `scipy/integrate/_ivp/dop853_coefficients.py`.
- Source SHA-256: `3ab62f5b41eee97eec3a1dfb154e7c80d9200bbeca56961c95ec2ffc04463001`.
- Controller reference: `scipy/integrate/_ivp/rk.py`.
- Controller SHA-256: `fa5d6300917f4f949e699b116f59ae147159d5c6147d55d940dc9d3303891056`.
- Retrieved from the project environment on 2026-09-19.
- Licenses: `LICENSE-SCIPY` and `LICENSE-DOP` in this directory.

Vendored artifact checksums:

- `dop853_coefficients.h`: `089697c6ea5749a7ae6c9fb9bcd002a2796ec9aef3481e52f0b2c85164df03d0`.
- `LICENSE-DOP`: `ed9bf58c6d74d3fad9d92d1d67d9bff8141d8ab60de784516b0711364fd43357`.
- `LICENSE-SCIPY`: `221e59f5e910fd7f94e44f0dac77436a11338c285c6346232e4a850a50da0e94`.

## Local modifications

Only the twelve integration stages and the `B`, `E3`, and `E5` vectors needed
for endpoint integration were translated to a C header.  SciPy's NumPy array
construction, Python solver classes, extra dense-output stages, and dense
interpolator coefficients were not copied.  The RelatiPy controller is an
independent C implementation with caller-owned state, a native context, typed
errors, no allocation, and no mutable global data.

The numeric values in `dop853_coefficients.h` must be checked against the
recorded upstream file whenever the upstream source is changed.
