# Vendored native material

Third-party-derived native artifacts are isolated in this directory. Each
entry must retain its complete applicable license notices, upstream identity,
version or commit, source hashes, local modifications, and checksums.

| Directory | Upstream | Purpose | License records |
| --- | --- | --- | --- |
| `scipy-dop853/` | SciPy 1.17.1 and the Hairer/Wanner DOP distribution | DOP853 tableau and embedded error coefficients | `LICENSE-SCIPY`, `LICENSE-DOP` |

See each directory's `PROVENANCE.md` before updating or redistributing its
contents. RelatiPy implementation sources remain outside `vendor/` and must
not silently absorb copied third-party code.
