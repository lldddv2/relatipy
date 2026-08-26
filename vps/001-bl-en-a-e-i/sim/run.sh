#!/usr/bin/env bash
# Bootstrap the relatipy environment on the VPS and hand over to python.
#
# The venv lives outside the job directory (~/envs/relatipy-001) so that a
# resumed or relaunched job reuses it instead of rebuilding relatipy — which
# means compiling radau_core.c and resolving the aarch64 wheel set again.
# It is rebuilt automatically whenever the shipped source changes.
#
#   ./run.sh main.py --n-orbits 100000
set -euo pipefail

export PATH="$HOME/.local/bin:$HOME/.cargo/bin:$PATH"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONUNBUFFERED=1

ENV_DIR="${RELATIPY_ENV:-$HOME/envs/relatipy-001}"
STAMP="$ENV_DIR/.source-hash"

command -v uv >/dev/null 2>&1 || {
  echo "run.sh: uv is missing on the VPS; run 'vpsjob setup --with-uv'" >&2
  exit 3
}

hash_source() {
  find relatipy-src -type f \( -name '*.py' -o -name '*.c' -o -name '*.toml' \) \
    -exec sha256sum {} + | sort | sha256sum | cut -d' ' -f1
}

want="$(hash_source)"
if [ ! -x "$ENV_DIR/bin/python" ] || [ "$(cat "$STAMP" 2>/dev/null || true)" != "$want" ]; then
  echo "[run.sh] building environment in $ENV_DIR" >&2
  rm -rf "$ENV_DIR"
  uv venv --python 3.11 "$ENV_DIR"
  uv pip install --python "$ENV_DIR/bin/python" ./relatipy-src
  "$ENV_DIR/bin/python" -c "from relatipy.numeric.metrics import Kerr; Kerr(1.0, 0.9)"
  printf '%s' "$want" > "$STAMP"
  echo "[run.sh] environment ready" >&2
fi

exec "$ENV_DIR/bin/python" -u "$@"
