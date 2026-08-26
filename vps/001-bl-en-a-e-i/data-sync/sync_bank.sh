#!/usr/bin/env bash
# Pull the Kerr orbit bank from the VPS until it is complete, then stop.
#
# `vpsjob sync-all` (the stock 2-minute timer) already mirrors partials, but it
# marks a finished job as synced regardless of whether the transfer actually
# succeeded — a download cut halfway would never be retried. This script keeps
# pulling until the local copy has every shard the run announced, and only then
# writes `.sync-complete`.
#
# Driven by relatipy-datasync.timer; safe to run by hand at any time.
set -uo pipefail

DS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATH="$HOME/.local/bin:$PATH"

JOB="$(cat "$DS/.job-id" 2>/dev/null || true)"
[ -n "$JOB" ] || { echo "sync_bank: $DS/.job-id is missing" >&2; exit 1; }

MARK="$DS/.sync-complete"
[ -f "$MARK" ] && exit 0

OUT="$HOME/vps/results/$JOB/out"

vpsjob pull "$JOB" >/dev/null 2>&1 || exit 0   # no network: try again next tick

[ -f "$OUT/run.json" ] || exit 0
want="$(jq -r '.n_shards // empty' "$OUT/run.json" 2>/dev/null || true)"
[ -n "$want" ] || exit 0

done_shards="$(jq -r '.next_shard // 0' "$OUT/checkpoint.json" 2>/dev/null || echo 0)"
have="$(find "$OUT/shards" -name 'orbits_*.npz' 2>/dev/null | wc -l)"

if [ "$done_shards" = "$want" ] && [ "$have" = "$want" ]; then
  date -u +%Y-%m-%dT%H:%M:%SZ > "$MARK"
  echo "sync_bank: bank complete ($have/$want shards) -> $DS/data"
else
  echo "sync_bank: $have/$want shards on disk (run at $done_shards/$want)"
fi
