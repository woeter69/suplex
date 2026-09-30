#!/usr/bin/env bash
set -euo pipefail

runner="$1"
source_root="$2"
work_dir="$3"
mkdir -p "$work_dir/problems"
cp "$source_root/data/examples/tiny.mps" "$work_dir/problems/tiny.mps"
printf 'problem,known_optimal\ntiny,3\n' > "$work_dir/reference.csv"

"$runner" "$work_dir/problems" "$work_dir/reference.csv" "$work_dir/results" >/dev/null
grep -q 'tiny,2,2,4,OPTIMAL' "$work_dir/results.csv"
grep -q ',1$' "$work_dir/results.csv"
grep -q '"problem":"tiny"' "$work_dir/results.json"
grep -q '"passed":true' "$work_dir/results.json"
