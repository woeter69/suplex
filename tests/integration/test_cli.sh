#!/usr/bin/env bash
set -euo pipefail

cli="$1"
source_root="$2"
work_dir="$3"
mkdir -p "$work_dir"

"$cli" --help | grep -q "Usage: suplex"
"$cli" --version | grep -q "Suplex 0.1.0"
"$cli" --algorithm primal "$source_root/data/examples/tiny.mps" | grep -q "Objective: 3"
"$cli" --algorithm primal "$source_root/data/examples/tiny.lp" | grep -q "Status: OPTIMAL"
"$cli" --algorithm primal --solution "$work_dir/tiny.sol" \
  "$source_root/data/examples/tiny.mps" >/dev/null
grep -q "Status: OPTIMAL" "$work_dir/tiny.sol"

if "$cli" --not-an-option value >/dev/null 2>&1; then
  echo "CLI accepted an unknown option" >&2
  exit 1
fi
if "$cli" "$source_root/data/examples/invalid.txt" >/dev/null 2>&1; then
  echo "CLI accepted an unsupported extension" >&2
  exit 1
fi
