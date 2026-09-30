#!/usr/bin/env bash
set -euo pipefail

destination="${1:-tests/benchmarks/netlib}"
base_url="${NETLIB_BASE_URL:-https://www.netlib.org/lp/data}"
mkdir -p "$destination"

work_dir="$(mktemp -d)"
trap 'rm -rf "$work_dir"' EXIT
curl --fail --silent --show-error --location "$base_url/" --output "$work_dir/index.html"
curl --fail --silent --show-error --location "$base_url/emps.c" --output "$work_dir/emps.c"
"${CC:-cc}" -O2 "$work_dir/emps.c" -o "$work_dir/emps"

# Netlib stores these as extensionless, EMPS-compressed files. The official
# index labels each downloadable model with "lang: compressed MPS".
awk '
  /file:/ && match($0, /href="[^"]+"/) {
    candidate = substr($0, RSTART + 6, RLENGTH - 7)
  }
  /compressed MPS/ && candidate != "" {
    print candidate
    candidate = ""
  }
' "$work_dir/index.html" | sort -u > "$work_dir/problems.txt"

while IFS= read -r problem; do
  [ -n "$problem" ] || continue
  echo "Downloading Netlib $problem"
  curl --fail --silent --show-error --location "$base_url/$problem" --output "$work_dir/$problem"
  "$work_dir/emps" "$work_dir/$problem" > "$destination/$problem.mps"
done < "$work_dir/problems.txt"

echo "Downloaded $(wc -l < "$work_dir/problems.txt") Netlib MPS files to $destination"
