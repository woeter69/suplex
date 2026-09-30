#!/usr/bin/env bash
set -euo pipefail

destination="${1:-tests/benchmarks/miplib}"
base_url="${MIPLIB_BASE_URL:-https://miplib.zib.de/downloads}"
archive_url="${MIPLIB_URL:-$base_url/benchmark.zip}"
page_url="${MIPLIB_PAGE_URL:-https://miplib.zib.de/download}"
mkdir -p "$destination"

archive="$destination/benchmark.zip"
curl --fail --location --continue-at - --output "$archive" "$archive_url"
unzip -o "$archive" -d "$destination"

# Easy-set membership and best-known solutions evolve. Resolve the first
# (latest) version advertised by the official download page.
index_file="$(mktemp)"
trap 'rm -f "$index_file"' EXIT
curl --fail --silent --show-error --location "$page_url" --output "$index_file"
easy_name="$(grep -Eo 'easy-v[0-9]+\.test' "$index_file" | head -n 1 || true)"
solution_name="$(grep -Eo 'miplib2017-v[0-9]+\.solu' "$index_file" | head -n 1 || true)"
if [ -n "$easy_name" ]; then
  curl --fail --location --output "$destination/easy.test" "$base_url/$easy_name"
fi
if [ -n "$solution_name" ]; then
  curl --fail --location --output "$destination/miplib2017.solu" "$base_url/$solution_name"
fi

echo "Downloaded the MIPLIB 2017 benchmark set to $destination"
