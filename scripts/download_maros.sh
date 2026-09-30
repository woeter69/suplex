#!/usr/bin/env bash
set -euo pipefail

destination="${1:-tests/benchmarks/maros}"
page_url="${MAROS_PAGE_URL:-https://www.doc.ic.ac.uk/~im/}"
mkdir -p "$destination"

work_dir="$(mktemp -d)"
trap 'rm -rf "$work_dir"' EXIT
curl --fail --silent --show-error --location "$page_url" --output "$work_dir/index.html"

# Resolve the site's current QPDATA1/2/3 links instead of assuming a legacy
# archive path. Links may be absolute or relative to the home page.
grep -Eio 'href="[^"]*QPDATA[123][^"]*"' "$work_dir/index.html" |
  sed -E 's/^href="//I; s/"$//' | sort -u > "$work_dir/archives.txt"

count=0
while IFS= read -r href; do
  [ -n "$href" ] || continue
  case "$href" in
    http://*|https://*|ftp://*) url="$href" ;;
    /*) url="https://www.doc.ic.ac.uk$href" ;;
    *) url="${page_url%/}/${href#./}" ;;
  esac
  archive="$work_dir/qpdata-$count.zip"
  curl --fail --silent --show-error --location "$url" --output "$archive"
  unzip -aa -o "$archive" -d "$destination"
  count=$((count + 1))
done < "$work_dir/archives.txt"

if [ "$count" -ne 3 ]; then
  echo "Expected three QPDATA archives, found $count" >&2
  exit 1
fi
echo "Downloaded the Maros-Meszaros QPS collection to $destination"
