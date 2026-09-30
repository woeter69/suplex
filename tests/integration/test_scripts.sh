#!/usr/bin/env bash
set -euo pipefail

source_root="$1"
bash -n "$source_root/scripts/download_netlib.sh"
bash -n "$source_root/scripts/download_miplib.sh"
bash -n "$source_root/scripts/download_maros.sh"
grep -q 'emps.c' "$source_root/scripts/download_netlib.sh"
grep -q 'benchmark.zip' "$source_root/scripts/download_miplib.sh"
grep -q 'QPDATA' "$source_root/scripts/download_maros.sh"
