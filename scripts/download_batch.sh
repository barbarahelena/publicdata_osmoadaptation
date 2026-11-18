#!/usr/bin/env bash
# Download one batch of FASTQ files from ENA, resumable and checksum-verified.
# ena_batches.R copies this file into every batch_NN/ directory as download.sh,
# next to that batch's manifest.tsv.
#
#   ./download.sh [target_dir] [parallel_jobs]
#
#   target_dir     where the FASTQ files go (default: ./fastq next to this script)
#   parallel_jobs  simultaneous downloads (default 4; ENA throttles above ~8)
#
# Safe to rerun: finished files (size and md5 verified) are skipped, partial files
# are resumed, and files that fail the md5 check are deleted and retried.
# Exit status is non-zero if any file could not be downloaded and verified.
#
# Needs: curl, awk, and md5sum or md5 (macOS).

set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="$HERE/manifest.tsv"
OUT="${1:-$HERE/fastq}"
JOBS="${2:-4}"
MAX_ATTEMPTS=5

[ -f "$MANIFEST" ] || { echo "manifest not found: $MANIFEST" >&2; exit 2; }
mkdir -p "$OUT" || exit 2
OUT="$(cd "$OUT" && pwd)"

# manifest columns: 1 run, 2 sample, 3 subject, 4 study, 5 income, 6 url, 7 md5, 8 bytes
NFILES=$(($(wc -l < "$MANIFEST") - 1))
NEED=$(awk -F'\t' 'NR>1 {s+=$8} END {printf "%.0f", s}' "$MANIFEST")
AVAIL=$(df -k "$OUT" | awk 'NR==2 {printf "%.0f", $4*1024}')
echo "batch: $(basename "$HERE")  files: $NFILES  size: $((NEED/1000000000)) GB  free: $((AVAIL/1000000000)) GB"
echo "target: $OUT  parallel jobs: $JOBS"
if [ "$AVAIL" -lt "$NEED" ]; then
  echo "WARNING: free space is less than the batch size; finished files will not fit." >&2
  [ -t 0 ] && { read -r -p "continue anyway? [y/N] " a; [ "$a" = "y" ] || exit 1; }
fi

md5_of() {
  if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | cut -d' ' -f1
  else md5 -q "$1"; fi
}

size_of() { wc -c < "$1" | tr -d ' '; }

fetch_one() {
  url="$1"; md5="$2"; bytes="$3"
  f="$OUT/$(basename "$url")"
  attempt=1
  while [ "$attempt" -le "$MAX_ATTEMPTS" ]; do
    if [ -f "$f" ] && [ "$(size_of "$f")" = "$bytes" ]; then
      if [ "$(md5_of "$f")" = "$md5" ]; then echo "ok   $(basename "$f")"; return 0; fi
      echo "md5 mismatch, retrying: $(basename "$f")" >&2
      rm -f "$f"
    fi
    # -C - resumes a partial file; 416 means the file is already complete
    curl -sS -f -L -C - --retry 3 --retry-delay 5 --connect-timeout 30 \
         -o "$f" "https://$url" 2>/dev/null
    rc=$?
    if [ "$rc" -ne 0 ] && [ "$rc" -ne 33 ] && [ "$(size_of "$f" 2>/dev/null || echo 0)" != "$bytes" ]; then
      echo "curl exit $rc, attempt $attempt/$MAX_ATTEMPTS: $(basename "$f")" >&2
      sleep $((attempt * 10))
    fi
    # a file longer than expected cannot be resumed: start over
    if [ -f "$f" ] && [ "$(size_of "$f")" -gt "$bytes" ]; then rm -f "$f"; fi
    attempt=$((attempt + 1))
  done
  if [ -f "$f" ] && [ "$(size_of "$f")" = "$bytes" ] && [ "$(md5_of "$f")" = "$md5" ]; then
    echo "ok   $(basename "$f")"; return 0
  fi
  echo "FAIL $(basename "$f")" >&2
  echo "$url" >> "$OUT/failed_$(basename "$HERE").txt"
  return 1
}
export -f fetch_one md5_of size_of
export OUT MAX_ATTEMPTS HERE

rm -f "$OUT/failed_$(basename "$HERE").txt"
START=$(date +%s)

awk -F'\t' 'NR>1 {print $6, $7, $8}' "$MANIFEST" \
  | xargs -P "$JOBS" -L 1 bash -c 'fetch_one "$@"' _

# final check, independent of what the workers reported
MISSING=0
while IFS=$'\t' read -r run sample subject study income url md5 bytes; do
  f="$OUT/$(basename "$url")"
  if [ ! -f "$f" ] || [ "$(size_of "$f")" != "$bytes" ]; then
    echo "missing or incomplete: $(basename "$url")" >&2
    MISSING=$((MISSING + 1))
  fi
done < <(tail -n +2 "$MANIFEST")

echo "done in $(( $(date +%s) - START )) s; files missing or incomplete: $MISSING"
[ "$MISSING" -eq 0 ] || { echo "rerun this script to retry; see $OUT/failed_$(basename "$HERE").txt" >&2; exit 1; }
