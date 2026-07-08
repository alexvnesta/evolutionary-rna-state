#!/usr/bin/env bash
# Parallel prefetch of all subset FASTQs with per-mate gzip-integrity check + resume.
# Network-bound; runs 4 concurrent. Align step is separate (serial, RAM-bound).
set -uo pipefail
REPO=/Users/alex/OrchestratedBiosciences/evolutionary-rna-state
FQ=$REPO/data/raw/fastq; mkdir -p "$FQ"
MAN=$REPO/editing_subset_manifest.csv
LOG=$FQ/prefetch.log; echo "=== prefetch $(date) ===" > "$LOG"

fetch_mate() {  # $1=url $2=outfile
  local u="$1" f="$2"
  for attempt in 1 2 3 4 5; do
    if [[ -s "$f" ]] && gzip -t "$f" 2>/dev/null; then echo "[ok] $f" >>"$LOG"; return 0; fi
    curl -sS -L --retry 5 --retry-delay 5 -C - "https://${u}" -o "$f" 2>>"$LOG" || true
  done
  if [[ -s "$f" ]] && gzip -t "$f" 2>/dev/null; then return 0; fi
  echo "[FAIL] $f" >>"$LOG"; rm -f "$f"; return 1
}
export -f fetch_mate; export LOG

# build a job list: one line per mate "url outfile"
JOBS=$(mktemp)
tail -n +2 "$MAN" | while IFS=, read -r line; do
  samp=$(echo "$line" | cut -d, -f1)
  ftp=$(echo "$line"  | awk -F, '{print $NF}')
  u1=$(echo "$ftp" | cut -d';' -f1); u2=$(echo "$ftp" | cut -d';' -f2)
  echo "$u1 $FQ/${samp}_1.fastq.gz" >> "$JOBS"
  echo "$u2 $FQ/${samp}_2.fastq.gz" >> "$JOBS"
done
echo "total mate-files: $(wc -l < "$JOBS")" >>"$LOG"

# 4 concurrent downloads
xargs -P 4 -n 2 bash -c 'fetch_mate "$0" "$1"' < "$JOBS"
rm -f "$JOBS"
echo "=== prefetch done $(date) ===" >>"$LOG"
# integrity summary
ok=0; bad=0
tail -n +2 "$MAN" | cut -d, -f1 | while read samp; do :; done
for f in "$FQ"/*.fastq.gz; do if gzip -t "$f" 2>/dev/null; then ok=$((ok+1)); else bad=$((bad+1)); fi; done
echo "valid gz: $ok  bad: $bad" >>"$LOG"
cat "$LOG" | tail -5
