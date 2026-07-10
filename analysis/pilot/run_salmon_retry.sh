#!/usr/bin/env bash
# analysis/pilot/run_salmon_pilot.sh
# De-novo transcript/gene quantification pilot — the modeling session's
# independent, lightweight raw-read arm (complements the sibling STAR/nf-core
# pipeline). Stream-align-delete so peak disk stays ~1 sample (~7 GB).
#
# For each run in the pilot manifest:
#   download paired FASTQ from ENA  ->  verify MD5  ->  salmon quant
#   ->  keep quant.sf + logs  ->  DELETE fastq
#
# Usage: run_salmon_pilot.sh <manifest.csv> <index_dir> <out_dir>
set -euo pipefail

MANIFEST="${1:?manifest csv}"
INDEX="${2:?salmon index dir}"
OUT="${3:?output dir}"
SALMON="${SALMON:-salmon}"
THREADS="${THREADS:-8}"

mkdir -p "$OUT" "$OUT/_fastq_tmp"
LOG="$OUT/pilot_run.log"; : > "$LOG"
echo "run_accession,status,mapping_rate,n_processed,quant_path" > "$OUT/pilot_index.csv"

md5of() { if command -v md5sum >/dev/null; then md5sum "$1" | awk '{print $1}'; else md5 -q "$1"; fi; }

# stream-subsample one mate: curl | zcat | head N*4 | gzip, robustly.
# Returns 0 and writes $2 iff the result is a VALID gzip with exactly N reads.
# (head closing the pipe raises SIGPIPE upstream — expected; we judge success
# by the OUTPUT integrity, not the pipe exit status.)
stream_subsample() {
  local url="$1" dest="$2" nreads="$3"
  ( set +o pipefail
    curl -sS --retry 3 --retry-delay 2 "$url" | zcat 2>/dev/null | head -n $((nreads*4)) | gzip > "$dest" )
  gzip -t "$dest" 2>/dev/null || return 1                 # complete gzip trailer?
  local n; n=$(zcat < "$dest" 2>/dev/null | wc -l); n=$((n/4))
  [ "$n" -ge "${MINREADS:-1500000}" ] || return 2                       # got the full N reads?
  return 0
}

# Read the manifest on FD 3 — NOT stdin — so curl/salmon inside the loop cannot
# consume manifest lines (the classic `while read` stdin-stealing bug that
# truncated a mate and duplicated an iteration in the first run).
NREADS="${NREADS:-3000000}"
while IFS=, read -r cohort run sample title patient arm tp recist resp gb ftp md5 rest <&3; do
  run=$(echo "$run" | tr -d '"')
  [ -z "$run" ] || [ "$run" = "run_accession" ] && continue
  sdir="$OUT/$run"
  if [ -f "$sdir/quant.sf" ]; then
    echo "[$run] already quantified — skip" | tee -a "$LOG"; continue
  fi
  echo "==== [$run] ($resp, $arm, $recist) ====" | tee -a "$LOG"

  IFS=';' read -r f1 f2 <<< "$(echo "$ftp" | tr -d '\"')"
  IFS=';' read -r m1 m2 <<< "$(echo "$md5" | tr -d '\"')"
  r1="$OUT/_fastq_tmp/${run}_1.fastq.gz"
  r2="$OUT/_fastq_tmp/${run}_2.fastq.gz"

  ok=1
  for pair in "$f1|$r1|$m1" "$f2|$r2|$m2"; do
    IFS='|' read -r url dest expmd5 <<< "$pair"
    url="${url#ftp://}"; url="${url#https://}"
    if [ "$NREADS" = "0" ]; then
      echo "  downloading (full) https://$url" | tee -a "$LOG"
      curl -sS --retry 3 --max-time 7200 -o "$dest" "https://$url" </dev/null || { ok=0; break; }
      got=$(md5of "$dest")
      if [ -n "$expmd5" ] && [ "$got" != "$expmd5" ]; then
        echo "  MD5 MISMATCH $dest exp=$expmd5 got=$got" | tee -a "$LOG"; ok=0; break
      fi
    else
      echo "  streaming first ${NREADS} reads https://$url" | tee -a "$LOG"
      if ! stream_subsample "https://$url" "$dest" "$NREADS"; then
        echo "  SUBSAMPLE FAILED (truncated/short) $dest" | tee -a "$LOG"; ok=0; break
      fi
    fi
  done
  if [ "$ok" -ne 1 ]; then
    echo "$run,dl_or_md5_fail,,," >> "$OUT/pilot_index.csv"; rm -f "$r1" "$r2"; continue
  fi

  GENEMAP_ARG=""; [ -n "${GENEMAP:-}" ] && [ -f "${GENEMAP:-}" ] && GENEMAP_ARG="-g $GENEMAP"
  # selective alignment is salmon 2.0's default; --validateMappings is a no-op there.
  "$SALMON" quant -i "$INDEX" -l A -1 "$r1" -2 "$r2" \
      -p "$THREADS" --gcBias --seqBias $GENEMAP_ARG \
      -o "$sdir" </dev/null >> "$LOG" 2>&1 || { echo "$run,quant_fail,,," >> "$OUT/pilot_index.csv"; rm -f "$r1" "$r2"; continue; }

  # robust: pull mapping stats from aux_info/meta_info.json (salmon 2.0)
  meta="$sdir/aux_info/meta_info.json"
  mr=$(grep -o '"percent_mapped"[^,]*' "$meta" 2>/dev/null | grep -o '[0-9.]*' | head -1 || echo "")
  np=$(grep -o '"num_processed"[^,]*' "$meta" 2>/dev/null | grep -o '[0-9]*' | head -1 || echo "")
  echo "$run,ok,${mr:-NA},${np:-NA},$sdir/quant.sf" >> "$OUT/pilot_index.csv"
  echo "  [$run] done mapping_rate=${mr:-NA}%" | tee -a "$LOG"
  rm -f "$r1" "$r2"
done 3< <(tail -n +2 "$MANIFEST")

rmdir "$OUT/_fastq_tmp" 2>/dev/null || true
echo "ALL DONE" | tee -a "$LOG"
