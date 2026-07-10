#!/usr/bin/env bash
# Full-depth salmon quant via the AWS SRA mirror (sra-pub-run-odp) — ~5x faster
# than ENA for these runs. Validated equivalent to the ENA/FASTQ path on
# SRR3184280 (identical read count 50,084,443; identical mapping 96.9348%;
# gene-TPM Pearson 1.0, max abs diff 2.06/1e6).
#
# Per run: curl .sra from AWS mirror -> fasterq-dump --split-files -> salmon quant -> delete.
# The AWS mirror serves .sra (not FASTQ), so the manifest FASTQ-MD5 check does
# not apply; integrity is enforced by (a) HTTP full-size download, (b) fasterq-dump
# refusing a corrupt .sra, (c) recording read count + mapping rate for QC.
# Peak disk ~ .sra (~4GB) + 2 FASTQ (~28GB) ~= 32GB/run — higher than the ENA
# stream path, so MIN_FREE_MB default is raised to 60GB.
#
# Usage: MANIFEST INDEX OUTDIR   (env: SALMON FASTERQ THREADS GENEMAP MIN_FREE_MB)
set -uo pipefail
MANIFEST="$1"; INDEX="$2"; OUT="$3"
SALMON="${SALMON:-salmon}"; FASTERQ="${FASTERQ:-fasterq-dump}"
THREADS="${THREADS:-8}"; GENEMAP="${GENEMAP:-}"
TMP="$OUT/_sra_tmp"; mkdir -p "$TMP"
LOG="$OUT/fulldepth_aws.log"; : > "$LOG"
[ -f "$OUT/pilot_index.csv" ] || echo "run_accession,status,mapping_rate,n_processed,quant_path" > "$OUT/pilot_index.csv"
MIN_FREE_MB="${MIN_FREE_MB:-60000}"
free_mb() { df -m "$OUT" | tail -1 | awk '{print $4}'; }

while IFS=, read -r cohort run sample title patient arm tp recist resp gb ftp md5 rest <&3; do
  sdir="$OUT/$run"
  if [ -s "$sdir/quant.sf" ]; then echo "  [$run] exists, skip" | tee -a "$LOG"; continue; fi
  waited=0
  while [ "$(free_mb)" -lt "$MIN_FREE_MB" ]; do
    echo "  [disk-guard] free $(free_mb)MB < ${MIN_FREE_MB}MB — waiting 120s (waited ${waited}s)" | tee -a "$LOG"
    sleep 120; waited=$((waited+120))
    if [ "$waited" -ge 7200 ]; then echo "  [disk-guard] gave up after 2h; stopping" | tee -a "$LOG"; break 2; fi
  done
  echo "==== [$run] ($resp) full-depth via AWS mirror ====" | tee -a "$LOG"
  sra="$TMP/$run.sra"
  aws_url="https://sra-pub-run-odp.s3.amazonaws.com/sra/$run/$run"
  echo "  downloading .sra $aws_url" | tee -a "$LOG"
  if ! curl -sS --retry 5 --retry-delay 3 --retry-all-errors --max-time 3600 \
        --speed-time 120 --speed-limit 2000 -o "$sra" "$aws_url" 2>>"$LOG"; then
    echo "  DL FAIL $sra" | tee -a "$LOG"; echo "$run,dl_fail,,," >> "$OUT/pilot_index.csv"; rm -f "$sra"; continue
  fi
  echo "  .sra $(stat -f%z "$sra" 2>/dev/null || stat -c%s "$sra") bytes; converting" | tee -a "$LOG"
  if ! "$FASTERQ" --split-files -e "$THREADS" -O "$TMP" "$sra" >>"$LOG" 2>&1; then
    echo "  FASTERQ FAIL $run" | tee -a "$LOG"; echo "$run,fasterq_fail,,," >> "$OUT/pilot_index.csv"; rm -f "$sra" "$TMP/${run}"*.fastq; continue
  fi
  r1="$TMP/${run}_1.fastq"; r2="$TMP/${run}_2.fastq"
  rm -f "$sra"   # free the .sra as soon as FASTQ exists
  if [ ! -s "$r1" ] || [ ! -s "$r2" ]; then
    echo "  MISSING MATES $run" | tee -a "$LOG"; echo "$run,fasterq_fail,,," >> "$OUT/pilot_index.csv"; rm -f "$TMP/${run}"*.fastq; continue
  fi
  # memory guard: salmon selective-alignment needs several GB; wait for a RAM window
  # (avoids OOM kills when sibling sessions spike memory on this shared box)
  MIN_FREE_MEM_MB="${MIN_FREE_MEM_MB:-8000}"
  free_mem_mb() { vm_stat 2>/dev/null | awk '/Pages free/{f=$3} /Pages inactive/{i=$3} END{gsub(/\./,"",f); gsub(/\./,"",i); print int((f+i)*16384/1048576)}'; }
  mwaited=0
  while [ "$(free_mem_mb)" -lt "$MIN_FREE_MEM_MB" ]; do
    echo "  [mem-guard] free $(free_mem_mb)MB < ${MIN_FREE_MEM_MB}MB — waiting 90s (waited ${mwaited}s)" | tee -a "$LOG"
    sleep 90; mwaited=$((mwaited+90))
    if [ "$mwaited" -ge 5400 ]; then echo "  [mem-guard] gave up after 90min" | tee -a "$LOG"; break; fi
  done
  GENEMAP_ARG=""; [ -n "$GENEMAP" ] && [ -f "$GENEMAP" ] && GENEMAP_ARG="-g $GENEMAP"
  "$SALMON" quant -i "$INDEX" -l A -1 "$r1" -2 "$r2" -p "$THREADS" --gcBias --seqBias $GENEMAP_ARG \
      -o "$sdir" </dev/null >> "$LOG" 2>&1 || { echo "$run,quant_fail,,," >> "$OUT/pilot_index.csv"; rm -f "$TMP/${run}"*.fastq; continue; }
  meta="$sdir/aux_info/meta_info.json"
  mr=$(grep -o '"percent_mapped"[^,]*' "$meta" 2>/dev/null | grep -o '[0-9.]*' | head -1 || echo "")
  np=$(grep -o '"num_processed"[^,]*' "$meta" 2>/dev/null | grep -o '[0-9]*' | head -1 || echo "")
  echo "$run,ok,${mr:-NA},${np:-NA},$sdir/quant.sf" >> "$OUT/pilot_index.csv"
  echo "  [$run] done mapping_rate=${mr:-NA}% n=${np:-NA}" | tee -a "$LOG"
  rm -f "$TMP/${run}"*.fastq
done 3< <(tail -n +2 "$MANIFEST")
rmdir "$TMP" 2>/dev/null || true
echo "ALL DONE" | tee -a "$LOG"
