#!/usr/bin/env bash
# Align prefetched FASTQs with HISAT2 (STAR arm64 build is broken on this host — reads 0 input).
# Serial; unique reads (MAPQ 60) kept; coordinate-sorted + indexed for AEI. Deletes FASTQ on success.
set -uo pipefail
REPO=/Users/alex/OrchestratedBiosciences/evolutionary-rna-state
IDX=$REPO/results/rnaseq_pilot_hisat2/genome/index/hisat2/GRCh38.primary_assembly.genome
OUT=$REPO/results/editing_bams; FQ=$REPO/data/raw/fastq; MAN=$REPO/editing_subset_manifest.csv
mkdir -p "$OUT"; LOG=$OUT/align_hisat2.log; echo "=== hisat2 align $(date) ===" >> "$LOG"
tail -n +2 "$MAN" | cut -d, -f1 | while read samp; do
  bam="$OUT/${samp}.hisat2.sorted.bam"
  f1="$FQ/${samp}_1.fastq.gz"; f2="$FQ/${samp}_2.fastq.gz"
  if [[ -s "$bam" ]] && [[ "$(samtools view -c "$bam" 2>/dev/null)" -gt 100000 ]]; then
    echo "[$samp] BAM ok, skip" >>"$LOG"; continue; fi
  waited=0
  while { [[ ! -s "$f1" ]] || [[ ! -s "$f2" ]] || ! gzip -t "$f1" 2>/dev/null || ! gzip -t "$f2" 2>/dev/null; } && [[ $waited -lt 3000 ]]; do
    sleep 30; waited=$((waited+30))
  done
  if [[ ! -s "$f1" ]] || [[ ! -s "$f2" ]] || ! gzip -t "$f1" 2>/dev/null || ! gzip -t "$f2" 2>/dev/null; then
    echo "[$samp] FASTQ not ready, skip" >>"$LOG"; continue; fi
  echo "[$samp] HISAT2 $(date)" >>"$LOG"
  tmpu="$OUT/${samp}.uniq.tmp.bam"
  # Stage 1: align + MAPQ60 filter -> unsorted temp (2-stage pipe, low contention)
  hisat2 -p 14 -x "$IDX" -1 "$f1" -2 "$f2" --no-unal 2>>"$LOG" \
    | samtools view -b -q 60 - > "$tmpu" 2>>"$LOG"
  # Stage 2: sort separately with generous memory
  samtools sort -@ 6 -m 3G -o "$bam" "$tmpu" 2>>"$LOG"
  n=$(samtools view -c "$bam" 2>/dev/null || echo 0)
  if [[ -s "$bam" ]] && [[ "$n" -gt 100000 ]]; then
    samtools index "$bam"; echo "[$samp] DONE $(date) $n uniq reads" >>"$LOG"
    rm -f "$f1" "$f2" "$tmpu"
  else echo "[$samp] HISAT2 FAILED ($n reads)" >>"$LOG"; rm -f "$tmpu"; fi
done
echo "=== hisat2 align done $(date) ===" >>"$LOG"
ls "$OUT"/*.hisat2.sorted.bam 2>/dev/null | wc -l
