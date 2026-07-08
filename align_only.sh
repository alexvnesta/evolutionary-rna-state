#!/usr/bin/env bash
# Align prefetched FASTQs serially with STAR (RAM-bound: one genome load at a time).
# Deletes each sample's FASTQ after a successful BAM to bound disk.
set -uo pipefail
REPO=/Users/alex/OrchestratedBiosciences/evolutionary-rna-state
IDX=$REPO/reference/GRCh38/star_index
OUT=$REPO/results/editing_bams; FQ=$REPO/data/raw/fastq; MAN=$REPO/editing_subset_manifest.csv
mkdir -p "$OUT"; LOG=$OUT/align.log; echo "=== align-only $(date) ===" >> "$LOG"
tail -n +2 "$MAN" | cut -d, -f1 | while read samp; do
  bam="$OUT/${samp}.Aligned.sortedByCoord.out.bam"
  f1="$FQ/${samp}_1.fastq.gz"; f2="$FQ/${samp}_2.fastq.gz"
  if [[ -s "$bam" ]] && [[ "$(samtools view -c "$bam" 2>/dev/null)" -gt 100000 ]]; then
    echo "[$samp] BAM ok, skip" >>"$LOG"; continue; fi
  # wait for prefetch to deliver both valid mates (up to ~40 min)
  waited=0
  while { [[ ! -s "$f1" ]] || [[ ! -s "$f2" ]] || ! gzip -t "$f1" 2>/dev/null || ! gzip -t "$f2" 2>/dev/null; } && [[ $waited -lt 2400 ]]; do
    sleep 30; waited=$((waited+30))
  done
  if [[ ! -s "$f1" ]] || [[ ! -s "$f2" ]] || ! gzip -t "$f1" 2>/dev/null || ! gzip -t "$f2" 2>/dev/null; then
    echo "[$samp] FASTQ not ready after wait, skip" >>"$LOG"; continue; fi
  echo "[$samp] STAR $(date)" >>"$LOG"
  STAR --runThreadN 14 --genomeDir "$IDX" --readFilesIn "$f1" "$f2" --readFilesCommand gunzip -c \
       --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$OUT/${samp}." \
       --outSAMattributes NH HI AS nM MD --outFilterMultimapNmax 1 \
       --outSAMprimaryFlag AllBestScore --limitBAMsortRAM 20000000000 >>"$LOG" 2>&1
  n=$(samtools view -c "$bam" 2>/dev/null || echo 0)
  if [[ -s "$bam" ]] && [[ "$n" -gt 100000 ]]; then
    samtools index "$bam"; echo "[$samp] DONE $(date) $n reads" >>"$LOG"
    rm -f "$f1" "$f2"; rm -rf "$OUT/${samp}._STARtmp"
  else echo "[$samp] STAR FAILED ($n reads)" >>"$LOG"; fi
done
echo "=== align-only done $(date) ===" >>"$LOG"
ls -la "$OUT"/*.Aligned.sortedByCoord.out.bam 2>/dev/null | wc -l
