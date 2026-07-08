#!/usr/bin/env bash
# Compute genome-wide Alu Editing Index for every STAR BAM in results/editing_bams.
# STAR unique reads have MAPQ 255 -> default --min-mapq 60 is correct here.
set -uo pipefail
REPO=/Users/alex/OrchestratedBiosciences/evolutionary-rna-state
BAMDIR=$REPO/results/editing_bams
FASTA=$REPO/reference/GRCh38/GRCh38.primary_assembly.genome.fa
ALU=$REPO/reference/GRCh38/repeats/alu.hg38.bed6
OUT=$BAMDIR/aei; mkdir -p "$OUT"
LOG=$OUT/aei_batch.log; echo "=== AEI batch $(date) ===" > "$LOG"
for bam in "$BAMDIR"/*.hisat2.sorted.bam; do
  [[ -f "$bam" ]] || continue
  samp=$(basename "$bam" .hisat2.sorted.bam)
  out="$OUT/${samp}.aei.tsv"
  [[ -s "$out" ]] && { echo "[$samp] exists, skip" >>"$LOG"; continue; }
  [[ -f "${bam}.bai" ]] || samtools index "$bam"
  nreads=$(samtools view -c "$bam" 2>/dev/null)
  if [[ "${nreads:-0}" -lt 100000 ]]; then echo "[$samp] TOO FEW READS ($nreads), skip" >>"$LOG"; continue; fi
  echo "[$samp] AEI start $(date) ($nreads reads)" >>"$LOG"
  python "$REPO/pipelines/rna_editing/bin/compute_aei.py" \
     --bam "$bam" --fasta "$FASTA" --alu "$ALU" --sample "$samp" \
     --min-mapq 60 --min-baseq 25 --out "$out" >>"$LOG" 2>&1
  echo "[$samp] AEI done $(date)" >>"$LOG"
done
# merge
python "$REPO/pipelines/rna_editing/bin/merge_aei.py" "$OUT"/*.aei.tsv > "$OUT/cohort_aei.tsv" 2>>"$LOG" || \
  { head -1 "$(ls $OUT/*.aei.tsv | head -1)"; tail -n +2 -q "$OUT"/*.aei.tsv; } > "$OUT/cohort_aei.tsv"
echo "=== AEI batch done $(date) ===" >>"$LOG"
cat "$OUT/cohort_aei.tsv"
