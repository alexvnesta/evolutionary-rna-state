#!/usr/bin/env bash
# Parallel genome-wide AEI over all HISAT2 BAMs (4 samples concurrent; each pysam pileup single-thread).
set -uo pipefail
REPO=/Users/alex/OrchestratedBiosciences/evolutionary-rna-state
BAMDIR=$REPO/results/editing_bams
FASTA=$REPO/reference/GRCh38/GRCh38.primary_assembly.genome.fa
ALU=$REPO/reference/GRCh38/repeats/alu.hg38.bed6
OUT=$BAMDIR/aei; mkdir -p "$OUT"
LOG=$OUT/aei_batch.log; echo "=== AEI parallel $(date) ===" > "$LOG"

do_aei() {
  local bam="$1"
  local samp=$(basename "$bam" .hisat2.sorted.bam)
  local out="$OUT/${samp}.aei.tsv"
  [[ -s "$out" ]] && { echo "[$samp] exists" >>"$LOG"; return; }
  echo "[$samp] start $(date '+%T')" >>"$LOG"
  python "$REPO/pipelines/rna_editing/bin/compute_aei.py" \
     --bam "$bam" --fasta "$FASTA" --alu "$ALU" --sample "$samp" \
     --min-mapq 60 --min-baseq 25 --out "$out" >>"$LOG" 2>&1 \
     && echo "[$samp] done $(date '+%T')" >>"$LOG" || echo "[$samp] FAILED" >>"$LOG"
}
export -f do_aei; export REPO FASTA ALU OUT LOG

ls "$BAMDIR"/*.hisat2.sorted.bam | grep -v tmp | xargs -P 4 -I {} bash -c 'do_aei "{}"'

# merge
hdr=$(head -1 "$(ls $OUT/*.aei.tsv | head -1)")
{ echo "$hdr"; for f in "$OUT"/*.aei.tsv; do tail -n +2 "$f"; done; } > "$OUT/cohort_aei.tsv"
echo "=== AEI parallel done $(date) ===" >>"$LOG"
echo "=== cohort_aei.tsv ==="; cat "$OUT/cohort_aei.tsv" | cut -f1-6
