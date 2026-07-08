#!/usr/bin/env bash
set -uo pipefail
REPO=/Users/alex/OrchestratedBiosciences/evolutionary-rna-state
BAMDIR=$REPO/results/editing_bams; OUT=$BAMDIR/aei; mkdir -p "$OUT"
FASTA=$REPO/reference/GRCh38/GRCh38.primary_assembly.genome.fa
ALU=$REPO/reference/GRCh38/repeats/alu.panel_chr1_19.bed6
LOG=$OUT/aei_panel.log; echo "=== AEI panel $(date) ===" > "$LOG"
do_aei(){ local bam="$1"; local s=$(basename "$bam" .hisat2.sorted.bam)
  local o="$OUT/${s}.aei.tsv"; [[ -s "$o" ]] && return
  echo "[$s] start $(date '+%T')" >>"$LOG"
  python "$REPO/compute_aei_fast.py" --bam "$bam" --fasta "$FASTA" --alu "$ALU" \
     --sample "$s" --min-mapq 60 --min-baseq 25 --out "$o" >>"$LOG" 2>&1 \
     && echo "[$s] done $(date '+%T')" >>"$LOG" || echo "[$s] FAILED" >>"$LOG"; }
export -f do_aei; export REPO FASTA ALU OUT LOG
ls "$BAMDIR"/*.hisat2.sorted.bam | grep -v tmp | xargs -P 6 -I {} bash -c 'do_aei "{}"'
hdr=$(head -1 "$(ls $OUT/*.aei.tsv|head -1)")
{ echo "$hdr"; for f in "$OUT"/*.aei.tsv; do tail -n +2 "$f"; done; } > "$OUT/cohort_aei.tsv"
echo "=== done $(date) ===" >>"$LOG"; cat "$OUT/cohort_aei.tsv"
