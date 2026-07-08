#!/usr/bin/env bash
# AEI panel batch. Reads the complete lossless CRAMs (results/editing_crams/*.hisat2.cram); the
# MAPQ-60 unique-read filter that AEI needs is applied at READ TIME by compute_aei_fast.py
# (samtools mpileup -q 60), not baked into the alignment. samtools reads CRAM natively given -f ref.
# Falls back to the legacy filtered BAMs (results/editing_bams/*.hisat2.sorted.bam) if CRAMs absent.
set -uo pipefail
REPO=/Users/alex/OrchestratedBiosciences/evolutionary-rna-state
CRAMDIR=$REPO/results/editing_crams; BAMDIR=$REPO/results/editing_bams
OUT=$REPO/results/editing_aei; mkdir -p "$OUT"
FASTA=$REPO/reference/GRCh38/GRCh38.primary_assembly.genome.fa
ALU=$REPO/reference/GRCh38/repeats/alu.panel_chr1_19.bed6
LOG=$OUT/aei_panel.log; echo "=== AEI panel $(date) ===" > "$LOG"
do_aei(){ local aln="$1"; local s=$(basename "$aln" | sed -E 's/\.(hisat2\.cram|hisat2\.sorted\.bam)$//')
  local o="$OUT/${s}.aei.tsv"; [[ -s "$o" ]] && return
  echo "[$s] start $(date '+%T')" >>"$LOG"
  python "$REPO/compute_aei_fast.py" --bam "$aln" --fasta "$FASTA" --alu "$ALU" \
     --sample "$s" --min-mapq 60 --min-baseq 25 --out "$o" >>"$LOG" 2>&1 \
     && echo "[$s] done $(date '+%T')" >>"$LOG" || echo "[$s] FAILED" >>"$LOG"; }
export -f do_aei; export REPO FASTA ALU OUT LOG
if ls "$CRAMDIR"/*.hisat2.cram >/dev/null 2>&1; then
  ls "$CRAMDIR"/*.hisat2.cram | xargs -P 6 -I {} bash -c 'do_aei "{}"'
else
  ls "$BAMDIR"/*.hisat2.sorted.bam | grep -v tmp | xargs -P 6 -I {} bash -c 'do_aei "{}"'
fi
hdr=$(head -1 "$(ls $OUT/*.aei.tsv|head -1)")
{ echo "$hdr"; for f in "$OUT"/*.aei.tsv; do tail -n +2 "$f"; done; } > "$OUT/cohort_aei.tsv"
echo "=== done $(date) ===" >>"$LOG"; cat "$OUT/cohort_aei.tsv"
# NOTE: AEI is depth/format-invariant (it's a rate). Values from complete CRAMs equal those from the
# legacy filtered BAMs because compute_aei_fast.py applies -q 60 at read time either way.
