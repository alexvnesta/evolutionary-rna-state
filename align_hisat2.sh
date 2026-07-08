#!/usr/bin/env bash
# Align prefetched FASTQs with HISAT2 (STAR arm64 build is broken on this host — reads 0 input).
#
# OUTPUT FORMAT POLICY (durable): COMPLETE, LOSSLESS, reference-based CRAM.
#   - Keeps EVERY read: no --no-unal, no MAPQ filter. Unmapped + multi-mapped reads are retained
#     because that is where TE activation, fusions, and viral/repeat signal live (project phenotypes).
#   - A complete CRAM is a full FASTQ replacement: `samtools fastq` regenerates the reads, and it is
#     ~30% of FASTQ size (reference-based, default lossless quality — NO quality binning).
#   - Downstream tools filter at READ TIME:  AEI/quantification -> `samtools view -q 60 <cram>`.
#   - FASTQs are NOT deleted here (raw reads preserved); prune them separately once CRAMs are verified.
# Coordinate-sorted + .crai indexed. Serial (low memory contention).
set -uo pipefail
REPO=/Users/alex/OrchestratedBiosciences/evolutionary-rna-state
IDX=$REPO/results/rnaseq_pilot_hisat2/genome/index/hisat2/GRCh38.primary_assembly.genome
FASTA=$REPO/reference/GRCh38/GRCh38.primary_assembly.genome.fa
OUT=$REPO/results/editing_crams; FQ=$REPO/data/raw/fastq; MAN=$REPO/editing_subset_manifest.csv
mkdir -p "$OUT"; LOG=$OUT/align_hisat2.log; echo "=== hisat2 align -> lossless CRAM $(date) ===" >> "$LOG"
tail -n +2 "$MAN" | cut -d, -f1 | while read samp; do
  cram="$OUT/${samp}.hisat2.cram"
  f1="$FQ/${samp}_1.fastq.gz"; f2="$FQ/${samp}_2.fastq.gz"
  if [[ -s "$cram" ]] && [[ "$(samtools view -c "$cram" 2>/dev/null)" -gt 100000 ]]; then
    echo "[$samp] CRAM ok, skip" >>"$LOG"; continue; fi
  waited=0
  while { [[ ! -s "$f1" ]] || [[ ! -s "$f2" ]] || ! gzip -t "$f1" 2>/dev/null || ! gzip -t "$f2" 2>/dev/null; } && [[ $waited -lt 3000 ]]; do
    sleep 30; waited=$((waited+30))
  done
  if [[ ! -s "$f1" ]] || [[ ! -s "$f2" ]] || ! gzip -t "$f1" 2>/dev/null || ! gzip -t "$f2" 2>/dev/null; then
    echo "[$samp] FASTQ not ready, skip" >>"$LOG"; continue; fi
  echo "[$samp] HISAT2 $(date)" >>"$LOG"
  tmpu="$OUT/${samp}.all.tmp.bam"
  # Stage 1: align, keep ALL reads (mapped, unmapped, multi-mapped) -> unsorted temp BAM
  hisat2 -p 14 -x "$IDX" -1 "$f1" -2 "$f2" 2>>"$LOG" \
    | samtools view -b - > "$tmpu" 2>>"$LOG"
  # Stage 2: coordinate-sort, emit reference-based lossless CRAM, index
  samtools sort -@ 6 -m 3G -O cram --reference "$FASTA" -o "$cram" "$tmpu" 2>>"$LOG"
  n=$(samtools view -c "$cram" 2>/dev/null || echo 0)
  if [[ -s "$cram" ]] && [[ "$n" -gt 100000 ]]; then
    samtools index "$cram"; echo "[$samp] DONE $(date) $n total reads (complete CRAM)" >>"$LOG"
    rm -f "$tmpu"   # keep FASTQs (raw reads); only drop the intermediate BAM
  else echo "[$samp] HISAT2 FAILED ($n reads)" >>"$LOG"; rm -f "$tmpu"; fi
done
echo "=== hisat2 align done $(date) ===" >>"$LOG"
ls "$OUT"/*.hisat2.cram 2>/dev/null | wc -l
