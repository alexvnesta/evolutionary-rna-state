#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_tecount_uniq_cohort.sh — CONSERVATIVE (biased-low) subfamily-level TE/ERV
# quantification directly from the existing HISAT2 unique-mapper BAMs, with NO
# FASTQ refetch and NO multi-mapping realignment.
#
# WHY --mode uniq: our spine BAMs are already unique-only (MAPQ>=60, -F 256), so
# the multi-mapping reads Telescope/TEtranscripts-multi need are gone. TEcount
# --mode uniq counts only unique reads -> subfamily/family-level TE counts. Per
# the 2026 methods review (pipelines/docs/te_erv_quant_methods_2026.md) this is
# a DEFENSIBLE lower-bound feature that SYSTEMATICALLY UNDERCOUNTS the young,
# active HERV-K/HERVH/L1 elements (their reads are multimappers). It is NOT a
# substitute for the locus-level Telescope readout (that needs the refetch).
#
# Usage: run_tecount_uniq_cohort.sh <samplesheet.csv> <outdir> [NPAR]
#   samplesheet cols: sample,fastq_1,fastq_2,bam,strandedness  (only sample+bam used)
# ---------------------------------------------------------------------------
set -uo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
SS="${1:?samplesheet}"; OUT="${2:?outdir}"; NPAR="${3:-6}"
ENV="$REPO/.nextflow_home/conda_cache/env-026b6c24470ddd246b2ec5531ffe713d"
GENE_GTF="$REPO/reference/GRCh38/gencode.v46.primary_assembly.annotation.gtf"
TE_GTF="$REPO/reference/te/GRCh38_rmsk_TE.gtf"
BIN="$ENV/bin"
mkdir -p "$OUT"
# make OUT absolute: run_one does `cd "$OUT"` then checks "$OUT/<sample>.cntTable";
# with a relative OUT that check resolves to $OUT/$OUT/... after the cd, giving a
# false [FAIL] even though the cntTable was written correctly. Absolute OUT fixes it.
OUT="$(cd "$OUT" && pwd)"

# reverse-stranded library -> TEcount --stranded reverse
run_one() {
  local sample="$1" bam="$2"
  local out="$OUT/${sample}.cntTable"
  if [ -s "$out" ]; then echo "[skip] $sample (exists)"; return 0; fi
  cd "$OUT"
  "$BIN/TEcount" --mode uniq --stranded reverse --sortByPos \
    --GTF "$GENE_GTF" --TE "$TE_GTF" \
    -b "$bam" --project "$sample" \
    > "${sample}.tecount.log" 2>&1
  if [ -s "$OUT/${sample}.cntTable" ]; then
    echo "[ok  ] $sample ($(wc -l < "$OUT/${sample}.cntTable") features)"
  else
    echo "[FAIL] $sample"; return 1
  fi
}
export -f run_one; export OUT BIN GENE_GTF TE_GTF

echo "[tecount] $(date '+%H:%M:%S') start uniq-mode, NPAR=$NPAR, out=$OUT"
tail -n +2 "$SS" | tr -d '\r' | while IFS=, read -r sample fq1 fq2 bam strand; do
  while [ "$(jobs -rp | wc -l)" -ge "$NPAR" ]; do sleep 3; done
  run_one "$sample" "$bam" &
done
wait
NDONE=$(ls "$OUT"/*.cntTable 2>/dev/null | wc -l | tr -d ' ')
echo "[tecount] $(date '+%H:%M:%S') done: $NDONE/40 cntTables"
