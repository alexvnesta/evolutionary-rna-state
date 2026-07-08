#!/usr/bin/env bash
# Run nf-core/rnafusion on this Mac via conda profile.
# TWO phases: (1) build_references (LARGE, ~1x), (2) fusion detection.
# arm64 note: Arriba has a native osx-arm64 conda build (primary caller).
# STAR-Fusion is perl+noarch (usable; CTAT ref ~30GB). FusionCatcher is
# amd64-oriented and heavy -> DEFERRED (skipped here; enable on cluster).
# Usage:
#   run_rnafusion.sh build   <refdir>
#   run_rnafusion.sh detect  <samplesheet.csv> <outdir> <refdir> [extra args]
set -euo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
source "$REPO/pipelines/env.sh"
MODE="${1:?mode: build|detect}"
REF="$REPO/reference/GRCh38"
COMMON=(
  -profile conda
  -c "$REPO/pipelines/conf/mac_arm64.config"
  --genomes_base_path "$REF"
  --fasta "$REF/GRCh38.primary_assembly.genome.fa"
  --gtf   "$REF/gencode.v46.primary_assembly.annotation.gtf"
  # arm64-viable callers only; skip fusioncatcher (amd64/heavy)
  --arriba --starfusion --stringtie
  --skip_qc false
)
if [ "$MODE" = "build" ]; then
  REFDIR="${2:?refdir}"; WORK="$REPO/results/large/nf_work_rnafusion_ref"
  mkdir -p "$REFDIR" "$WORK"
  echo ">> Building rnafusion references into $REFDIR (large download; hours)"
  nextflow run "$REPO/pipelines/nfcore/rnafusion-4.1.3/main.nf" "${COMMON[@]}" \
    -work-dir "$WORK" --build_references --genomes_base "$REFDIR" \
    --outdir "$REFDIR" -resume
elif [ "$MODE" = "detect" ]; then
  SS="${2:?samplesheet}"; OUT="${3:?outdir}"; REFDIR="${4:?refdir}"; shift 4 || true
  WORK="$REPO/results/large/nf_work_rnafusion"; mkdir -p "$OUT" "$WORK"
  nextflow run "$REPO/pipelines/nfcore/rnafusion-4.1.3/main.nf" "${COMMON[@]}" \
    -work-dir "$WORK" --input "$SS" --outdir "$OUT" --genomes_base "$REFDIR" \
    -resume "$@"
else echo "unknown mode $MODE"; exit 1; fi
