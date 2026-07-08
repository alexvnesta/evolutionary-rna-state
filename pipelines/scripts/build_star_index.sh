#!/usr/bin/env bash
# Build a shared STAR genome index once (~35 GB RAM, ~40 min). Reused everywhere.
set -euo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
source "$REPO/pipelines/env.sh"
REF="$REPO/reference/GRCh38"
IDX="$REF/star_index"; mkdir -p "$IDX"
# sjdbOverhang 100 suits ~100-150bp reads (Gide/Riaz HiSeq)
STAR --runMode genomeGenerate --genomeDir "$IDX" \
     --genomeFastaFiles "$REF/GRCh38.primary_assembly.genome.fa" \
     --sjdbGTFfile "$REF/gencode.v46.primary_assembly.annotation.gtf" \
     --sjdbOverhang 100 --runThreadN 12 \
     --limitGenomeGenerateRAM 40000000000
echo "[done] STAR index at $IDX"
