#!/usr/bin/env bash
# Fetch the reference inputs the RNA-editing subworkflow needs:
#   1) UCSC hg38 RepeatMasker table (rmsk.txt.gz) -> Alu source for the AEI
#   2) samtools faidx of the genome FASTA (JACUSA2 + pysam need the .fai)
# Optional (commented): REDIportal known editing sites, and a dbSNP BED for masking.
#
# Usage: pipelines/rna_editing/fetch_alu.sh
set -euo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
REF="$REPO/reference/GRCh38"
REPEATS="$REF/repeats"
mkdir -p "$REPEATS"

# --- 1) UCSC RepeatMasker (hg38). Confirmed reachable (200 OK, 2022-10-19). ---
RMSK_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"
RMSK="$REPEATS/rmsk.hg38.txt.gz"
if [ ! -s "$RMSK" ]; then
  echo "[fetch_alu] downloading $RMSK_URL"
  curl -fL --retry 3 -o "$RMSK" "$RMSK_URL"
fi
echo "[fetch_alu] rmsk: $(du -h "$RMSK" | cut -f1)  $RMSK"

# --- 2) FASTA index (.fai) for JACUSA2 / pysam ---
source "$REPO/pipelines/env.sh" >/dev/null 2>&1 || true
FASTA="$REF/GRCh38.primary_assembly.genome.fa"
SAMTOOLS="${CONDA_ROOT:-/Users/alex/.claude-science/conda}/envs/rnaio/bin/samtools"
if [ ! -s "${FASTA}.fai" ]; then
  echo "[fetch_alu] samtools faidx $FASTA"
  "$SAMTOOLS" faidx "$FASTA"
fi
echo "[fetch_alu] fai ready: ${FASTA}.fai"

# --- 3) (OPTIONAL) REDIportal known A-to-I editing sites -----------------------
# REDIportal host (srv00.recas.ba.infn.it) is NOT on the sandbox network
# allowlist and was NOT reachable at authoring time. To use REDIportal
# known-sites mode (REDItools) or to annotate JACUSA2 calls, download the TABLE1
# (hg38) from http://srv00.recas.ba.infn.it/atlas/download.html on a networked
# machine and place it here, then pass --editing_known_sites:
#   REDIPORTAL="$REPEATS/REDIportal_hg38.txt.gz"
echo "[fetch_alu] NOTE: REDIportal known sites are OPTIONAL and must be fetched"
echo "            manually (host not on sandbox allowlist). AEI does not need them."

# --- 4) (OPTIONAL) dbSNP BED for masking known SNPs ---------------------------
# Editing calls should exclude germline SNPs. Provide a bgzipped+tabixed BED via
# --editing_snp_bed. A common source is dbSNP (NCBI) or the UCSC snpXXX track,
# reduced to BED and indexed with:  bgzip snp.bed && tabix -p bed snp.bed.gz
echo "[fetch_alu] done."
