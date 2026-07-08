#!/usr/bin/env bash
# Build the persistent tximeta conda env for nf-core/rnaseq on osx-arm64.
# Why this exists: bioconda's bioconductor-genomeinfodbdata installs a post-link
# script that downloads the data payload from the Bioconductor archive CDN
# (mghp.osn.xsede.org), which is blocked in the sandbox. A solver-built env
# therefore has GenomeInfoDbData registered but its R library dir empty, so
# GenomeInfoDb (and thus tximeta) fails to load. We build the env once and
# R CMD INSTALL GenomeInfoDbData from the current Bioconductor release (served
# directly from bioconductor.org, not the archive CDN).
set -euo pipefail
_REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "$_REPO_ROOT/pipelines/env.sh"
MM="${MAMBA_EXE:-/Users/alex/.claude-science/conda/bin/micromamba}"
TXENV="$_REPO_ROOT/.nextflow_home/persistent_envs/tximeta"

if [ ! -x "$TXENV/bin/R" ]; then
  "$MM" create --yes --prefix "$TXENV" -c bioconda -c conda-forge \
    "bioconductor-tximeta=1.24.0" "bioconductor-genomeinfodbdata=1.2.13"
fi

# Seed GenomeInfoDbData if its library dir is empty
if [ ! -f "$TXENV/lib/R/library/GenomeInfoDbData/DESCRIPTION" ]; then
  TARBALL=/tmp/GenomeInfoDbData.tar.gz
  # current release (served from bioconductor.org, not the blocked archive CDN)
  curl -sL --max-time 120 -o "$TARBALL" \
    "https://bioconductor.org/packages/release/data/annotation/src/contrib/GenomeInfoDbData_1.2.15.tar.gz"
  "$TXENV/bin/R" CMD INSTALL --library="$TXENV/lib/R/library" "$TARBALL"
fi

# sanity check
"$TXENV/bin/R" --vanilla -e 'suppressMessages(library(tximeta)); cat("tximeta env OK\n")'
