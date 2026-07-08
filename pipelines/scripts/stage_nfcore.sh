#!/usr/bin/env bash
# Stage nf-core pipelines as release TARBALLS (not git clones).
# This sandbox blocks creation of any `.git` directory (COARSE git-protection),
# so `nextflow pull nf-core/x` (a git clone) fails with "Operation not permitted".
# Tarballs have no .git and run identically via a local path.
# On an unrestricted machine, instead just: nextflow run nf-core/x -r VER
set -euo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
DEST="$REPO/pipelines/nfcore"; mkdir -p "$DEST"
PIPES="rnaseq:3.26.0 rnasplice:1.0.4 rnafusion:4.1.3"
for spec in $PIPES; do
  name="${spec%%:*}"; ver="${spec##*:}"; d="$DEST/${name}-${ver}"
  if [ -f "$d/main.nf" ]; then echo "[skip] $name $ver present"; continue; fi
  echo "[get ] nf-core/$name $ver"
  curl -fsSL --retry 3 -o "/tmp/${name}.tgz" \
    "https://github.com/nf-core/${name}/archive/refs/tags/${ver}.tar.gz"
  tar xzf "/tmp/${name}.tgz" -C "$DEST" --exclude=".vscode" --exclude=".git*" --exclude=".devcontainer" || true
  [ -f "$d/main.nf" ] && echo "  -> $d/main.nf" || echo "  !! main.nf missing"
done
echo "[done] staged under $DEST"
