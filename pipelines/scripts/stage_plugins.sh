#!/usr/bin/env bash
# Manually stage Nextflow plugins into $NXF_HOME/plugins.
# Why: Nextflow 26.x routes plugin resolution through registry.nextflow.io's
# /api/v1/... path, which this sandbox's proxy blocks ("Operation not permitted")
# even with the host allowlisted. The plugin ZIPs are GitHub release assets
# (reachable), so we install them offline and run with NXF_OFFLINE=true.
# On an unrestricted machine this whole script is unnecessary — Nextflow
# auto-downloads plugins on first run.
set -euo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
source "$REPO/pipelines/env.sh" >/dev/null
PLUGDIR="$NXF_HOME/plugins"; mkdir -p "$PLUGDIR"
# plugin:version:github_repo
SPECS="nf-schema:2.5.1:nextflow-io/nf-schema nf-validation:1.1.3:nextflow-io/nf-validation"
for spec in $SPECS; do
  name="${spec%%:*}"; rest="${spec#*:}"; ver="${rest%%:*}"; ghrepo="${rest##*:}"
  target="$PLUGDIR/${name}-${ver}"
  if [ -f "$target/classes/META-INF/MANIFEST.MF" ]; then echo "[skip] $name-$ver present"; continue; fi
  echo "[get ] $name-$ver from $ghrepo"
  url="https://github.com/${ghrepo}/releases/download/${ver}/${name}-${ver}.zip"
  curl -fsSL --retry 3 --max-time 180 -o "/tmp/${name}-${ver}.zip" "$url"
  tmp="$(mktemp -d)"; unzip -oq "/tmp/${name}-${ver}.zip" -d "$tmp"
  rm -rf "$target"; mkdir -p "$target"
  # zip contains classes/ + lib/ at top level -> nest under versioned dir
  mv "$tmp"/classes "$tmp"/lib "$target/" 2>/dev/null || mv "$tmp"/* "$target/"
  rm -rf "$tmp"
  [ -f "$target/classes/META-INF/MANIFEST.MF" ] && echo "  -> $target OK" || echo "  !! manifest missing in $target"
done
echo "[done] plugins staged in $PLUGDIR"; ls -d "$PLUGDIR"/*/ 2>/dev/null
