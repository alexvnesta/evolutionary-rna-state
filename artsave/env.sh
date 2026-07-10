#!/usr/bin/env bash
_REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)"
# Source before any Nextflow invocation on this Apple Silicon Mac.
#   source pipelines/env.sh
# Sets JAVA_HOME to the conda-provided JDK (the managed env doesn't symlink
# java into bin/), and puts the conda `nextflow` env tools on PATH.
_CONDA_ROOT="${CONDA_ROOT:-/Users/alex/.claude-science/conda}"
_NF_ENV="$_CONDA_ROOT/envs/nextflow"
export JAVA_HOME="$_NF_ENV/lib/jvm"
export PATH="$JAVA_HOME/bin:$_NF_ENV/bin:$PATH"
# Nextflow tuning for a 64 GB / 18-core single machine
export NXF_OPTS='-Xms1g -Xmx4g'          # cap Nextflow's own JVM heap
export NXF_ANSI_LOG=false
echo "[env] JAVA_HOME=$JAVA_HOME"
echo "[env] nextflow=$(nextflow -version 2>&1 | awk '/version/{print $2; exit}')"
echo "[env] STAR=$(STAR --version 2>/dev/null)"

# Nextflow needs a writable home for its assets/plugins cache; ~ may be restricted.
export NXF_HOME="${NXF_HOME:-$_REPO_ROOT/.nextflow_home}"

# Run offline: this sandbox blocks registry.nextflow.io's /api/v1 plugin path and
# the remote nf-core custom-config include. Plugins are pre-staged (stage_plugins.sh),
# so offline mode uses them locally. On an unrestricted machine, unset NXF_OFFLINE.
export NXF_OFFLINE=true

# --- Conda/micromamba plumbing for Nextflow's conda profile ---
# This machine's package manager is micromamba; expose it as conda+mamba via shims,
# and redirect HOME so ~/.condarc (sandbox-blocked) resolves into the writable NXF_HOME.
export PATH="$_REPO_ROOT/pipelines/bin:$PATH"     # conda/mamba shims -> micromamba
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$_CONDA_ROOT}"
export MAMBA_EXE="$_CONDA_ROOT/bin/micromamba"
mkdir -p "$NXF_HOME"
export HOME="$NXF_HOME"                            # so ~/.condarc is writable
export CONDARC="$NXF_HOME/.condarc"
[ -f "$CONDARC" ] || printf 'channels: [conda-forge, bioconda]\nchannel_priority: flexible\n' > "$CONDARC"
# Cache conda envs in one place so pipelines reuse them
export NXF_CONDA_CACHEDIR="$_REPO_ROOT/.nextflow_home/conda_cache"
mkdir -p "$NXF_CONDA_CACHEDIR"

# micromamba's default pkgs cache is the managed conda root (read-only here) —
# redirect it to the writable NXF_HOME so env creation can write/lock.
export CONDA_PKGS_DIRS="$_REPO_ROOT/.nextflow_home/pkgs_cache"
mkdir -p "$CONDA_PKGS_DIRS"

# GNU coreutils env (sed/grep/awk) for nf-core version-capture snippets that
# assume GNU syntax — macOS ships BSD variants. Prepended per-task via the
# process.beforeScript in mac_arm64.config. Created by:
#   micromamba create -p .nextflow_home/gnu_tools -c conda-forge sed grep coreutils gawk
export NXF_GNU_TOOLS="$_REPO_ROOT/.nextflow_home/gnu_tools/bin"

# The conda/mamba shim dir must be on each TASK's PATH too (Nextflow's launcher
# runs `conda info --json` to find the activate script). Exposed to the
# process.beforeScript in mac_arm64.config.
export NXF_SHIM_BIN="$_REPO_ROOT/pipelines/bin"

# Persistent tximeta env (pre-seeded with GenomeInfoDbData, whose bioconda
# post-link data download from the Bioconductor archive CDN is sandbox-blocked).
# Built by pipelines/scripts/build_tximeta_env.sh. Referenced by the
# TXIMETA_TXIMPORT withName override in arm64_module_overrides.config.
export NXF_TXIMETA_ENV="$_REPO_ROOT/.nextflow_home/persistent_envs/tximeta"

# Local conda channel holding repackaged, sandbox-safe builds (currently the
# `referencing` package with its .git test-fixture stripped — that .git dir
# breaks micromamba extraction under the sandbox's git-protection). Built by
# pipelines/scripts/build_local_channel.sh. Prepended to conda.channels in
# mac_arm64.config at highest priority.
export NXF_LOCAL_CHANNEL="file://$_REPO_ROOT/.nextflow_home/local_channel"

# pip config for conda-env pip installs (e.g. te_erv Telescope/TEtranscripts).
# The sandbox's TLS-terminating proxy presents a cert the conda envs' pip
# cannot verify -> SSLCertVerificationError (macOS OSStatus -26276), even for
# the allowlisted PyPI. pip.conf trusts pypi.org + the wheel host to skip that
# proxy-cert check (the hosts are already network-allowlisted, so reach is not
# widened). Exported so Nextflow task shells inherit it.
export PIP_CONFIG_FILE="$_REPO_ROOT/pipelines/conf/pip.conf"
