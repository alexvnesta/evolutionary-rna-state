# RNA-seq annotation pipelines — RUNBOOK (Apple Silicon, no Docker)

Processing + annotation of ICB melanoma RNA-seq for the six **evolutionary
RNA-state** phenotypes, on a single Apple Silicon Mac (M-series, 64 GB RAM,
arm64) using **Nextflow + `-profile conda`** — **no Docker, no Rosetta**.

## 0. TL;DR

```bash
cd <repo>
source pipelines/env.sh          # JDK17, Nextflow, STAR, conda plumbing, offline mode
# one-time setup (already done in the build session):
#   pipelines/scripts/fetch_reference.sh      # GRCh38 + GENCODE v46
#   pipelines/scripts/stage_nfcore.sh         # rnaseq/rnasplice/rnafusion tarballs
#   pipelines/scripts/stage_plugins.sh        # nf-schema, nf-validation plugins
# NOTE: run_rnaseq.sh lets nf-core BUILD its own STAR index (~35 GB RAM, ~27 min
#   on first run) and persists it via --save_reference; -resume reuses it after.
#   A separately prebuilt index (build_star_index.sh) is NOT used by default —
#   passing --star_index + --gtf trips STAR's geneInfo.tab error on the
#   TranscriptomeSAM 2-pass step. To force an external index anyway:
#   export STAR_INDEX_OVERRIDE=/path/to/star_index

# base spine (STAR BAMs + Salmon counts):
pipelines/scripts/run_rnaseq.sh   pipelines/rnaseq/samplesheet_pilot.csv  results/rnaseq

# then the phenotype layers (consume the spine BAMs):
pipelines/scripts/run_rnasplice.sh pipelines/rnasplice/samplesheet_selection.csv \
    pipelines/rnasplice/contrastsheet.csv results/rnasplice
nextflow run pipelines/te_erv/main.nf            -profile conda -c pipelines/conf/mac_arm64.config ...
nextflow run pipelines/intron_retention/main.nf  -profile conda -c pipelines/conf/mac_arm64.config ...
nextflow run pipelines/rna_editing/main.nf       -profile conda -c pipelines/conf/mac_arm64.config ...
```

## 1. Toolchain (installed)

| Component | Version | How |
|-----------|---------|-----|
| Nextflow  | 26.04.4 | conda env `nextflow` |
| JDK       | OpenJDK 17 (Zulu, arm64) | conda env `nextflow`; `JAVA_HOME` set by `env.sh` |
| STAR      | 2.7.11b (arm64) | conda env `nextflow` |
| micromamba| 2.5.0 | system; Nextflow's conda solver (via shims) |
| Salmon    | 2.3.1 (arm64) | conda env `rnaio` |
| samtools / fastp | (arm64) | conda env `rnaio` |

**Everything is native arm64 or noarch.** See `pipelines/docs/ARM64_TOOLCHAIN.md`
for the per-tool arm64 availability map.

## 2. Reference

- `reference/GRCh38/GRCh38.primary_assembly.genome.fa` — GENCODE-distributed GRCh38 primary assembly (194 seqs, md5 `49bdb80d21a64dcb16acfc941843356e`)
- `reference/GRCh38/gencode.v46.primary_assembly.annotation.gtf` — GENCODE v46 (3.47 M feature lines)
- `reference/GRCh38/star_index/` — a standalone STAR index (Genome 3 GB + SA 24 GB + SAindex 1.5 GB, `--sjdbOverhang 100`). **Not consumed by run_rnaseq.sh by default** (see §0 note); the pipeline builds and caches its own under the results genome dir. Kept for the TE/ERV subworkflow's dedicated multimap-aware STAR pass and as an optional `STAR_INDEX_OVERRIDE`.
- `reference/te/retro.hg38.v1.transcripts.gtf` — Telescope HERV/L1 locus annotation (28,513 loci)
- All git-ignored (`reference/**`). Checksums in `reference/GRCh38/checksums.json`.

## 3. The six phenotypes → pipeline map

| Phenotype | Pipeline | Status |
|-----------|----------|--------|
| Base align + quant | **nf-core/rnaseq 3.26.0** (STAR + Salmon) | mature |
| Alternative splicing | **nf-core/rnasplice 1.0.4** (rMATS, DEXSeq, edgeR, SUPPA2) | mature |
| Fusion transcripts | **nf-core/rnafusion 4.1.3** (Arriba + STAR-Fusion; FusionCatcher deferred) | mature |
| TE / ERV activation | **custom** `pipelines/te_erv/` (Telescope locus + TEcount family) | authored here |
| Intron retention | **custom** `pipelines/intron_retention/` (featureCounts IR-ratio) | authored here |
| RNA editing | **custom** `pipelines/rna_editing/` (JACUSA2 sites + Alu Editing Index) | authored here |

The three nf-core pipelines all consume the **STAR BAMs from the rnaseq spine**,
as do the three custom subworkflows — so **run rnaseq first**, then fan out.
(TE/ERV additionally does one dedicated multimap-aware STAR pass for the locus branch.)

## 4. Data & cohorts

- Curated subset (`data/manifests/selection_manifest.csv`): **40 pre-treatment samples** (30 Gide + 10 Riaz), balanced response (21 R / 19 N), ~219 GB raw.
- Pilot (`data/manifests/pilot_manifest.json`): **4 Gide samples** (~22 GB) — the end-to-end validation target.
- FASTQs are fetched md5-verified from ENA into `data/raw/fastq/` (git-ignored) by `pipelines/scripts/fetch_fastq.sh`.
- **Disk strategy:** the subset (219 GB) fits the free disk (~697 GB), so bulk download is fine. `run_cohort_batched.sh` bounds peak FASTQ footprint (fetch batch → align → keep BAM+counts → delete FASTQ) and becomes essential if you add the full Gide+Riaz sets (~1.3 TB).

## 5. Apple Silicon / sandbox specifics (why this repo has extra plumbing)

`source pipelines/env.sh` handles all of the below. On an unrestricted Mac with
Docker you would not need most of it — you could `nextflow run nf-core/rnaseq -r 3.26.0 -profile docker`.

1. **No Docker →** `-profile conda`. All tools have native arm64/noarch conda builds (verified), so conda is not a downgrade here.
2. **JDK not on PATH:** the conda `openjdk` doesn't symlink `java` into `bin/`; `env.sh` sets `JAVA_HOME` to `.../envs/nextflow/lib/jvm`.
3. **`~` is restricted:** `env.sh` points `NXF_HOME`, `HOME`, and `CONDARC` into the writable `.nextflow_home/` so Nextflow's assets/plugins cache and micromamba's config resolve.
4. **Plugin registry blocked:** Nextflow 26 resolves plugins via `registry.nextflow.io/api/v1/...`, which the sandbox proxy blocks. Plugins (`nf-schema@2.5.1`, `nf-validation@1.1.3`) are pre-downloaded from GitHub release assets into `.nextflow_home/plugins/` by `stage_plugins.sh`, and `NXF_OFFLINE=true` makes Nextflow use them locally.
5. **`nextflow pull` blocked:** the sandbox forbids creating any `.git` directory, so a git clone of a pipeline fails. Pipelines are downloaded as **release tarballs** (`stage_nfcore.sh`) and run via their local `main.nf` path.
6. **micromamba, not conda:** exposed as `conda`/`mamba` shims (`pipelines/bin/`). The `conda` shim special-cases `conda config --show channels` to emit conda-style YAML (nf-core pipelines probe this at startup; micromamba's output shape differs and would NPE).
7. **rnasplice 1.0.4 config:** uses the legacy `check_max()`/`max_memory` pattern that Nextflow 26's v2 config parser rejects — `run_rnasplice.sh` sets `NXF_SYNTAX_PARSER=v1`.

## 6. Resource model (`pipelines/conf/mac_arm64.config`)

- `executor.local`: 16 CPU / 58 GB ceiling (leaves headroom on the 18-core / 64 GB machine).
- STAR align/index: `memory 40.GB`, **`maxForks 1`** (never two genome-loads at once).
- Salmon/fastp/QC: light, 2–3 parallel forks.
- `resourceLimits = [cpus:16, memory:58.GB]` caps any per-process request.
