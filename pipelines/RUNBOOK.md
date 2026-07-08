# RNA-seq annotation pipelines — RUNBOOK (Apple Silicon, no Docker)

Processing + annotation of ICB melanoma RNA-seq for the six **evolutionary
RNA-state** phenotypes, on a single Apple Silicon Mac (M-series, 64 GB RAM,
arm64) using **Nextflow + `-profile conda`** — **no Docker, no Rosetta**.

> **Aligner note (read first):** the base spine uses **HISAT2**, not STAR. The
> only osx-arm64 conda build of STAR 2.7.11b silently ingests **0 reads** (writes
> empty BAMs, exits 0) — reproduced across all three arm64 builds and the osx-64
> build under Rosetta 2. HISAT2 is arm64-native and verified. Salmon
> (`--pseudo_aligner salmon`) provides transcript/gene quantification. STAR paths
> below are retained only as the amd64/Docker fallback. See §7 for the evidence.

## 0. TL;DR

```bash
cd <repo>
source pipelines/env.sh          # JDK17, Nextflow, conda plumbing, offline mode, GNU tools
# one-time setup (already done in the build session):
#   pipelines/scripts/fetch_reference.sh      # GRCh38 + GENCODE v46
#   pipelines/scripts/stage_nfcore.sh         # rnaseq/rnasplice/rnafusion tarballs
#   pipelines/scripts/stage_plugins.sh        # nf-schema, nf-validation plugins
# run_rnaseq.sh runs nf-core/rnaseq with --aligner hisat2 --pseudo_aligner salmon.
# HISAT2 auto-builds its index on first run; on 64 GB RAM it builds a PLAIN index
# (no splice sites — the splice-aware build needs ~200 GB) and persists it via
# --save_reference; -resume reuses it.

# base spine (HISAT2 genome BAMs + Salmon counts):
pipelines/scripts/run_rnaseq.sh   pipelines/rnaseq/samplesheet_pilot.csv  results/rnaseq

# then the phenotype layers (consume the spine BAMs in results/rnaseq/hisat2/):
pipelines/scripts/run_rnasplice.sh pipelines/rnasplice/samplesheet_selection.csv \
    pipelines/rnasplice/contrastsheet.csv results/rnasplice
# custom subworkflows set conda.enabled in their own nextflow.config — do NOT pass
# -profile conda; pass their config with a second -c instead:
nextflow run pipelines/te_erv/main.nf \
    -c pipelines/conf/mac_arm64.config -c pipelines/te_erv/nextflow.config \
    --input <csv> --genome_fasta <fa> --gene_gtf <gtf> \
    --te_gtf_locus <locus.gtf> --te_gtf_family <family.gtf> --outdir results/te_erv
nextflow run pipelines/intron_retention/main.nf \
    -c pipelines/conf/mac_arm64.config -c pipelines/intron_retention/nextflow.config \
    --input <csv> --gtf <gtf> --outdir results/intron_retention
nextflow run pipelines/rna_editing/main.nf \
    -c pipelines/conf/mac_arm64.config -c pipelines/rna_editing/conf/editing.config \
    --bam_glob 'results/rnaseq/hisat2/*.markdup.sorted.bam' \
    --fasta <fa> --rmsk <rmsk.txt.gz> --outdir results/rna_editing
```

## 1. Toolchain (installed)

| Component | Version | How |
|-----------|---------|-----|
| Nextflow  | 26.04.4 | conda env `nextflow` |
| JDK       | OpenJDK 17 (Zulu, arm64) | conda env `nextflow`; `JAVA_HOME` set by `env.sh` |
| HISAT2    | 2.2.1 (arm64) | per-process conda env (base spine aligner) |
| bowtie2   | 2.5.4 (arm64) | per-process conda env (TE/ERV Telescope aligner) |
| STAR      | 2.7.11b (arm64) | conda env `nextflow` — **broken, unused** (see §7) |
| micromamba| 2.5.0 | system; Nextflow's conda solver (via shims) |
| Salmon    | (arm64) | nf-core module env; `--pseudo_aligner salmon` |
| samtools / fastp | (arm64) | conda env `rnaio` + module envs |

**Everything is native arm64 or noarch.** See `pipelines/docs/ARM64_TOOLCHAIN.md`
for the per-tool arm64 availability map.

## 2. Reference

- `reference/GRCh38/GRCh38.primary_assembly.genome.fa` — GENCODE-distributed GRCh38 primary assembly (194 seqs, md5 `49bdb80d21a64dcb16acfc941843356e`)
- `reference/GRCh38/gencode.v46.primary_assembly.annotation.gtf` — GENCODE v46 (3.47 M feature lines)
- `reference/GRCh38/star_index/` — a standalone STAR index, **unused** on arm64 (STAR is broken; see §7). Kept only for an amd64/Docker fallback. HISAT2 builds and caches its own index under the results genome dir on first run.
- `reference/te/retro.hg38.v1.transcripts.gtf` — Telescope HERV/L1 **locus** annotation (28,513 loci).
- `reference/te/GRCh38_rmsk_TE.gtf` — TEtranscripts **family**-level GTF, built from the UCSC hg38 RepeatMasker table (`pipelines/te_erv/bin/fetch_te_annotation.sh BUILD_FAMILY=1`) when the Hammell-lab curated GTF is unavailable.
- `reference/GRCh38/repeats/rmsk.hg38.txt.gz` + `alu.hg38.bed6` — UCSC RepeatMasker; Alu source for the RNA-editing Alu Editing Index.
- All git-ignored (`reference/**`). Checksums in `reference/GRCh38/checksums.json`.

## 3. The six phenotypes → pipeline map

| Phenotype | Pipeline | Status |
|-----------|----------|--------|
| Base align + quant | **nf-core/rnaseq 3.26.0** (HISAT2 + Salmon on arm64) | mature; validated |
| Alternative splicing | **nf-core/rnasplice 1.0.4** (rMATS, DEXSeq, edgeR, SUPPA2) | mature |
| Fusion transcripts | **nf-core/rnafusion 4.1.3** (Arriba + STAR-Fusion; FusionCatcher deferred) | mature |
| TE / ERV activation | **custom** `pipelines/te_erv/` (Telescope locus via bowtie2 + TEcount family) | authored; family branch (TEcount) run on real data, Telescope locus branch in progress — see status section |
| Intron retention | **custom** `pipelines/intron_retention/` (featureCounts IR-ratio) | authored + integration-tested |
| RNA editing | **custom** `pipelines/rna_editing/` (JACUSA2 sites + Alu Editing Index) | authored + integration-tested |

The three nf-core pipelines and the three custom subworkflows all consume the
**HISAT2 genome BAMs from the rnaseq spine** (`results/rnaseq/hisat2/*.markdup.sorted.bam`),
so **run rnaseq first**, then fan out. Two exceptions in the TE/ERV subworkflow:
the **family** branch (TEcount) uses the spine BAM, but the **locus** branch
(Telescope) needs multimapper-permissive alignment, so it runs its own **bowtie2**
pass (`-k 100 --very-sensitive-local`) from raw FASTQs — Telescope's documented
aligner (Bendall et al. 2019, PLOS Comput Biol; DOI 10.1371/journal.pcbi.1006453).

### Integration status (Gide pilot, PD1_35_PRE, 79.08M-record HISAT2 BAM)

Two of the three custom subworkflows have been run end-to-end against the real
pilot BAM; te_erv is code-complete and parse-verified but not yet executed on
real data.

- **rnaseq spine (HISAT2+Salmon):** RUN, completed; genome BAM 79,080,188 records, Salmon 85,756 nonzero genes, StringTie/featureCounts/bigWig/RSeQC/MultiQC all produced.
- **intron_retention:** RUN end-to-end (exit 0); 289,429 introns, 203,683 evaluated, median IR 0.0076, 26,977 introns IR>0.1; cohort matrix produced.
- **rna_editing:** RUN end-to-end (exit 0, chr21-restricted); 116 A-to-I sites (61 A>G, 55 T>C), Alu Editing Index 0.148% (24 A>G / 16,254 A-cov), cohort_aei.tsv produced. MAPQ floor corrected to 60 (HISAT2 unique-mapper MAPQ; STAR uses 255); optional `--editing_region_bed` confines JACUSA2 + the Alu set to a panel/chromosome.
- **te_erv:** NOT yet run on real data. The Telescope branch was reworked off the broken arm64 STAR onto bowtie2 and the DSL2 parses clean, but BOWTIE2_BUILD / BOWTIE2_ALIGN_MULTI / TELESCOPE_ASSIGN have not been executed against the real FASTQs/BAM. Remaining integration test: bowtie2 index build (~4 GB RAM) + Telescope EM.

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
3. **`~` is restricted:** `env.sh` points `NXF_HOME`, `HOME`, `CONDARC`, and `CONDA_PKGS_DIRS` into the writable `.nextflow_home/` so Nextflow's assets/plugins cache, micromamba's config, and the pkgs cache/lockfile resolve.
4. **Plugin registry blocked:** Nextflow 26 resolves plugins via `registry.nextflow.io/api/v1/...`, which the sandbox proxy blocks. Plugins (`nf-schema@2.5.1`, `nf-validation@1.1.3`) are pre-downloaded from GitHub release assets into `.nextflow_home/plugins/<name>-<ver>/` (the `classes/`+`lib/` must be nested under the versioned dir) by `stage_plugins.sh`, and `NXF_OFFLINE=true` makes Nextflow use them locally.
5. **`nextflow pull` blocked:** the sandbox forbids creating any `.git` directory, so a git clone of a pipeline fails. Pipelines are downloaded as **release tarballs** (`stage_nfcore.sh`), extracted with `--exclude=".git*" --exclude=".vscode" --exclude=".devcontainer"`, and run via their local `main.nf` path. (This repo's own git metadata lives in `.gitmeta/`, not `.git/`, for the same reason — see the note in `.gitignore`.)
6. **micromamba, not conda:** exposed as `conda`/`mamba` shims (`pipelines/bin/`). The `conda` shim special-cases `conda config --show channels` (emit conda-style YAML — nf-core probes this at startup and micromamba's shape would NPE) and `conda info --json` (emit a `conda_prefix`); `pipelines/conda_shim/bin/activate` stands in for the per-env activate script micromamba doesn't create.
7. **macOS BSD sed/awk break nf-core version-capture:** nf-core modules capture tool versions with GNU-syntax `sed`/`grep`. `env.sh` builds a GNU-tools env (`.nextflow_home/gnu_tools`: GNU sed + grep + gzip/zcat **only, not awk** — a GNU gawk shadowed system awk and silently broke conda_prefix extraction), and `process.beforeScript` prepends it plus the conda shims onto every task PATH.
8. **rnasplice 1.0.4 config:** uses the legacy `check_max()`/`max_memory` pattern that Nextflow 26's v2 config parser rejects — `run_rnasplice.sh` sets `NXF_SYNTAX_PARSER=v1`.
9. **CLI boolean flags rejected:** nf-schema strict validation rejects `--flag true` on the command line (parsed as a string). Typed booleans (`save_reference`, `skip_bbsplit`, `skip_linting`, `gencode`, `star_ignore_sjdbgtf`) are set in `params{}` in `mac_arm64.config` instead.
10. **`fq` (FASTQ linter) has no osx-arm64 build:** `skip_linting=true`, and each samplesheet row must set an EXPLICIT strandedness (the `auto` path also needs `fq`). The production sheet records `auto` as intent — resolve it before the real run (see `pipelines/rnaseq/STRANDEDNESS_NOTE.md`).
11. **Qualimap aborts on arm64:** `skip_qualimap=true` (JVM `Abort trap: 6` during plot creation, after fully analysing the BAM — non-essential QC).
12. **Bioconductor arm64 gaps:** `arm64_module_overrides.config` bumps module pins that lack osx-arm64 builds (tximeta 1.20→1.24, summarizedexperiment 1.32→1.36, perl 5.26→5.32, subread 2.0.6→2.0.8, stringtie 2.2.1→2.2.3, ucsc-* 377/469→482). `GenomeInfoDbData` is pre-seeded into a persistent tximeta env (its bioconda post-link CDN download is blocked).
13. **MultiQC `referencing`/`rich` gaps:** the `referencing` conda pkg bundles a `.git` dir that trips git-protection, so a `.git`-stripped rebuild is served from a local file:// channel (`NXF_LOCAL_CHANNEL`); and MULTIQC is pinned to `python=3.12 rich=13 rich-click=1.7.4 multiqc=1.33` (rich 15 broke `rich.panel` submodule access).
14. **te_erv pip installs (Telescope + TEtranscripts):** three layered blocks, all handled in `pipelines/te_erv/environment.yml` + `env.sh`:
    - *PyPI TLS:* pip inside a conda env can't verify the sandbox proxy's cert (`SSLCertVerificationError`, macOS `OSStatus -26276`), even for the allowlisted PyPI. `pipelines/conf/pip.conf` (wired via `PIP_CONFIG_FILE` in `env.sh`) trusts `pypi.org`/`files.pythonhosted.org` to skip the proxy-cert check — the hosts are already network-allowlisted, so reach is not widened. This lets `TEtranscripts==2.2.4` resolve from PyPI.
    - *GitHub blocked for Telescope:* `git+https` fails (no `.git` dir allowed) and a remote tarball URL fails the same TLS check. Telescope v1.0.3 is therefore vendored (`pipelines/te_erv/vendor/`) and installed from a **pre-built arm64 wheel**.
    - *Cython won't compile v1.0.3 as-is:* `calignment.pyx` does `from calignment cimport AlignedPair` — a self-cimport modern Cython (0.29.36) rejects (the sibling `.pxd` is auto-applied). `pipelines/te_erv/vendor/build_telescope_wheel.sh` strips that line and builds the wheel with `--no-build-isolation` (so the conda numpy/cython are used). Rebuild it if the env's Python minor version changes.
15. **Nextflow null-param interpolation:** an unset process param (e.g. `params.bowtie2_extra_args`) interpolates into a task script as the literal string `"null"`, which a tool then treats as a positional argument (bowtie2 wrote its SAM to a file named `null`, leaving the samtools pipe empty). All optional `*_extra_args` in `te_erv.nf` are coalesced (`def x = params.x ?: ""`) before interpolation.

## 6. Resource model (`pipelines/conf/mac_arm64.config`)

- `executor.local`: 16 CPU / 58 GB ceiling (leaves headroom on the 18-core / 64 GB machine).
- HISAT2 align: memory-bounded; on 64 GB it builds a **plain** index (no splice sites — the splice-aware build needs ~200 GB, so it auto-degrades), ~5–6 GB RAM to align.
- bowtie2 (TE/ERV Telescope branch): `memory 16.GB`, **`maxForks 1`**.
- Salmon/fastp/QC: light, 2–3 parallel forks.
- `resourceLimits = [cpus:16, memory:58.GB]` caps any per-process request.

## 7. Why HISAT2, not STAR (the arm64 STAR finding)

The base spine originally used STAR. During pilot validation, STAR was found to
**silently ingest 0 reads** on this machine: `Log.final.out` reported "Number of
input reads: 0" and the output BAM had 0 records (~3 KB), yet STAR exited 0
"finished successfully". Reproduced exhaustively:

- all three osx-arm64 conda builds of STAR 2.7.11b (`haf7d672_6/_7/_8`; 2.7.11b is the latest release per the GitHub releases API);
- the osx-64 build under Rosetta 2 (Mach-O x86_64) — also 0 reads;
- prebuilt, freshly-built, and tiny self-built indexes; gzipped and plain FASTQ; single- and paired-end; `--readFilesCommand zcat`/`gunzip -c`/none; stdin input; all filters disabled.

The same input reads are ingested correctly by **HISAT2** (verified: 100k real
reads → 100% paired) and by **bowtie2**. Salmon pseudo-alignment is also
unaffected (real: 42.2M processed, 85.14% mapped). The spine therefore uses
`--aligner hisat2 --pseudo_aligner salmon`; the TE/ERV Telescope branch uses
bowtie2. STAR remains wired only as the amd64/Docker fallback. UV/pip cannot help
— STAR is a compiled C++ binary, not a Python package.

## 8. Disk: archival CRAM (post-processing)

The persistent output of the spine is the coordinate-sorted markdup BAM
(~3.3 GB/sample). To reclaim disk we archive these to CRAM **after** the
subworkflows have consumed them — NOT as the working format. Rationale:

- Every consumer (nf-core/rnaseq, rnasplice, rnafusion, and the three custom
  subworkflows) reads BAM. Making CRAM the working format would require threading
  the reference FASTA through featureCounts + TEcount and adds CRAM encode/decode
  CPU on every pass. Since the BAM is transient either way, converting *after* the
  fan-out gives the same disk savings with no downstream changes.
- `pipelines/scripts/archive_bam_to_cram.sh <ref.fa> <bam|dir>`:
  `samtools view -C -T <ref>` → index → `quickcheck` → verify → delete BAM
  (unless `--keep`). Verify is a **deep** core-field checksum by default
  (QNAME/FLAG/RNAME/POS/MAPQ/CIGAR/RNEXT/PNEXT/SEQ/QUAL + `|TLEN|`); `--quick`
  drops to record-count + flagstat.
- **Losslessness:** the CRAM round-trip is lossless for everything that carries
  alignment meaning. The only per-record difference is the **TLEN sign** for mate
  pairs at the *same position* (POS==MPOS) — htslib recomputes it on decode with a
  spec-permitted tie-break; magnitude is preserved and no downstream tool depends
  on it. Verified on the pilot BAM: 79,080,188 records, deep-verify identical
  after TLEN-sign normalization, **~53% smaller** (3.34 GB → 1.56 GB).
- Quality scores are retained in full (no binning), so the conversion is
  information-preserving, not lossy-CRAM.
