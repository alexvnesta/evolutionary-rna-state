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
| Alternative splicing | **nf-core/rnasplice 1.0.4** (DEXSeq + edgeR on arm64; rMATS/SUPPA amd64-only) | mature; validated (genome_bam) |
| Fusion transcripts | **nf-core/rnafusion 4.1.3** (Arriba + STAR-Fusion) | mature; **arm64-blocked** (needs STAR) |
| TE / ERV activation | **custom** `pipelines/te_erv/` (Telescope locus via bowtie2 + TEcount family) | authored + integration-tested |
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

All three custom subworkflows have now been run end-to-end against the real
pilot data (Gide PD1_35_PRE).

- **rnaseq spine (HISAT2+Salmon):** RUN, completed; genome BAM 79,080,188 records, Salmon 85,756 nonzero genes, StringTie/featureCounts/bigWig/RSeQC/MultiQC all produced.
- **intron_retention:** RUN end-to-end (exit 0); 289,429 introns, 203,683 evaluated, median IR 0.0076, 26,977 introns IR>0.1; cohort matrix produced.
- **rna_editing:** RUN end-to-end (exit 0, chr21-restricted); 116 A-to-I sites (61 A>G, 55 T>C), Alu Editing Index 0.148% (24 A>G / 16,254 A-cov), cohort_aei.tsv produced. MAPQ floor corrected to 60 (HISAT2 unique-mapper MAPQ; STAR uses 255); optional `--editing_region_bed` confines JACUSA2 + the Alu set to a panel/chromosome.
- **te_erv:** RUN end-to-end (exit 0). bowtie2 multimap align (42.25M pairs, 75.4% overall), Telescope EM converged (200 iters, log-likelihood ~2.04e7, 28,513-feature locus report) → **locus matrix 15,209 loci**; TEcount family branch → **family matrix 1,328 TE families + 63,140 genes**. Outputs in `results/te_erv_integration_test/matrices/`. Telescope v1.0.3 is a patched arm64 wheel (self-cimport + `np.int` fixes, see §5.14); bowtie2 replaces the broken arm64 STAR multimap pass.

The two mature nf-core splice/fusion pipelines were also pilot-run:

- **rnasplice (DEXSeq + edgeR):** RUN end-to-end (exit 0) on a 3-responder vs
  3-nonresponder Gide contrast via **`--source genome_bam`** — the HISAT2 spine
  BAMs feed rnasplice directly, bypassing the arm64-broken STAR. DEXSeq: 701,811
  exon bins tested, 74 significant (padj<0.05), 9,535 genes q-valued. edgeR: 16,730
  genes tested for differential exon usage (top hit TTN/ENSG00000155657, FDR 7e-82).
  Outputs in `results/rnasplice_pilot/star_salmon/{dexseq_exon,edger}/`. Run it with
  `SOURCE=genome_bam bash pipelines/scripts/run_rnasplice.sh <ss> <contrasts> <out>`;
  `pipelines/conf/rnasplice_arm64_overrides.config` remaps ~10 conda pins with no
  osx-arm64 build (samtools, htseq, subread, edgeR, gffread, dexseq→1.56, rsem
  minus STAR, sed/grep). **rMATS and SUPPA are opt-in / amd64-only** — rMATS's
  `rmats + r-pairadise` env has an unsatisfiable R-toolchain conflict on arm64
  (r-base 4.2 vs 4.5, gsl, libgfortran); enable with `RMATS=1` only on amd64/Docker.
  SUPPA needs a `salmon_results` source (fastq path), not genome_bam.
- **rnafusion:** NOT runnable natively on arm64. Both fusion callers wired in
  (Arriba, STAR-Fusion) route through `STAR_ALIGN`, and STAR-Fusion additionally
  consumes STAR's `Chimeric.out.junction` — there is no non-STAR fusion path in
  4.1.3. arm64 STAR ingests 0 reads (§7), and osx-64 STAR under Rosetta does the
  same, so both callers would see an empty chimeric BAM. Fusion detection needs an
  amd64 host (native Linux or a Linux VM/container with a working STAR); the runner
  and references are ready for that target. See §7.

## 4. Data & cohorts

- Curated subset (`data/manifests/selection_manifest.csv`): **40 pre-treatment samples** (30 Gide + 10 Riaz), balanced response (21 R / 19 N), ~219 GB raw.
- Pilot (`data/manifests/pilot_manifest.json`): **4 Gide samples** (~22 GB) — the end-to-end validation target.
- FASTQs are fetched into `data/raw/fastq/` (git-ignored) by `pipelines/scripts/fetch_fastq.sh`. **Fetch source: AWS Open Data SRA mirror first, ENA HTTPS fallback** (see §5.16).
- **Disk strategy:** the subset (219 GB) fits the free disk, so bulk download is fine. `run_cohort_batched.sh` bounds peak FASTQ footprint (fetch batch → align → keep BAM+counts → delete FASTQ) and becomes essential if you add the full Gide+Riaz sets (~1.3 TB). **Disk is the binding constraint at batch-size 4, not a comfortable one.** When all 4 samples hit the concurrent alignment→sort→markdup peak together, free disk fell to **~4 GB** in the batch-3 run (an active near-`No space left` crisis, not a graceful trough), and it does **not** self-recover until the first BAM of the batch finalizes. Two manual levers were needed to survive it: (a) delete a batch's FASTQs the moment HISAT2 alignment finishes (they're consumed by the sort step, so waiting for full publication is too late), and (b) in extremis, prune stale `results/large/nf_work_<name>` task scratch older than the current batch (reclaimed 1.39 TB of already-published BAM scratch mid-run). **Prefer batch-size 2–3, or free ≥300 GB before launching, if you want the run to survive unattended.**

### Cohort spine status — COMPLETE (40/40)

All 40 curated samples are aligned and quantified. Published under `results/rnaseq_cohort/`:
- `hisat2/*.markdup.sorted.bam` (+ `.bai`) — **40/40**, ~142 GB total, all `samtools quickcheck` OK.
- `salmon/<sample>/quant.sf` — **40/40**, 253,221 transcripts each.
- `cohort_gene_counts_40samples.tsv` — **unified 40-sample gene-count matrix** (62,679 genes × 40 samples; 30 Gide + 10 Riaz; 21 R / 19 N; per-sample library size 34.6–47.7 M, median 39.7 M). Built by summing salmon `NumReads` per gene via `salmon.merged.tx2gene.tsv`. Saved as an artifact (results/ is git-ignored). This is the matrix that feeds DE + the phenotype subworkflows.
- The batched runner produces a **per-batch** `salmon.merged.gene_counts.tsv`; the single 40-sample matrix above is assembled from the per-sample quants, not from those partial merges.
- The final batch's aggregate MultiQC did not complete (the background run was killed by a host-mount drop mid-QC, §5.17) — cosmetic only; every per-sample BAM + quant is intact.

### Cohort phenotype results — 40 samples (this is what feeds the model)

All results below were run on the 40-sample spine and saved as artifacts (`results/` is git-ignored). Runner scripts + method-choice reviews are committed under `pipelines/scripts/` and `pipelines/docs/`.

- **Alternative splicing — SUPPA2 (primary, DONE).** Event-level dPSI from the existing 40 Salmon quants via `pipelines/scripts/run_suppa_cohort.sh` (arm64-native, no realignment). 306,514 events; PSI matrix `suppa_cohort_psi_40samples.psi`. Differential R-vs-N computed by per-event **Mann-Whitney + BH-FDR** (`suppa_diffsplice_mannwhitney.csv`) because SUPPA's own `diffSplice -gc` compresses p-values on this cohort — **3,427 events at p<0.05 & |dPSI|≥0.10** (AF 1679, SE 546, AL 504, A5 235, A3 234, RI 120, MX 109); nothing survives genome-wide FDR<0.05 (expected at n=40 — the value is the dPSI feature layer). ADAR surfaces as a top differential AF event. 2026 tool-landscape review: `pipelines/docs/splicing_tool_landscape_2026.md` (SUPPA2 = qualified-yes primary; add **SpliceWiz** for novel/cryptic junctions — no osx-arm64 conda build, needs from-source build or Linux; rMATS-turbo Linux-only via r-pairadise).
- **RNA editing — Alu Editing Index (DONE, 40/40).** `pipelines/scripts/run_aei_standalone.sh` — pooled A>G/(A+G) over Alu, strand-aware (reverse lib), MAPQ 60, base-Q 25. Full 1.18M-Alu genome-wide is ~5.4 h/sample (per-interval pysam pileup); ran a **deterministic 200k-Alu representative subset** (seed 42, all 24 chroms; `AEI_FULL=1` for the full set) for ~comparable cross-sample estimates. Result `cohort_aei_40samples.tsv`: mean AEI 0.23% (0.14–0.57%), S/N 2.28. **No R-vs-N burden difference (p=0.65)** — global editing burden does not stratify response; motivates per-site calling. 2026 review `pipelines/docs/aei_editing_methods_2026.md`: our estimator IS the reference AEI algorithm; production upgrade = single `samtools mpileup -l alu.bed` pass (~10× faster, full genome-wide) validated once vs RNAEditingIndexer (Docker/Linux); complement with **REDItools3** (pip, arm64) for per-site + REDIportal; add a SNP mask (dbSNP/gnomAD) and confirm BAQ-off.
- **TE/ERV — TEcount --mode uniq (INTERIM floor, DONE, 40/40).** `pipelines/scripts/run_tecount_uniq_cohort.sh` — `TEcount --mode uniq --stranded reverse --sortByPos` on the existing unique-mapper BAMs, **no FASTQ refetch**. Result `cohort_te_erv_counts_40samples.tsv` (1,328 TE subfamilies × 40; 616 LTR/ERV families) + `te_erv_diff_R_vs_N.tsv`. Coherent signal: **14/20 top differential families are LTR/ERV and all up in non-responders** (of those top-20, ~6 reach nominal p<0.05, lowest ≈0.026; none survive FDR at n=40 — the value is the directional consistency, not per-family significance). This is a deliberately CONSERVATIVE (biased-low) feature — unique-mapper counting undercounts exactly the young HERV-K/HERVH/L1 elements that drive ERV antigenicity (their reads are multimappers, discarded by the spine's `-F 256`/MAPQ≥60). 2026 review `pipelines/docs/te_erv_quant_methods_2026.md`: **the locus-level readout the hypothesis needs requires the ~127 GB FASTQ refetch + bowtie2/STAR `-k 100` multimap realign → Telescope** (still the 2025 benchmark #1; our validated `pipelines/te_erv/` path). Ship the uniq floor now; schedule Telescope on refetch/Linux.
- **Intron retention (DONE, earlier).** featureCounts IR-ratio, 40-sample matrix.
- **Exon-level splicing — DEXSeq/edgeR (rnasplice) (INCOMPLETE).** All 40 DEXSEQ_COUNT + featureCounts finished (preserved in `results/rnasplice_cohort/dexseq_inputs/`), but the final cohort-level `DEXSEQ_EXON` R test stalled/died twice at `estimateDispersions` (759k exon bins × 40) under the host-mount instability (§5.17) — both the nextflow task and a standalone rerun froze at the same point. SUPPA2 is the primary cohort splicing deliverable; DEXSeq is complementary and pilot-validated (74 sig bins, TTN top). Re-run the standalone `run_dexseq_exon.R` on a stable host from the preserved counts, or treat as Linux work.
- **Fusion transcripts (DEFERRED).** rnafusion arm64-blocked (§7b); Linux/amd64 work.

**Deferred-to-Linux/amd64 bucket:** rnafusion (STAR-Fusion), rMATS-turbo, SpliceWiz, the definitive Telescope locus-ERV run, DEXSeq cohort exon test, and the one-time RNAEditingIndexer AEI concordance check.

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
10. **`fq` (FASTQ linter) has no osx-arm64 build → strandedness must be EXPLICIT:** `skip_linting=true`, and every samplesheet row sets `strandedness=reverse` (NOT `auto`). nf-core's `auto` runs `FASTQ_SUBSAMPLE_FQ_SALMON`, which needs `fq=0.12.0` — no osx-arm64 build, so the env solve fails and the whole run dies. RESOLVED: the pilot (ERR2208952) empirically inferred **reverse** (Salmon `expected_format` ISR, RSeQC 92.3% antisense, strand-check status `pass`); Gide (PRJEB23709) + Riaz (PRJNA356761) are TruSeq stranded, all reverse. `_ss_from_manifest.py` writes `reverse` for the cohort; the per-sample strand-check in MultiQC still flags any deviant sample.
11. **Qualimap aborts on arm64:** `skip_qualimap=true` (JVM `Abort trap: 6` during plot creation, after fully analysing the BAM — non-essential QC).
12. **Bioconductor arm64 gaps:** `arm64_module_overrides.config` bumps module pins that lack osx-arm64 builds (tximeta 1.20→1.24, summarizedexperiment 1.32→1.36, perl 5.26→5.32, subread 2.0.6→2.0.8, stringtie 2.2.1→2.2.3, ucsc-* 377/469→482). `GenomeInfoDbData` is pre-seeded into a persistent tximeta env (its bioconda post-link CDN download is blocked).
13. **MultiQC `referencing`/`rich` gaps:** the `referencing` conda pkg bundles a `.git` dir that trips git-protection, so a `.git`-stripped rebuild is served from a local file:// channel (`NXF_LOCAL_CHANNEL`); and MULTIQC is pinned to `python=3.12 rich=13 rich-click=1.7.4 multiqc=1.33` (rich 15 broke `rich.panel` submodule access).
14. **te_erv pip installs (Telescope + TEtranscripts):** three layered blocks, all handled in `pipelines/te_erv/environment.yml` + `env.sh`:
    - *PyPI TLS:* pip inside a conda env can't verify the sandbox proxy's cert (`SSLCertVerificationError`, macOS `OSStatus -26276`), even for the allowlisted PyPI. `pipelines/conf/pip.conf` (wired via `PIP_CONFIG_FILE` in `env.sh`) trusts `pypi.org`/`files.pythonhosted.org` to skip the proxy-cert check — the hosts are already network-allowlisted, so reach is not widened. This lets `TEtranscripts==2.2.4` resolve from PyPI.
    - *GitHub blocked for Telescope:* `git+https` fails (no `.git` dir allowed) and a remote tarball URL fails the same TLS check. Telescope v1.0.3 is therefore vendored (`pipelines/te_erv/vendor/`) and installed from a **pre-built arm64 wheel**.
    - *Cython won't compile v1.0.3 as-is:* `calignment.pyx` does `from calignment cimport AlignedPair` — a self-cimport modern Cython (0.29.36) rejects (the sibling `.pxd` is auto-applied). `pipelines/te_erv/vendor/build_telescope_wheel.sh` strips that line and builds the wheel with `--no-build-isolation` (so the conda numpy/cython are used). Rebuild it if the env's Python minor version changes.
15. **Nextflow null-param interpolation:** an unset process param (e.g. `params.bowtie2_extra_args`) interpolates into a task script as the literal string `"null"`, which a tool then treats as a positional argument (bowtie2 wrote its SAM to a file named `null`, leaving the samtools pipe empty). All optional `*_extra_args` in `te_erv.nf` are coalesced (`def x = params.x ?: ""`) before interpolation.
16. **Fetch source — AWS Open Data first, ENA fallback:** ENA HTTPS (`ftp.sra.ebi.ac.uk`) repeatedly stalls mid-stream on this network (`curl (56) OpenSSL SSL_read: unexpected eof`), turning a multi-hour cohort fetch into a retry crawl. `fetch_fastq.sh` now tries the **AWS Open Data SRA mirror** first: `curl` the `.sra` object from `https://sra-pub-run-odp.s3.amazonaws.com/sra/<ACC>/<ACC>` (anonymous, ~25–30 MB/s steady, no drops), then `fasterq-dump --split-files --threads 6` regenerates the paired FASTQs and they're gzipped to `<ACC>_1/_2.fastq.gz`. Measured ~5–8× faster end-to-end (a 4.3 GB `.sra` → paired FASTQ in ~3.5 min vs 20–40 min on ENA). If a run has no `.sra` on ODP, or fasterq-dump fails, it falls back to the original resilient ENA per-file curl loop. `FETCH_FORCE_ENA=1` skips AWS. **Verification asymmetry:** AWS-regenerated gzip does NOT reproduce ENA's `fastq_md5` (different compressor), so the AWS path verifies by presence + non-empty + size (trusting fasterq-dump's spot accounting); the ENA path still md5-checks against the manifest. Requires `sra-tools` (already in the `nextflow` env). **Transient disk cost:** fasterq-dump writes uncompressed FASTQ first (~18.9 GB/file, ~38 GB/sample) before gzip — the driver processes runs sequentially so the peak is one sample's worth.
17. **Host-mount instability on long background runs:** the granted mount `/Users/alex/OrchestratedBiosciences` intermittently drops (`Operation not permitted` on every path under it, while `/Users/alex` stays readable) during multi-hour background cells, killing the run partway. Fresh foreground shell invocations re-acquire the mount; in-kernel persistent processes keep the broken handle. This (not disk or a pipeline error) is what has repeatedly ended long cohort runs mid-QC — the log simply stops with no error banner and the `nextflow`/`java` process is gone. **Mitigation:** the batched runner + `-resume` makes every relaunch resumable from published BAMs, and the AWS fast-fetch means more samples complete per session-window before any drop. If a run dies, check `results/rnaseq_cohort/hisat2/*.bam` count — published BAMs are safe; just relaunch on a remaining-samples manifest.

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
bowtie2. UV/pip cannot help — STAR is a compiled C++ binary, not a Python package.

### 7a. Root cause FOUND — it is libc++, not arm64, and it is fixable (one line)

The "0 reads" bug is **not** arm64-specific, **not** a conda packaging defect, and
**not** an x86-SIMD/codegen problem. A fresh native arm64 recompile (Apple clang,
SIMDe→NEON, libomp) *still* ingests 0 reads — identical to the conda binary. The
actual cause is a **libc++ vs libstdc++ standard-library difference**:

STAR wires each read-chunk buffer to an `istringstream` via `pubsetbuf`
(`readInStream[ii]->rdbuf()->pubsetbuf(chunkIn[ii], …)` in `ReadAlignChunk.cpp`).
`std::basic_stringbuf::setbuf` (what `pubsetbuf` calls) is a **documented no-op in
libc++** (the only practical C++ stdlib on macOS), whereas **libstdc++ (Linux)
honors it**. Under libc++ the stream stays empty → first `peek()` returns EOF →
"0 input reads", exit 0. HISAT2/bowtie2 don't use this buffer trick. `Log.out`
fingerprint: `Thread #N end of input stream, nextChar=-1` before any read.

**The fix** (one line in `source/ReadAlignChunk_mapChunk.cpp`, inside
`ReadAlignChunk::mapChunk()`): load the chunk buffer via `.str()` (which libc++
honors) instead of relying on the no-op `pubsetbuf`:
```cpp
readInStream[ii]->str(string(chunkIn[ii], chunkInSizeBytesTotal[ii]));  // added, before clear()/seekg
```

A patched arm64 STAR 2.7.11b built this way **ingests reads correctly** —
independently verified in this repo (chr21 test: 15,056/15,056 input reads,
13,221 uniquely mapped, 37,978 BAM records, and a non-empty
`Chimeric.out.junction` — the input STAR-Fusion consumes).

**Installed binary:** `pipelines/bin/star-arm64/STAR-arm64-native` (self-contained;
bundles `lib/libomp.dylib` via `@loader_path/lib` rpath, code-signed). Build recipe
+ patch + STAR-Fusion setup notes: `pipelines/bin/star-arm64/BUILD-NOTE.md`.

**Spine decision unchanged:** the 40-sample spine stays on HISAT2 — it is done,
validated, and there is no reason to re-align. The patched STAR exists to **unblock
nf-core/rnafusion** (see §7b), which categorically needs STAR's chimeric junctions.

### 7b. rnafusion — now unblockable, still needs the CTAT resource library

With the patched arm64 STAR, the only remaining rnafusion requirement is the
**CTAT genome resource library** (`GRCh38_gencode_v??_CTAT_lib_*.plug-n-play.tar.gz`,
~30–40 GB) — it bundles the full-genome STAR index, fusion-annotation databases,
and coding-effect data. Use the prebuilt plug-n-play lib (building the full-genome
STAR index directly needs ~30 GB RAM). Point rnafusion at
`pipelines/bin/star-arm64/STAR-arm64-native` + a CTAT lib matching GENCODE. STAR-Fusion
v1.15.1 itself installs on system Perl (bundled `Set::IntervalTree` XS patched for
libc++/tr1; see BUILD-NOTE). This supersedes the earlier "rnafusion is arm64-blocked"
ruling — it is blocked only on the CTAT download, not on the toolchain.

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

### 8a. Disk hygiene during runs

Disk is at a premium on this machine. Two practices keep it bounded:

- **Prune Nextflow work dirs after a pipeline completes.** Each pipeline stages
  intermediates under `results/large/nf_work_<name>/`; once a run finishes (exit 0)
  and its results are published under `results/<name>/`, that work dir is disposable
  temp — `rm -rf results/large/nf_work_<name>`. These can be very large (a single
  botched te_erv run left a 274 GB orphaned `null` file from the pre-fix
  null-param bug — see §5). Keep the work dir only while you still intend to
  `-resume` (it holds the cached task outputs). Superseded reference builds are
  reclaimable too: the abandoned STAR pilot left a 44 GB duplicate STAR index under
  `results/rnaseq_pilot/genome/` (the live index is `reference/GRCh38/star_index`).
- **Archive spine BAMs to CRAM in the cohort loop.** `run_cohort_batched.sh` takes
  `ARCHIVE_CRAM=1` to convert each batch's spine BAMs to CRAM (via
  `archive_bam_to_cram.sh`) right after the FASTQ delete. Leave it UNSET while the
  custom subworkflows still need to read the BAMs (run them first), then re-invoke
  with `ARCHIVE_CRAM=1`, or archive by hand. `REF_FASTA` and `CRAM_THREADS`
  (default 6) are overridable via env.
