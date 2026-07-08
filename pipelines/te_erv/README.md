# TE / ERV activation subworkflow

Quantifies transposable-element (TE) and endogenous-retrovirus (ERV)
expression for the evolutionary-RNA-state melanoma/ICB study. TE/ERV
reactivation is one of the six RNA phenotypes and is central to the project
hypothesis (TE-derived transcripts as a source of tumour antigenicity). No
mature nf-core pipeline exists for this, so it is authored here as a
self-contained DSL2 subworkflow.

## What it quantifies

Two complementary resolutions, from two alignment inputs:

| Level  | Tool          | Input BAM                                   | Output |
|--------|---------------|---------------------------------------------|--------|
| **Locus** | **Telescope** (EM reassignment) | dedicated multimap-aware STAR run | per-locus HERV/L1 counts |
| **Family** | **TEtranscripts / TEcount** | the coordinate-sorted genome BAM from the nf-core/rnaseq `hisat2` spine (HISAT2 aligner on arm64; STAR on amd64) | per-family/subfamily counts (Alu, L1, ERVK, …) |

**Why two alignments.** Telescope resolves *which specific TE/ERV locus* a
multimapped read came from via an EM model, and for that it needs an alignment
that *retains* multimappers. The standard rnaseq spine BAM discards them, so
the locus branch runs its own `STAR_ALIGN_MULTI` with
`--outFilterMultimapNmax 100 --winAnchorMultimapNmax 200`, unsorted output, and
`--outSAMprimaryFlag AllBestScore` (all alignments of a read kept together — the
layout Telescope expects). TEcount, in `--mode multi`, distributes multimappers
statistically at the *family* level and works directly from the existing sorted
spine BAM, so that branch re-uses it and adds no alignment cost.

## Processes

```
                 ┌─ STAR_ALIGN_MULTI (FASTQ) ─→ TELESCOPE_ASSIGN ─→ MERGE_TELESCOPE ─→ telescope_locus_counts.tsv
  samplesheet ──┤
                 └─ (rnaseq spine BAM) ───────→ TETRANSCRIPTS_COUNT ─→ MERGE_TETRANSCRIPTS ─→ tetranscripts_family_counts.tsv
```

- `STAR_ALIGN_MULTI` — 12 cpu / 40 GB / `maxForks 1` (the memory hog; serialized).
- `TELESCOPE_ASSIGN`, `TETRANSCRIPTS_COUNT` — 4 cpu / 10 GB.
- `MERGE_*` — 2 cpu / 8 GB; `bin/merge_matrices.py` builds the union count matrices.

## Inputs

Samplesheet CSV (`assets/samplesheet_test.csv` is a template):

```
sample,fastq_1,fastq_2,bam,strandedness
pilot01,…/pilot01_R1.fastq.gz,…/pilot01_R2.fastq.gz,…/pilot01.markdup.sorted.bam,reverse
```

- `fastq_1`/`fastq_2` feed the locus (Telescope) branch. Leave `fastq_2` empty for single-end.
- `bam` is the rnaseq-spine coordinate-sorted genome BAM; feeds the family (TEcount) branch.
- A row may populate either or both branches; empty columns are simply skipped.

Reference parameters:

| Param | Meaning |
|-------|---------|
| `--star_index`    | prebuilt shared STAR index (`reference/GRCh38/star_index`) |
| `--gene_gtf`      | GENCODE genic GTF (`reference/GRCh38/gencode.v46.primary_assembly.annotation.gtf`) |
| `--te_gtf_locus`  | Telescope locus GTF (`reference/te/retro.hg38.v1.transcripts.gtf`) |
| `--te_gtf_family` | TEtranscripts family GTF (see **Annotation** below) |

## Outputs

```
<outdir>/
  star_multi/<id>.Aligned.out.bam        multimap-aware BAM (per sample)
  telescope/<id>-telescope_report.tsv    per-sample locus report
  tetranscripts/<id>.cntTable            per-sample family counts
  matrices/telescope_locus_counts.tsv    LOCUS-level count matrix (samples × loci)
  matrices/tetranscripts_family_counts.tsv  FAMILY-level count matrix
```

## Annotation

Run `bin/fetch_te_annotation.sh` to stage annotations under `reference/te/`.

- **Locus (confirmed, auto-fetched):** `retro.hg38.v1` (HERV_rmsk + L1Base) from
  the `mlbendall/telescope_annotation_db` repo (Git LFS →
  `media.githubusercontent.com`). Verified reachable: HTTP 200, 18.9 MB,
  **28,513 distinct HERV/L1 loci** across 72,169 exon-feature lines, UCSC `chr`
  naming (matches the GENCODE primary assembly).
- **Family (parameterized):** the canonical `GRCh38_GENCODE_rmsk_TE.gtf` is
  distributed from the Hammell lab (`mghlab.org` / `labshare.cshl.edu`); that
  download path did **not** resolve through this environment at authoring time,
  so `--te_gtf_family` is a required parameter. If you have the file, drop it at
  `reference/te/GRCh38_GENCODE_rmsk_TE.gtf`. Otherwise re-run the fetch script
  with `BUILD_FAMILY=1` to build a reproducible fallback GTF from the
  allowlisted UCSC hg38 RepeatMasker table (`hgdownload.soe.ucsc.edu`).

## Run command

```bash
source pipelines/env.sh          # Nextflow 26.04.4 + JDK17 + STAR 2.7.11b (arm64)

nextflow run pipelines/te_erv/main.nf \
  -profile conda \
  -c pipelines/conf/mac_arm64.config \
  --input        pipelines/te_erv/assets/samplesheet_test.csv \
  --outdir       results/te_erv \
  --star_index   reference/GRCh38/star_index \
  --gene_gtf     reference/GRCh38/gencode.v46.primary_assembly.annotation.gtf \
  --te_gtf_locus  reference/te/retro.hg38.v1.transcripts.gtf \
  --te_gtf_family reference/te/GRCh38_GENCODE_rmsk_TE.gtf \
  -work-dir      results/large/nf_work_te_erv \
  -resume
```

Dry validation (no alignment; proves the DSL2 parses and the DAG wires):

```bash
nextflow run pipelines/te_erv/main.nf -c pipelines/conf/mac_arm64.config \
  -stub-run -process.conda=false --input <stub.csv> --star_index <dir> \
  --gene_gtf <gtf> --te_gtf_locus <gtf> --te_gtf_family <gtf>
```

## arm64 caveats

- **Native / noarch (run directly):** STAR 2.7.11b, samtools, subread, and
  **TEtranscripts 2.2.4** (noarch) all have osx-arm64/noarch conda builds.
- **Telescope has no arm64 conda build.** bioconda `telescope` is osx-64/linux-64
  only, and the PyPI `telescope-ngs` wheel is x86/py36. The `environment.yml`
  therefore installs it via `pip` from source
  (`git+https://github.com/mlbendall/telescope.git@main`) — pure Python plus one
  Cython extension — compiled against the conda-provided
  numpy/scipy/pysam/cython/c-compiler (all arm64). First env creation compiles
  that extension, so allow extra time (`conda.createTimeout = '2 h'` is already
  set in the machine config).
- No Docker/Singularity; execution is `-profile conda` only, per the repo
  toolchain map.

## Files

- `te_erv.nf` — module: processes + `TE_ERV` subworkflow.
- `main.nf` — standalone entrypoint + samplesheet parsing + `--help`.
- `nextflow.config` — binds the conda env and per-process resources.
- `environment.yml` — arm64-installable tool env.
- `bin/merge_matrices.py` — per-sample → count-matrix merge.
- `bin/fetch_te_annotation.sh` — stage locus GTF (+ optional family fallback).
- `assets/samplesheet_test.csv` — template samplesheet.
