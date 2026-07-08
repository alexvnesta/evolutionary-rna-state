# RNA-editing subworkflow (A-to-I sites + Alu Editing Index)

Quantifies adenosine-to-inosine (A-to-I) RNA editing — one of the six RNA
phenotypes in the evolutionary-RNA-state project. It produces two
complementary readouts from the standard STAR genome BAMs:

1. **Per-site A-to-I calls** (JACUSA2 `call-1`) — a table of edited positions
   with editing frequency and depth, strand-aware (A>G on `+`, T>C on `-`),
   optionally masked against dbSNP.
2. **Alu Editing Index (AEI)** — the robust, cohort-comparable summary of
   *global* editing: the A>G mismatch rate pooled over **all** Alu elements
   genome-wide (Roth SH, Levanon EY, Eisenberg E, *Nat Methods* 2019). Because it
   aggregates millions of Alu adenosines it is stable at modest depth and is the
   recommended starting point for cross-sample / cohort comparison. The script
   also reports a non-A>G "noise floor" and a signal-to-noise ratio as QC.

## Why these tools on arm64

- **JACUSA2** is a **noarch Java jar** (bioconda) → runs natively on Apple
  Silicon under `-profile conda`, no Docker, no Rosetta. It is the primary
  per-site caller (single-BAM `call-1` mode detects A>G mismatches).
- **AEI** is computed by a small **pysam** script (`bin/compute_aei.py`); pysam,
  samtools, bedtools all have native `osx-arm64` conda builds.
- **REDItools** (pip pure-python) is the documented alternative per-site caller
  (known-sites mode against REDIportal); JACUSA2 is preferred here because it is
  a clean noarch jar with no extra dependencies.

## Inputs

| Input | Source | Notes |
|-------|--------|-------|
| Genome BAM (`[meta, bam, bai]`) | nf-core/rnaseq `hisat2` (arm64) / `star_salmon` (amd64) | coordinate-sorted; consumed, **not** re-aligned |
| `--fasta` | `reference/GRCh38/GRCh38.primary_assembly.genome.fa` | chr-prefixed (matches UCSC) |
| `.fai` | `samtools faidx` (via `fetch_alu.sh`) | required by JACUSA2 + pysam |
| `--rmsk` | UCSC hg38 `rmsk.txt.gz` | Alu source for the AEI |
| `--editing_snp_bed` (optional) | dbSNP / UCSC snp track, bgzipped+tabixed | mask known SNPs |

### Annotation sources (probed at authoring time)

- **UCSC RepeatMasker** `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz`
  — **confirmed reachable (HTTP 200)**. Columns parsed: `genoName, genoStart,
  genoEnd, strand, repName, repClass(=SINE), repFamily(=Alu)`. `bin/make_alu_bed.py`
  extracts the Alu-only BED6.
- **REDIportal** (known editing sites) `http://srv00.recas.ba.infn.it/atlas/`
  — **NOT reachable** from the build sandbox (not on the network allowlist).
  It is **optional** (AEI does not need it) and is therefore **parameterized**:
  download `TABLE1` (hg38) manually and pass `--editing_known_sites`. See
  `fetch_alu.sh` for the exact note.

Run `pipelines/rna_editing/fetch_alu.sh` once to download `rmsk.txt.gz` and build
the FASTA `.fai`.

## Outputs (`--outdir`, default `results/rna_editing/`)

```
rna_editing/
├── refs/alu.hg38.bed6              # Alu intervals (built once)
├── sites/<sample>.editing_sites.tsv   # per-site A-to-I calls (JACUSA2, filtered)
├── sites/<sample>.jacusa2.out         # raw JACUSA2 call-1 output
└── aei/
    ├── <sample>.aei.tsv            # per-sample AEI + mismatch-type breakdown + S/N
    └── cohort_aei.tsv             # all samples merged (the cohort-level summary)
```

`cohort_aei.tsv` columns: `sample, AEI_percent, AG_mismatches, A_coverage,
signal_to_noise, noise_floor_percent, cov_{A,C,G,T}, n_{12 mismatch types}`.

## Run command

```bash
source pipelines/env.sh
pipelines/rna_editing/fetch_alu.sh          # one-time: rmsk.txt.gz + genome .fai

nextflow run pipelines/rna_editing/main.nf \
  -profile conda \
  -c pipelines/conf/mac_arm64.config \
  -c pipelines/rna_editing/conf/editing.config \
  --bam_glob 'results/rnaseq/hisat2/*.markdup.sorted.bam' \
  --fasta reference/GRCh38/GRCh38.primary_assembly.genome.fa \
  --rmsk  reference/GRCh38/repeats/rmsk.hg38.txt.gz \
  --outdir results/rna_editing \
  -work-dir results/large/nf_work_rna_editing \
  -resume
```

Samplesheet mode instead of a glob: `--input editing_samplesheet.csv` with
columns `sample,bam,bai`.

Disable the (heavier) per-site caller and compute only the AEI:
`--editing_call_sites false`. Add dbSNP masking with `--editing_snp_bed dbsnp.bed.gz`.

Help / config parse check: `nextflow run pipelines/rna_editing/main.nf --help`.

## Key parameters (defaults)

| Param | Default | Meaning |
|-------|---------|---------|
| `--editing_call_sites` | `true` | run JACUSA2 per-site caller |
| `--editing_min_mapq` | `60` | HISAT2 unique-mapper MAPQ (arm64); STAR uses 255 |
| `--editing_min_baseq` | `25` | min base quality |
| `--editing_min_cov` | `10` | min coverage to call a site |
| `--editing_min_freq` | `0.10` | min editing frequency |
| `--editing_min_edit_reads` | `3` | min edited reads |
| `--editing_snp_bed` | `null` | dbSNP BED to mask (recommended for site calls) |
| `--editing_alu_standard_only` | `true` | Alu on chr1..22,X,Y,M only |

## arm64 caveats

- All tools resolve on `osx-arm64` via `[conda-forge, bioconda]` (see
  `environment.yml`): JACUSA2 (noarch jar), samtools/htslib/pysam/bedtools
  (native arm64). No Docker/Singularity.
- **dbSNP masking is strongly recommended** for the per-site calls: A>G
  mismatches at germline A/G SNPs are not editing. The AEI is comparatively
  robust to this (it pools genome-wide and reports a noise floor), but masking
  still improves it. dbSNP is a parameter because its size/version is a project
  choice.
- Resource directives (`conf/editing.config`): editing steps are light
  (≤10 GB); JACUSA2 is `maxForks 1`. This never triggers the 40 GB STAR path —
  BAMs are consumed, not re-aligned.
