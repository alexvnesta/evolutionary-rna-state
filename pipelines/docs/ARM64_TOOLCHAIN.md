# Apple Silicon (osx-arm64) toolchain map

Determines which tools run **natively** under Nextflow `-profile conda` (no
Docker, no Rosetta) vs. which need a source build or pip install, on this
M-series Mac (64 GB, arm64).

**Provenance of each row:**
- Rows marked *(API-probed)* were checked directly against the anaconda.org
  package API for an `osx-arm64` build: STAR, salmon, samtools, telescope-ngs,
  tetranscripts, reditools/reditools2, sprint, jacusa2, arriba, star-fusion,
  stringtie, subread, hisat2, rseqc, rmats/rmats2sashimiplot, sra-tools. `fq`
  (nf-core/rnaseq QC module) was later API-probed and confirmed **not** on
  osx-arm64 (linux-64/linux-aarch64/osx-64 only).
- Rows for **IRFinder-S**, **FusionCatcher**, and **SUPPA2** are **design
  decisions, not API probes** — IRFinder-S and FusionCatcher were not queried
  against the anaconda.org API; their handling reflects known packaging (no
  osx-arm64 bioconda build as of writing) and the project's deferral choice.
  SUPPA2 is pip/noarch (installable anywhere). Re-probe before relying on these.

## Native osx-arm64 conda builds (bioconda) — run natively
| Tool | Phenotype | Notes |
|------|-----------|-------|
| STAR | base align, TE, IR-input, editing-input | genome index build ~35 GB RAM |
| Salmon | base quant, SUPPA2 | already installed in `rnaio` (arm64) |
| samtools | all | already installed in `rnaio` (arm64) |
| fastp | trimming | already installed in `rnaio` (arm64) |
| Arriba | fusion | native arm64 — primary fusion caller |
| rMATS | splicing | native arm64 — primary splicing caller |
| StringTie | splicing/ORF assembly | native arm64 |
| subread (featureCounts) | TE/gene counts | native arm64 |
| HISAT2 | alt aligner | native arm64 |
| sra-tools | data acquisition | native arm64 |

## noarch (pure Python / Java) — run anywhere incl. arm64
| Tool | Phenotype | Notes |
|------|-----------|-------|
| TEtranscripts / TElocal | TE/ERV | noarch; primary TE quant |
| JACUSA2 | RNA editing | noarch Java; primary editing caller |
| RSeQC | QC | noarch |
| SUPPA2 | splicing | pip/noarch |

## pip-installable pure Python — arm64 fine
| Tool | Phenotype | Notes |
|------|-----------|-------|
| Telescope | TE/ERV (locus-level, EM reassignment) | `pip install telescope-ngs` |
| REDItools (2/3) | RNA editing | github/pip; python+samtools |

## Needs source build or deferral on arm64
| Tool | Phenotype | Plan |
|------|-----------|------|
| IRFinder-S | intron retention | RESOLVED: pipelines/intron_retention/ ships a featureCounts-based IR ratio (native arm64, primary). IRFinder-S source build documented as optional/deferred higher-fidelity path. |
| STAR-Fusion | fusion | perl+noarch; usable, but CTAT ref ~30 GB. Arriba is the native primary; STAR-Fusion secondary |
| FusionCatcher | fusion | heavy, amd64-oriented — DEFER to future cluster compute |

## Consequence
`-profile conda` is the correct execution mode on this Mac. Docker/Singularity
are NOT required and NOT installed. Core align+quant+splicing+fusion+TE+editing
all have native or noarch arm64 paths. Only IRFinder needs a source build and
FusionCatcher is deferred.

## nf-core/rnaseq 3.26.0 module-level arm64 gaps (discovered during pilot)

Auditing all 74 module conda specs against the anaconda.org API found 10 exact
version-pins with no osx-arm64 build. Handling (see `conf/arm64_module_overrides.config`):

| Module pin | osx-arm64 status | Fix |
|------------|------------------|-----|
| `perl=5.26.2` (gtf2bed) | only 5.32.x on arm64 | override → perl=5.32.1 |
| `subread=2.0.6` (featurecounts) | 2.0.8/2.1.1 | override → subread=2.0.8 |
| `stringtie=2.2.1` (merge) | 2.2.3 / 3.x | override → stringtie=2.2.3 |
| `ucsc-bedclip=377` | 482 | override → 482 |
| `ucsc-bedgraphtobigwig=469` | 482 | override → 482 |
| `sed=4.7` (gencode preprocess) | 4.8/4.9/4.10 | override → sed=4.8 |
| `fq=0.12.0` (lint/subsample) | only 0.15+ (incompatible) | **skip** via `--skip_linting` + explicit strandedness |
| `future=0.18.3`, `pysam=0.22.0` (umitools/rsem) | near versions exist | only used with `--with_umi`/`star_rsem`; not on our star_salmon path |

**RSeQC runs natively — do NOT skip it.** `rseqc=5.0.4` is a **noarch** conda
package (architecture-independent) and its `r-base=4.3` dependency resolves on
osx-arm64 (`micromamba create rseqc=5.0.4 r-base=4.3` solves). An earlier audit
of this file incorrectly flagged it as arm64-unavailable because the check only
looked at the `osx-arm64` package subdir and ignored `noarch`; corrected here.

Net: the STAR+Salmon spine runs natively on arm64 with 6 version-overrides and
**one** QC-step skip (fq linting only). Alignment, quantification, trimming,
dedup, bigWig, featureCounts, StringTie, RSeQC, and MultiQC all run.
