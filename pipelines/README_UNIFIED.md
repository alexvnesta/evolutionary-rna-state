# Unified non-reference RNA feature workflow (Apple Silicon)

`pipelines/main.nf` — one top-level DSL2 workflow that takes a spine BAM channel and fans it out to
the validated per-feature subworkflows, each emitting a cohort matrix. This is the single
"handles-everything" entry the project's north star calls for. It is a **wrapper over arms already
validated on osx-arm64**, not a reimplementation.

## Design
```
   spine BAMs (HISAT2 / nf-core rnaseq)   results/rnaseq_cohort/hisat2/*.markdup.sorted.bam
              │  one Channel: [ meta, bam, bai ]
     ┌────────┼─────────────────────────┐
     ▼        ▼                          ▼
 RNA_EDITING  INTRON_RETENTION        TE_ERV
 (AEI,        (featureCounts          (TEtranscripts family-level;
  unique-     intron/exon →            +optional Telescope locus-level
  mapper)     IR ratio)                when --te_locus)
     │        │                          │
     ▼        ▼                          ▼
 cohort_aei  intron_retention matrix   te family/locus matrix
```
Splicing (`rnasplice`) and fusion (`rnafusion`) are wrapped nf-core pipelines that re-align from
FASTQ rather than consuming the spine BAM, so they run as separate entries (`pipelines/scripts/`)
and their per-sample outputs are merged into the final matrix by `analysis/build_nonref_matrix.py`.

## Run
```bash
source pipelines/env.sh                       # arm64 conda/mamba shims + NXF_LOCAL_CHANNEL
nextflow run pipelines/main.nf -profile apple_silicon \
  --bam_glob 'results/rnaseq_cohort/hisat2/*.markdup.sorted.bam' \
  --outdir results/nonref_unified
```

## Options
| flag | default | meaning |
|---|---|---|
| `--bam_glob` / `--input` | — | spine BAM glob, or a `sample,bam[,bai]` samplesheet |
| `--te_locus` | `false` | family-level TE only (~1 min/sample). `true` adds Telescope locus (bowtie2 -k100, ~2 h/sample — use Modal/cloud) |
| `--editing_call_sites` | `false` | AEI only. `true` adds JACUSA2 per-site A-to-I calling (heavier) |
| `--outdir` | `results/nonref_unified` | output root |

## Profile
`-profile apple_silicon` (alias `arm64`) layers on `conf/mac_arm64.config`: conda arm64 envs,
micromamba solver, the local-channel repackaging that works around sandbox git-protection, and the
per-tool arm64 overrides (STAR skipped in favor of HISAT2, Qualimap/fq disabled, etc.). Verified to
compile and build a valid DAG on Nextflow 26.04 against the 32-sample cohort BAM set.

## Status
- **Compiles clean on Nextflow 26.04** under `-profile apple_silicon`, resolves every reference path,
  and the DSL for all three fan-out subworkflow calls parses with no errors (`-preview`). Note: NF 26's
  `-preview` does not emit process nodes to the DAG file, so this is a clean-compile verification, not a
  rendered fan-out graph — an actual `-resume` run on the 32 BAMs is the remaining end-to-end check.
- Each arm was independently validated end-to-end on real arm64 data in its own directory prior to
  unification; this wrapper composes them off a single BAM channel.
- **`--te_locus` (Telescope locus-level) is hard-gated OFF from this BAM-only entry** — it needs a FASTQ
  read source for bowtie2 -k100 re-alignment, which this entry does not consume. Run the standalone
  `te_erv/` pipeline with `--input <fastqs>` (or the cloud runner) for locus-level resolution. Family-level
  TE is what runs here.
