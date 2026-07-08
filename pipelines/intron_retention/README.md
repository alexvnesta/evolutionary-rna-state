# INTRON_RETENTION subworkflow

Per-intron, per-sample **intron-retention (IR) ratio** quantification for the
evolutionary-RNA-state pipeline. One of the six RNA phenotypes; **distinct from
splicing** (rnasplice / rMATS / SUPPA2 handle exon-skipping and PSI — this
subworkflow measures reads that remain *inside* introns).

Fully **arm64-native** — featureCounts (subread, native osx-arm64) + pandas.
No Docker, no source build, no Rosetta.

---

## What it quantifies

For every annotated intron, an IRFinder-like coverage ratio:

```
intron_density = intron_read_count / pure_intron_length
exon_density   = host_gene_exon_count / host_gene_exon_length
IR_ratio       = intron_density / (intron_density + exon_density)      # in [0, 1]
```

`IR_ratio → 0` = intron fully spliced out; `→ 1` = intron retained at (or above)
the depth of its host gene's exons.

**"Pure intronic" definition (why this is not naive gene-minus-exon):** intron
intervals are the gaps between a gene's union-exon blocks, then **masked against
the genome-wide union of all exons of every gene on both strands**. Any part of
a candidate intron that is exonic in some *other* (e.g. antisense or nested)
transcript is removed, so intronic counts are not contaminated by an overlapping
gene's spliced exons. Introns can be fragmented by this masking; fragments are
summed back into a single per-intron meta-feature (`GeneID = <gene_id>__intron_<n>`,
strand-aware numbering). See `bin/make_intron_saf.py`.

---

## Inputs

- **BAMs**: coordinate-sorted STAR **genome** BAMs from the nf-core/rnaseq
  `hisat2` spine (HISAT2 aligner on arm64; STAR on amd64) (`[meta, bam, bai]`). It consumes existing BAMs — it does
  **not** re-align.
- **GTF**: GENCODE annotation (`reference/GRCh38/gencode.v46.primary_assembly.annotation.gtf`).

Supply BAMs either via a samplesheet or a glob:

```
# samplesheet CSV (recommended) — columns:
run_accession,cohort,bam,bai,single_end
ERR2208952,gide2019,/abs/path/ERR2208952.markdup.sorted.bam,/abs/path/...bam.bai,false
```

`run_accession` is the hand-off join key (see `docs/HANDOFF_CONTRACT.md`).

---

## Outputs

Published under `--outdir` (default `results/intron_retention/`):

| File | Content |
|------|---------|
| `features/intron_retention.parquet` | **contract deliverable** — tidy-wide: rows = `run_accession`, cols `[run_accession, cohort, <intron_id…>]`, values = IR ratio (float32), NA where not evaluated |
| `features/intron_retention_summary.tsv` | per-sample summary: `n_introns_evaluated`, `median_IR`, `mean_IR`, `n_IR_gt_<thr>` |
| `<id>.ir_long.tsv` | per-sample long table (counts, lengths, IR ratio per intron) |
| `<id>.ir_ratio.parquet` | per-sample wide matrix |
| `<id>.ir_summary.tsv` | per-sample summary row |
| `reference/introns.saf`, `exons.saf`, `intron2gene.tsv` | derived intervals (built once, reusable) |
| `counts/*.featureCounts.txt` | raw featureCounts output (intron & exon passes) |

The cohort `intron_retention.parquet` matches `docs/HANDOFF_CONTRACT.md`
exactly (per-sample IR ratio per intron, keyed on `run_accession`).

---

## Run command

```bash
source pipelines/env.sh          # Nextflow 26.04.4 + JDK17 on PATH

nextflow run pipelines/intron_retention/main.nf \
    -profile conda \
    -c pipelines/conf/mac_arm64.config \
    -c pipelines/intron_retention/nextflow.config \
    --input   pipelines/intron_retention/samplesheet.csv \
    --gtf     reference/GRCh38/gencode.v46.primary_assembly.annotation.gtf \
    --outdir  results/intron_retention \
    -work-dir results/large/nf_work_ir \
    -resume
```

Key parameters (defaults in `nextflow.config`):

| Param | Default | Meaning |
|-------|---------|---------|
| `--ir_strandedness` | `2` | featureCounts `-s`: 0 unstranded, 1 fwd, 2 rev. **rev (2)** for Illumina dUTP libraries (nf-core/rnaseq default). Set to match your library. |
| `--ir_min_intron_len` | `50` | drop pure-intronic sub-intervals shorter than this (bp) |
| `--ir_min_gene_exon_count` | `20` | require this many host-gene exonic reads before an IR value is trusted (else NA) |
| `--ir_high_threshold` | `0.1` | IR ratio above which an intron is called "retained" in the summary |

`--input` accepts either the samplesheet above or `--bam_glob '/path/*.bam'`.

---

## Validation done

- **DSL2 parses**: `nextflow run main.nf --help` and a full **`-stub-run`**
  execute the whole DAG (`MAKE_INTRON_SAF → FEATURECOUNTS_IR[intron] +
  FEATURECOUNTS_EXON[exon] → COMPUTE_IR_RATIO → MERGE_IR_MATRIX`) across two
  samples, publishing every expected output. The exon pass is a second aliased
  instance of the featureCounts process (`FEATURECOUNTS_IR as FEATURECOUNTS_EXON`).
- **Conda solves on osx-arm64**: `environment.yml` resolves via micromamba
  `--platform osx-arm64` — subread 2.1.1 (native arm64), pandas, pyarrow,
  samtools, all with arm64/noarch builds.
- **Script logic unit-tested**: intron derivation (incl. antisense-exon masking
  and fragment summing) and IR-ratio math verified against hand-computed values;
  intron/exon featureCounts columns are paired **by sample name**, so the two
  passes align even if column order differs.

**Not yet run on real BAMs** — deferred until the rnaseq spine emits STAR genome
BAMs for the pilot cohort. Everything upstream of that (authoring + config parse
+ arm64 env solve + script logic) is validated.

---

## arm64 caveat & the higher-fidelity option

`bioconda::irfinder` (IRFinder-S) had **no osx-arm64 build** in our probe, so
this subworkflow uses the **featureCounts-based IR ratio** as the primary,
arm64-native deliverable — no source build, unblocks the Mac immediately.

**IRFinder-S is the higher-fidelity option** (it models intron-depth
distribution, splice-site overhang, and known-exon overlap more finely). To use
it, build from source for arm64 (deferred / optional):

```bash
# OPTIONAL — higher-fidelity path, requires a source build (deferred)
git clone https://github.com/RitchieLabIGH/IRFinder
cd IRFinder && ./install.sh          # C++/STAR-based; needs a working arm64 STAR
# IRFinder builds its own reference from GRCh38 + GTF, then emits per-intron
# IR ratios per sample (IRFinder-S 'IRratio' column) — same per-sample contract.
```

The featureCounts approach captures the same signal (intronic vs exonic read
density per intron, per sample) with a fully native toolchain; swap in IRFinder-S
later for cross-validation without changing the hand-off format.

## Files

```
pipelines/intron_retention/
├── main.nf                     # standalone entrypoint (samplesheet/glob → subworkflow)
├── subworkflow.nf              # INTRON_RETENTION workflow (include into a larger pipeline)
├── nextflow.config             # process labels + IR param defaults
├── environment.yml             # arm64-native conda env (subread, pandas, pyarrow, samtools)
├── modules/
│   ├── make_intron_saf.nf      # derive pure-intronic + exonic SAF from GTF
│   ├── featurecounts_ir.nf     # featureCounts (aliased for intron & exon passes)
│   ├── compute_ir_ratio.nf     # per-sample IR ratio
│   └── merge_ir_matrix.nf      # cohort matrix (contract format)
├── bin/                        # auto-added to PATH by Nextflow
│   ├── make_intron_saf.py      # pure-stdlib GTF → intron/exon intervals
│   ├── compute_ir_ratio.py     # IR ratio + per-sample summary
│   └── merge_ir_matrix.py      # pivot per-sample → cohort wide matrix
└── test/                       # stub-run samplesheet + mini GTF + dummy BAMs
```
