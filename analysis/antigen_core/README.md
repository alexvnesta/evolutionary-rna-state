# antigen_core — shared antigen-derivation keystone

The keystone every downstream antigen module (splicing, TE/ERV, intron
retention, fusion, SNV/indel) calls. It provides two things and nothing else:

1. **`mhc_binding.py`** — the MHC-I peptide→binder engine (MHCflurry 2.x).
2. **`hla_typing.py`** — RNA HLA-I genotyping (arcasHLA) + heterozygosity flag.

Siblings never touch MHCflurry or arcasHLA directly — they import from here so
every `*_neoantigen_burden` feature is defined identically and comparably.

## MHC binding engine

```python
from analysis.antigen_core.mhc_binding import (
    score_peptides,   # (peptides, hla_alleles) -> tidy per-peptide DataFrame
    count_binders,    # (peptides, hla_alleles, rank_threshold=2.0) -> int
    binder_counts,    # (peptides, hla_alleles) -> {n_strong_binders, n_weak_binders, n_scored}
    best_per_peptide, # (peptides, hla_alleles) -> one row per peptide (best allele)
    STRONG_BINDER_RANK,  # 0.5  (affinity percentile)
    WEAK_BINDER_RANK,    # 2.0
)

# a splicing module, per sample:
peps    = derive_peptides_from_splice_junctions(sample)   # 8-11mers
alleles = ["HLA-A*02:01", "HLA-A*01:01", "HLA-B*07:02", ...]  # from hla_typing
burden  = count_binders(peps, alleles)      # -> splice_neoantigen_burden
```

- **Binder definition (field standard, exposed as constants):** strong =
  affinity percentile ≤ 0.5, weak = ≤ 2.0 (weak is inclusive of strong).
  `count_binders` defaults to the weak threshold; pass `STRONG_BINDER_RANK`
  for strong-only.
- **Unit = unique peptide.** A peptide presented by several of the sample's
  alleles is one antigen. `count_binders` counts unique peptides.
- **Input hygiene handled for you:** peptides are uppercased, deduped, and
  filtered to standard-AA 8–11mers; alleles are normalized and unsupported
  ones dropped. Empty/invalid input → 0 / empty DataFrame, never an error.
- **Percentile rank is the batch-robust unit** (calibrated per allele against
  a fixed random-peptide background), so the burden counts are comparable
  across cohorts/platforms — the project's hard reproducibility constraint.

### Models (one-time fetch)

MHCflurry models are re-fetchable from GitHub releases (openvax). They must
land in a writable dir — the macOS default (`~/Library/Application Support`)
is not writable in the sandbox, so we fetch into the repo:

```bash
export MHCFLURRY_DATA_DIR=$REPO/reference/mhcflurry_models
mhcflurry-downloads fetch models_class1_presentation
```

`mhc_binding.py` auto-points `MHCFLURRY_DATA_DIR` at `reference/mhcflurry_models`
if it exists and the env var is unset, so imports "just work".

## HLA typing

```python
from analysis.antigen_core.hla_typing import type_sample, build_hla_table
row = type_sample(bam, run_accession="ERR2208952", cohort="gide2019",
                  outdir="results/hla/ERR2208952")
# -> {run_accession, cohort, HLA_A_1..HLA_C_2, HLA_I_heterozygous, n_het_loci, tool, tool_version}
```

- **Tool: arcasHLA** (RabadanLab/arcasHLA), the standard RNA-seq HLA-I typer.
  It genotypes directly from a STAR BAM (extract chr6 + unmapped → kallisto
  pseudoalign vs IPD-IMGT/HLA → genotype). Installs on osx-arm64 via bioconda.
  **OptiType** is the documented fallback (`parse_optitype_result` converges on
  the same schema).
- **`HLA_I_heterozygous`** = heterozygous at A **and** B **and** C (2-field
  resolution) — the Chowell 2017 favorable checkpoint-response feature.
- **Nextflow:** `hla_typing.nf` is the subworkflow stub consuming the pipeline
  session's `(meta, bam, bai)` BAM channel; `bin/merge_hla_table.py` merges
  per-sample genotype JSONs into the cohort tidy table.

### arm64 caveat (why full typing runs on the pilot host)

arcasHLA + kallisto build and launch on the arm64 dev box, but arcasHLA's
prebuilt IMGT/HLA kallisto index expects kallisto 0.44.0 while bioconda ships
0.52.0 on arm64 — the index must be rebuilt (`arcasHLA reference`) against the
installed kallisto before genotyping. Full typing therefore runs on the Linux
pipeline host with a matched kallisto. The **parsing + heterozygosity logic**
(the part that must be exactly right) is unit-tested here on synthetic
genotypes; **no real-sample allele calls are fabricated**.

## Tests

```bash
cd analysis/antigen_core
PYTHONPATH="$PWD" python test_mhc_binding.py   # known A*02:01 epitopes vs decoys
PYTHONPATH="$PWD" python test_hla_typing.py    # parse + heterozygosity logic
```
