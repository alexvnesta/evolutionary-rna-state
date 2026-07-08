# Differentiated bucket — retained-intron antigen module

`intron_retention.py` derives two **named, interpretable** per-sample features
from the pipeline session's intron-retention quantification, keyed on
`run_accession` + `cohort` (FEATURE_CONTRACT_v2):

| feature | dtype | meaning |
|---|---|---|
| `retained_intron_load` | int | # introns with IR ratio ≥ threshold (0.10) in the sample |
| `ir_neoantigen_burden` | int | MHC-I binder count over peptides translated from retained introns, via the **shared** `antigen_core` engine |

Companion / QC columns (batch-robust forms): `n_introns_evaluated`,
`retained_intron_fraction` (= load / evaluated; depth-robust),
`retained_intron_load_weighted` (Σ IR of retained introns),
`retained_intron_load_cohortz` (within-cohort robust z of the fraction),
`n_candidate_peptides`, `n_retained_introns_used`.

## Cryptic-ORF peptide derivation

Each retained intron (IR ≥ threshold, top-`MAX_INTRONS_FOR_BURDEN` by IR) is
translated into two candidate-peptide classes, tiled as 8–11mers:

1. **Junction-spanning read-through peptides** — the 5′ exon flank + intron are
   translated in all 3 frames; kmers that straddle the exon→intron boundary
   (the novel read-through sequence, pre-premature-stop) are kept.
2. **Intronic ORF peptides** — internal ATG…stop ORFs wholly inside the intron,
   all 3 frames.

Strand-aware (`-` strand introns are reverse-complemented; the exon flank is
taken on the transcriptional 5′ side). Pooled unique peptides + the sample's 6
HLA-I alleles → `antigen_core.mhc_binding.count_binders` → `ir_neoantigen_burden`.
Because every antigen module calls the same engine, the burden is directly
comparable to `splice_/te_/fusion_/snv_indel_neoantigen_burden`.

## Batch robustness (IR is coverage/depth sensitive)

The IR **ratio** is a within-sample intron/exon read-density ratio
(self-normalising for library size). The derived **count** of retained introns
is depth-sensitive, so the module: counts only introns evaluable in the sample
(non-NA IR = host gene had sufficient exonic coverage upstream); emits a
depth-robust `retained_intron_fraction` and within-cohort `_cohortz`; caps
candidate introns for the burden; and calls binders on MHCflurry percentile
rank (allele-calibrated, batch-invariant). **Report per clinical context /
z-scored within cohort; never pool raw counts across platforms.** Full text in
`intron_retention.BATCH_ROBUSTNESS_NOTE`.

## Files

- `intron_retention.py` — feature logic (importable module).
- `test_intron_retention.py` — 5 synthetic-input unit tests (load logic,
  seq helpers, cryptic-ORF peptide gen with a planted A\*02:01 epitope,
  burden through the shared engine, end-to-end). No real-cohort values fabricated.
- `bin/ir_antigen_features.py` — CLI wrapper (pipeline outputs → feature parquet).
- `intron_retention.nf` + `modules/ir_antigen_features.nf` — Nextflow subworkflow
  stub (`IR_NEOANTIGEN`) consuming `intron_retention.parquet` + `introns.saf` +
  genome FASTA + `hla_typing.parquet`. Validated with `nextflow -stub-run`.

## Run

```bash
# unit tests (needs the 'antigen' env: mhcflurry + models, pandas, pysam)
cd analysis/differentiated
PYTHONPATH="../antigen_core:." python test_intron_retention.py

# CLI on real pipeline outputs (pilot)
python bin/ir_antigen_features.py \
    --ir-matrix   results/features/intron_retention.parquet \
    --intron-saf  <make_intron_saf output>/introns.saf \
    --genome      reference/GRCh38/GRCh38.primary_assembly.genome.fa \
    --hla-table   results/features/hla_typing.parquet \
    --out         results/features/ir_antigen_features.parquet
```
