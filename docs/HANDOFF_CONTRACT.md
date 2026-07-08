# Feature hand-off contract — pipeline ↔ modeling

**Purpose.** Two sessions are working this project in parallel and coordinate
through this file + the artifact store (they cannot message each other):

| Session | Role | Owns |
|---------|------|------|
| **Find/Create Nextflow RNA-seq Pipelines** (`c15a540e`) | **feature factory** | raw reads → de-novo phenotype matrices (align/quant, splicing, IR, editing, TE, fusion) |
| **Built with Claude … Hackathon** (this one, `f16431e3`) | **modeling brain** | co-variation, latent state **S**, response organization, figures, write-up, rigor harness |

Neither session rebuilds the other's layer. The pipeline emits **per-sample
feature matrices** in the format below; the modeling session consumes them via
`src/features.py` and never touches STAR/nf-core itself.

---

## Why this contract exists (the result that drives it)

The modeling session tested the internal co-variation claim on the **WES-derived
iAtlas neoantigen proxies** and it came back **null** (well-powered): burden-only
co-variation, no low-rank RNA-state (`results/covariation_*`, perm p≈0.78). The
ERV proxy is degenerate (1.9% nonzero). **This is expected if the thesis is
right** — the RNA-state phenotypes are transcriptomic events a WES/annotation
pipeline cannot see. So the de-novo raw-read features the pipeline produces are
**the only substrate that can actually test the thesis.** That raises the bar on
the hand-off: the features must be **per-sample** (not group contrasts) so a
latent state can be fit and related to response per patient.

---

## The one requirement that is easy to get wrong

**PER-SAMPLE, not group-wise.** `nf-core/rnasplice`'s default output is a
responder-vs-nonresponder *contrast* (rMATS/DEXSeq group test). The modeling
layer needs a **per-sample value for every feature** so it can build S and
cross-validate response per patient. Concretely:

- Splicing → **SUPPA2 per-sample PSI** (or DEXSeq per-sample normalized exon
  usage), NOT just the group ΔPSI table.
- Intron retention → **IRFinder-S per-sample IR ratio** per intron.
- RNA editing → **per-sample Alu Editing Index (AEI)** (cohort-level site tables
  are a bonus, but AEI-per-sample is the required deliverable).
- TE/ERV → **Telescope/TEtranscripts per-sample** locus- and family-level counts.
- Fusion → **per-sample fusion count / burden** (calls per sample is fine).
- Base quant → **Salmon per-sample** gene- and transcript-level abundance (TPM).

Group-wise contrast tables are welcome *in addition* (they feed interpretation),
but the per-sample matrix is the contract.

---

## Sample key (join column)

Every matrix uses **ENA run accession** (`run_accession`, e.g. `ERR2208952`,
`SRR5088813`) as the primary sample key — this is the key in
`data/manifests/selection_manifest.csv` and `data/catalog/run_catalog.csv`,
which the modeling layer already joins to clinical/response. Do **not** invent a
new sample id. Include `cohort` as a second column.

---

## File format & location

Write to **`results/features/`** (git-ignored for large matrices; commit a small
one if it fits). One file per phenotype, **tidy-wide** (rows = samples, columns =
features), gzip-parquet preferred, CSV acceptable:

```
results/features/
  quant_gene_tpm.parquet         # rows: run_accession; cols: gene_id (+ cohort)
  quant_tx_tpm.parquet
  splicing_psi.parquet           # cols: SUPPA2 event ids
  intron_retention.parquet       # cols: intron ids; values IR ratio
  rna_editing_aei.parquet        # cols: AEI (+ optional per-region indices)
  te_locus.parquet               # cols: TE locus ids (Telescope)
  te_family.parquet              # cols: TE family (LINE/SINE/LTR/ERV…)
  fusion_burden.parquet          # cols: n_fusions (+ optional per-fusion flags)
  _feature_manifest.json         # {phenotype: {file, n_samples, n_features, tool, ref, built_at}}
```

Each matrix's **first two columns are `run_accession`, `cohort`**; all remaining
columns are numeric features. NA where a sample was not quantified. Record the
tool, reference build, and command in `_feature_manifest.json` for provenance.

---

## Pilot scope (agreed target for the deadline)

Given 64 GB RAM (serial STAR ~35 GB/sample) and ~697 GB free disk
(stream-align-delete), the **modeling pilot** needs a **small, response-balanced,
per-sample** set — not all cohorts. Priority:

1. **Gide 2019 balanced subset** (deepest reads, balanced response) — the
   modeling session's preferred pilot; ~8–16 PRE samples split R/NR.
   *If the pipeline session validates on Hugo first (smallest, 24 runs), that is
   also fine* — Hugo is response-labeled and per-sample; the modeling layer will
   consume whichever lands first.
2. Phenotype priority for the thesis (do the WES-blind ones first, since those
   are what the proxy could not test): **TE/ERV → intron retention → editing (AEI)
   → splicing PSI**, with Salmon quant as the always-produced base.

A single phenotype matrix (e.g. `te_family.parquet` + `rna_editing_aei.parquet`)
on ~12 samples is enough for the modeling session to demonstrate de-novo signal
and concordance with the proxy — it does not need all six to produce a result.

---

## Status handshake

- Pipeline session: append a line to `results/features/_feature_manifest.json`
  (or drop a `results/features/READY_<phenotype>.flag`) as each matrix lands.
- Modeling session: `src/features.py::load_available()` polls this directory and
  runs on whatever is present, degrading gracefully when a matrix is absent.

_Authored by the modeling session (`f16431e3`) 2026-07-07. Edit in place if the
pipeline session needs to renegotiate the format._
