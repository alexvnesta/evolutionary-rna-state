# Clonal-trajectory analysis

Testing whether a tumor's clonal-evolutionary trajectory, defined from single-cell
data, imprints its bulk transcriptome strongly enough to be read from bulk RNA alone
and used to organize immune-checkpoint-blockade response. This is the single-cell arm
of the project's latent evolutionary RNA-state thesis: the clone-defined state is one
concrete, measurable instance of the latent state S.

## Design

Supervised label-recovery, not de-novo phylogeny. Define discrete clonal-trajectory
classes from matched scRNA (+ genomic where available), transfer the class label into
bulk RNA as a classification target, then test whether the inferred class tracks the
immune / ICB axis. De-novo phylogeny reconstruction from a bulk clone mixture is
underdetermined and is not attempted.

## What is here

- `reports/` : the written record.
  - `clonal_evolution_bulkRNA_ICB_review.md` : literature review and precedent.
  - `clonal_evo_feasibility_protocol.md` : the 4-phase protocol and power analysis.
  - `phase0_report.md` : cSCC Phase 0 (clone inference, effect size, validation, classes, dilution).
  - `ddlps_report.md` : DDLPS cross-cohort replication and the matched-bulk reality check.
  - `investigation_summary.md` : integrated top-level synthesis.
- `figures/` : publication-style figures for each result.
- `tables/` : clone assignments, CNV summaries, validation z-scores, program classes, power grid, bulk scores.
- `code/` : the analysis scripts (inferCNV clone inference, effect size, validation, programs, dilution, compartment test; cSCC and DDLPS).

## Results in one paragraph

The clone partition imprints pseudobulk expression genome-wide (not a CNV-clustering
artifact: independent-space silhouette 6-63 SD above a permutation null), forms
recurrent cross-patient classes coupled to antigen presentation, and survives
tumor-purity dilution. This replicates in two unrelated tumor types (cSCC squamous,
DDLPS mesenchymal). The load-bearing caveat: on real matched bulk (DDLPS), naive
pathway-score transfer breaks because bulk immune signal is dominated by infiltration
(Simpson's paradox); restricting to malignant-compartment-specific features recovers
the low-purity signal in controlled mixtures. A properly powered supervised bulk test
needs N>=40 matched patients.

## Data

Public GEO cohorts, code-and-metadata-only in this repo (large scRNA matrices and
derived AnnData are gitignored and re-derivable from the cited accessions):
cSCC GSE144240/GSE144236/GSE144237; DDLPS GSE221492 (bulk) + GSE221493 (scRNA);
AML GSE116256 (van Galen scRNA).

## Provenance

Produced in Claude Science. Every figure and table is a saved artifact; the reports
were style-checked and audited for numeric fidelity.
