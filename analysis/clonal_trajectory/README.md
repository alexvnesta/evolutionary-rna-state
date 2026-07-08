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

## Status: complete (data-gated)

Every analysis runnable on public data is done, written up, and committed. The one
remaining scientific step is external, a matched scRNA+bulk+ICB cohort at N>=40 in one
tissue, which is a data-acquisition problem, not a method gap. Start with
`reports/capstone_abstract.md` and `figures/capstone_synthesis.png` for the whole story.

## What is here

- `reports/` : the written record (15 files). Entry points first:
  - `capstone_abstract.md` : one-page synthesis of the complete investigation.
  - `investigation_summary.md` : integrated top-level narrative (current: v10).
  - `methods.md` : full publication-grade methods, grounded in the committed code.
  - `clonal_evolution_bulkRNA_ICB_review.md` : literature review and precedent (Zhang 2022 Stem.Sig template).
  - `clonal_evo_feasibility_protocol.md` : the 4-phase protocol and power analysis.
  - `phase0_report.md` : cSCC Phase 0 (clone inference, effect size, validation, classes, dilution).
  - `ddlps_report.md`, `ddlps_deconvolution_report.md`, `ddlps_compartment_bulk_report.md` : DDLPS replication, deconvolution (negative), compartment features on real bulk.
  - `aml_report.md` : AML (van Galen) replication in a genomically quiet tumor.
  - `luad_report.md` : LUAD 4th-tissue replication and the continuous-immune-axis boundary condition.
  - `pooled_classifier_report.md`, `pooled4_classifier_report.md` : the powered supervised classifier (3- and 4-cohort).
  - `riaz_icb_report.md` : the ICB outcome endpoint (coherent negative).
  - `cohort_survey.md` : GEO survey establishing the matched-cohort gap.
- `figures/` : 18 publication-style figures, one per result plus `capstone_synthesis.png`.
- `tables/` : clone assignments, CNV summaries, validation z-scores, program classes, per-cohort clone pseudobulk (committed so the classifier reproduces without re-running inferCNV), power grid, pooled results.
- `code/` : analysis scripts by cohort (`aml/`, `luad/`, `icb/`) and phase; see `code/README.md`.

## Results in one paragraph

The clone partition imprints pseudobulk expression genome-wide (not a CNV-clustering
artifact: independent-space silhouette well above a permutation null), forms recurrent
cross-patient classes coupled to antigen presentation, and survives tumor-purity dilution.
This replicates in four unrelated tumor types (cSCC squamous, DDLPS mesenchymal, AML
hematopoietic, LUAD glandular), though LUAD's immune axis is continuous rather than the
discrete hot/cold split of the other three. A supervised classifier recovers the harmonized
immune-state label at grouped cross-validated balanced accuracy 0.91 (N=99, p=0.005) across
the three discrete-axis cohorts; adding LUAD (N=145) keeps it significant at 0.79, showing
the ceiling is set by how discretely a tissue's immune axis is structured, not by N. The
load-bearing caveat: on real bulk, naive pathway-score transfer breaks because bulk immune
signal is dominated by infiltration (Simpson's paradox), NNLS deconvolution does not recover
the within-malignant axis, and restricting to malignant-compartment-specific features
recovers the signal only partially (not significant at n=19). At the ICB endpoint (Riaz
melanoma bulk, 49 pre-treatment, 10 responders), neither the malignant-intrinsic IFN/MHC
signature (AUC 0.56) nor an infiltration proxy (AUC 0.62) predicts response and the two
correlate at r=0.84: a coherent negative that localizes the blocker to the bulk-transfer
step, not the biology.

## Data

Public GEO cohorts, code-and-metadata-only in this repo (large scRNA matrices and
derived AnnData are gitignored and re-derivable from the cited accessions):
cSCC GSE144240/GSE144236/GSE144237; DDLPS GSE221492 (bulk) + GSE221493 (scRNA);
AML GSE116256 (van Galen scRNA); LUAD GSE131907 (Kim); ICB bulk Riaz GSE91061,
Hugo GSE78220, Auslander GSE115821. Per-cohort clone pseudobulk `.npy` files under
`tables/` are committed (about 12 MB total) so the pooled classifier reproduces without
re-running inferCNV.

## Regeneration

The `evolutionary-rna-state-refresh` skill has a clonal-trajectory track documenting how to
fold in a new matched cohort when one appears. This analysis has fixed public inputs, so it
is deliberately excluded from the daily input-diff and regenerated on-demand only.

## Provenance

Produced in Claude Science. Every figure and table is a saved artifact; the reports
were style-checked and audited for numeric fidelity.
