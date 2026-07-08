# Clonal-trajectory-to-bulk-to-ICB: investigation summary

This is the consolidated record of the feasibility investigation for the hypothesis that a tumor's clonal-evolutionary trajectory, defined from single-cell data, imprints its bulk transcriptome strongly enough to be read from bulk RNA alone and used to predict immune-checkpoint-blockade response. It spans a literature review, a data-availability survey, a power analysis, and two experimental cohorts (cSCC and DDLPS liposarcoma).

## The idea, and where it stands

Design: define discrete clonal-trajectory classes from matched scRNA + genomic data, transfer the class label into bulk RNA as a supervised classification target, then test whether the inferred class predicts ICB response. The load-bearing and novel step is transferring a categorical trajectory label into bulk (not de-novo phylogeny reconstruction, which is underdetermined from a bulk clone mixture). Framed as supervised label-recovery, it is tractable in principle.

Verdict after this investigation: the underlying biology is real and replicates across two unrelated tumor types. The bulk-transfer engineering is materially harder than an optimistic pseudobulk proxy suggested, and the decisive test is currently gated by sample size, not by any demonstrated failure of the concept.

## What the literature established

No published method does exactly this. The closest templates: Zhang 2022 (Stem.Sig) derived a single-cell signature and validated it across bulk ICB cohorts, the direct precedent for scRNA-to-bulk-to-ICB signature transfer, but it transfers a continuous pathway score, not a categorical trajectory class. Liao 2020 and REVOLVER (Caravagna 2018) define recurrent evolutionary subtypes from multi-region TRACERx data. McGranahan 2016, Dijkstra 2025, and Alban 2024 establish the clonal-architecture-to-ICB link. Three of the four required links are precedented; the categorical-label transfer into bulk is the genuinely new piece. Full review in clonal_evolution_bulkRNA_ICB_review.md.

## What the data survey established

No single deposit carries all four layers (scRNA + genomic + bulk + ICB outcome). The workable route is federated: training cohorts with matched scRNA + genomic to define and label classes, then separate bulk-plus-ICB cohorts to test. Matched-cohort patient counts are the binding constraint: cSCC 9, DDLPS 16, with AML scDNA (~38) and small others. The calibrated power analysis (on real Hugo melanoma bulk covariance) set the usable floor: a moderate within-sample effect needs N~40-67 matched patients for balanced accuracy 0.62-0.82; a strong effect is usable from N~25. Manifest in clonal_evo_data_manifest.csv, power grid in clonal_evo_power_grid.csv.

## What the experiments established

### cSCC (GSE144240), Phase 0 gate — PASSED

Expression-based CNV clone inference on 48,164 cells resolved 2-4 malignant subclones per tumor (25 subclones across 9 patients; 8 of 10 tumors cleared the whole-tumor CNV threshold, with P1 and P3 contributing 0-1 clones and excluded from multi-clone analyses). Findings:
- The clone partition imprints pseudobulk strongly (within-patient subclone separation 1.10 SD/gene, ~3.5 SD on top signal genes, at or above the strong end of the power axis).
- The signal is not a CNV-clustering artifact: only ~10% of between-clone variance concentrates on one chromosome, and in independent expression space the clones separate 6-63 SD above a permutation null (median z 24) in every multi-clone tumor.
- Subclones fall into two recurrent, transferable classes (differentiated vs basal/stem-like) that couple to antigen presentation (basal/stem vs IFN/MHC r=+0.49, p=0.012).
- The imprint survives realistic tumor-purity dilution: 78% of separability retained at 60% purity, positive down to ~30%.

### DDLPS liposarcoma (GSE221493 scRNA + GSE221492 bulk) — biology replicates, bulk-transfer negative

De-novo clone inference on 102,041 cells (31 patients, no author labels) resolved >=2 clones in 19 patients. Findings:
- The clone-imprint result replicates: 17 of 19 tumors show clones separated above null (median z ~15) in a mechanistically unrelated tumor type.
- Recurrent classes replicate with a sharper, genetic immune axis: IFN/MHC anti-correlates with the MDM2/CDK4 12q-amplicon (the defining DDLPS driver) at r=-0.63, p<0.001. Immune-cold clones carry the amplified oncogene, a direct genotype-to-antigenicity coupling that matches the project thesis.
- The critical negative: on the 53 matched REAL bulk samples, naive program scoring breaks. The malignant anti-correlation flips to a weak positive (Simpson's paradox) because bulk IFN/MHC is almost entirely explained by immune infiltration (r=+0.94, ~88% of variance), not malignant antigen presentation. Patient-level scRNA-to-bulk program transfer is weak (r=0.08-0.33, n.s.).
- The supervised recovery (the actual proposed method) is untestable at n=19: balanced accuracy 0.55, p=0.27, but this n is below the power floor and the bulk is a different tumor piece than the scRNA.

### Constructive follow-up — compartment-specific features rescue the signal

The DDLPS negative pointed to a remedy: use malignant-compartment-specific features instead of whole-transcriptome pathway scores. Tested directly in the controlled cSCC dilution setup (where mixing is known, removing the piece-mismatch confound), restricting to the 224 genes specifically enriched in the malignant compartment makes the class signal markedly more dilution-robust. At low purity, where admixture masks the signal, compartment-specific features recover meaningfully more separability than all-gene features (silhouette 0.30 vs 0.18 at 10% purity, a 64% relative gain; 0.51 vs 0.43 at 30%, a 20% gain), crossing over to an advantage below ~65% purity, exactly the realistic-bulk regime. This is direct, controlled evidence that the deconvolution / compartment-restriction strategy is the right fix for the admixture confound. Detail in phase0_compartment.png.

## The honest synthesis

Two things are now established with real data. First, the core biological premise holds and generalizes: clone-defined states imprint the transcriptome genome-wide, form recurrent cross-patient classes, and those classes couple to antigen presentation, in two unrelated cancers. Second, reading that malignant-compartment state from real bulk is confounded by immune admixture and tissue sampling, and the naive pathway-score readout fails. The supervised, deconvolution-aware version of the method is the right next design, but it cannot be judged at the sample sizes reached here.

The investigation therefore converts an open feasibility question into a concrete, de-risked plan: the concept is sound and the biology is banked, the method must be deconvolution-based or supervised rather than pathway-score-based, and the binding requirement is N>=40 matched patients (ideally same-tissue scRNA+bulk) before the bulk classifier can be adjudicated.

## Recommended next steps

1. Deconvolution-aware features: rerun the bulk readout through CIBERSORTx or BayesPrism to isolate the malignant-compartment expression, then re-test class transfer. Directly addresses the immune-admixture confound, and the controlled cSCC test already shows compartment-restricted features recover the low-purity signal, so this is a validated direction, not a gamble.
2. Assemble N>=40 matched scRNA+bulk (pool cSCC + DDLPS + AML + additional GEO cohorts) to bring the supervised test above the power floor.
3. Prioritize genetic-amplicon / copy-number-driven features (compartment-specific, dosage-robust) over immune programs for the malignant-state readout.
4. Only after a supervised classifier clears cross-validation on matched data, apply it to the bulk ICB cohorts (Riaz, Hugo, Auslander) for the outcome test.

## Artifact index

Reviews and reports: clonal_evolution_bulkRNA_ICB_review.md, phase0_report.md, ddlps_report.md, clonal_evo_feasibility_protocol.md.
Figures: phase0_summary.png, phase0_validation.png, phase0_classes.png, phase0_dilution.png, phase0_compartment.png, ddlps_replication.png, ddlps_bulk_transfer.png, clonal_evo_power.png, clonal_evo_data_flow.png.
Data: clonal_evo_data_manifest.csv, clonal_evo_power_grid.csv, clonal_evo_matched_counts.csv, cscc_infercnv_summary.csv, cscc_tumor_clone_assignments.csv, phase0_validation.csv, clone_programs_classified.csv, phase0_dilution.csv, ddlps_infercnv_summary.csv, ddlps_validation.csv, ddlps_clone_programs_classified.csv, ddlps_bulk_scores.csv, ddlps_loo_preds.csv.
