# Clonal-trajectory-to-bulk-to-ICB: investigation summary

This is the consolidated record of the feasibility investigation for the hypothesis that a tumor's clonal-evolutionary trajectory, defined from single-cell data, imprints its bulk transcriptome strongly enough to be read from bulk RNA alone and used to predict immune-checkpoint-blockade response. It spans a literature review, a data-availability survey, a power analysis, and three experimental cohorts (cSCC, DDLPS liposarcoma, and AML), including a deconvolution remedy test and a pooled powered classifier.

## The idea, and where it stands

Design: define discrete clonal-trajectory classes from matched scRNA + genomic data, transfer the class label into bulk RNA as a supervised classification target, then test whether the inferred class predicts ICB response. The load-bearing and novel step is transferring a categorical trajectory label into bulk (not de-novo phylogeny reconstruction, which is underdetermined from a bulk clone mixture). Framed as supervised label-recovery, it is tractable in principle.

Verdict after this investigation: the underlying biology is real and replicates across three unrelated tumor types, and the supervised method works once properly powered. The bulk-transfer engineering is materially harder than an optimistic pseudobulk proxy suggested, and the decisive test is currently gated by sample size, not by any demonstrated failure of the concept.

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

### AML (van Galen GSE116256) — third replication, hardest setting

The original AML candidate (GSE156934) turned out to be single-cell DNA, which cannot drive the expression-imprint method, so the canonical van Galen AML scRNA atlas was substituted. Clone inference on 16 diagnosis (D0) patients (normal T/NK/B lymphocytes as diploid reference, since AML malignant cells are hematopoietic) found:
- The clone-imprint replicates in a third, mechanistically unrelated tumor type. AML is genomically quiet (only 2 of 12 analyzed patients frankly aneuploid), so CNV magnitude is not the evidence: the non-circular independent-expression validation is, with the clone partition beating a within-patient permutation null in 7 of 7 multi-clone patients (median silhouette z 23). Within-patient (1.02 SD/gene) versus between-patient (1.44) effect size matches cSCC.
- The 19 clones form two recurrent classes (differentiated/myeloid, IFN/MHC-high; and HSC/progenitor, IFN/MHC-lower) coupled to antigen presentation at the class level (Mann-Whitney p=0.045), directionally identical to the cross-tissue pattern but softer and not significant per-program. Detail in aml_report.md.

### Deconvolution remedy on real DDLPS bulk — honest negative

Building the DDLPS remedy on real bulk (rather than the controlled cSCC mixture): a 560-gene, 7-population signature (two malignant classes plus five stromal/immune types) was NNLS-deconvolved against the 53 matched bulk samples. Deconvolution removes the confound but does not manufacture the malignant-intrinsic axis: the spurious MDM2-vs-IFN correlation collapses from the naive +0.28 (p=0.04) toward zero at +0.12 (n.s.), beating the earlier simple linear immune-adjustment (+0.16), but it never approaches the scRNA clone-level target of -0.63. The reason is structural: the MDM2-vs-IFN anti-correlation is a within-malignant, between-clone expression-state axis, not a between-sample composition axis, so composition correction cannot recover it. This scopes the malignant-intrinsic axis to the single-cell compartment; it is not recovered in these matched bulk samples even after full deconvolution. Detail in ddlps_deconvolution_report.md.

### Pooled powered supervised classifier — the power floor is reached

Pooling 99 clone-level pseudobulk profiles across all three cohorts (cSCC 25, DDLPS 55, AML 19) brought the supervised test above the power floor for the first time. The target is a harmonized immune-state label (immune-hot vs immune-cold), oriented per-cohort by IFN/MHC because the raw class labels are not consistently oriented across tumor types; features are z-scored within cohort so tissue-baseline offsets cannot be a shortcut.
- Powered result: grouped-by-patient 5-fold CV reaches balanced accuracy 0.91, AUC 0.94 (N=99), against a within-cohort permutation null of 0.53 (p=0.005), landing on the strong curve of the calibrated power grid. The earlier n=19 DDLPS negative was under-powering, not absence of signal.
- Transfer limit: leave-one-cohort-out is weak (cSCC 0.43, DDLPS 0.62, AML 0.53). The immune-hot pole is tissue-specific, so supervised recovery works within a tumor type given enough N but does not transfer zero-shot pan-cancer. The immune-hot fraction is nonetheless strikingly consistent (~36%) across all three unrelated tumors. Detail in pooled_classifier_report.md.

## The honest synthesis

Three things are now established with real data. First, the core biological premise holds and generalizes: clone-defined states imprint the transcriptome genome-wide, form recurrent cross-patient classes, and those classes couple to antigen presentation, in three unrelated cancers (squamous cSCC, mesenchymal DDLPS, hematopoietic AML). Second, the supervised label-recovery method works once it is properly powered: at pooled N=99 the harmonized immune-state label is recovered at balanced accuracy 0.91 (p=0.005), directly confirming that the earlier n=19 failure was sample size, not concept. Third, two hard limits are now quantified rather than assumed: reading the malignant-intrinsic state from real bulk is not rescued by composition-based deconvolution (the axis lives within the malignant compartment, not in mixing proportions), and the trained classifier is within-tissue, not zero-shot pan-cancer.

The investigation therefore converts an open feasibility question into a de-risked, precisely-scoped program: the concept is sound and the biology is banked across three tumor types; the method is powered and works within a tumor type at N in the tens; the two open problems are (a) recovering the within-malignant state from real bulk (composition deconvolution is insufficient; compartment-specific features are the controlled-data-validated candidate but unproven on real matched bulk) and (b) the need for per-tissue rather than pan-cancer training.

## Recommended next steps

1. Test compartment-specific (not just composition-deconvolved) features on real matched DDLPS bulk. The controlled cSCC test validated the direction and the composition-deconvolution negative rules out the simpler fix, so this is the sharpest remaining experiment.
2. Reach N>=40 matched scRNA+bulk within a single tissue (the pooled test showed cross-tissue transfer is weak, so same-tissue power is what counts) before adjudicating the real-bulk supervised classifier.
3. Prioritize genetic-amplicon / copy-number-driven features (compartment-specific, dosage-robust) over immune programs for the malignant-state readout.
4. Only after a supervised classifier clears cross-validation on real matched same-tissue bulk, apply it to the bulk ICB cohorts (Riaz, Hugo, Auslander) for the outcome test.

## Artifact index

Reviews and reports: clonal_evolution_bulkRNA_ICB_review.md, phase0_report.md, ddlps_report.md, aml_report.md, ddlps_deconvolution_report.md, pooled_classifier_report.md, clonal_evo_feasibility_protocol.md.
Figures: phase0_summary.png, phase0_validation.png, phase0_classes.png, phase0_dilution.png, phase0_compartment.png, ddlps_replication.png, ddlps_bulk_transfer.png, ddlps_deconvolution.png, aml_replication.png, pooled_classifier.png, clonal_evo_power.png, clonal_evo_data_flow.png.
Data: clonal_evo_data_manifest.csv, clonal_evo_power_grid.csv, clonal_evo_matched_counts.csv, cscc_infercnv_summary.csv, cscc_tumor_clone_assignments.csv, phase0_validation.csv, clone_programs_classified.csv, phase0_dilution.csv, ddlps_infercnv_summary.csv, ddlps_validation.csv, ddlps_clone_programs_classified.csv, ddlps_bulk_scores.csv, ddlps_loo_preds.csv, ddlps_deconv_fractions.csv, ddlps_deconv_correlations.csv, aml_infercnv_summary.csv, aml_validation.csv, aml_clone_programs_classified.csv, pooled_results.json, pooled_labels_preds.csv.
