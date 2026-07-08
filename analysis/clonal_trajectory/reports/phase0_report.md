# Phase 0 results: does a clonal-trajectory class imprint bulk expression? (cSCC, executed)

## The gating question

The feasibility of the whole program hinged on one unmeasured number: how strongly a clone-defined trajectory class imprints whole-tumor (pseudobulk) expression. If clones are transcriptionally near-identical once mixed, a bulk classifier can never recover the class and the project stops. If the imprint is strong and recurs across patients, the supervised label-recovery design is viable. This phase measured it on real data.

## What was actually run

Cohort: cutaneous squamous cell carcinoma, GSE144240 (scRNA) + GSE144236 (processed counts) + GSE144237 (WES mutation table). 48,164 cells x 32,738 genes across 10 patients, with author-provided cell-type annotations (malignant keratinocyte vs immune/stromal). No BAMs are deposited, so SNV-anchored clone assignment (cardelino) is not possible for this cohort; clones were called from expression-based copy-number inference (infercnvpy) using immune and stromal cells as the diploid reference. Pipeline: normalize -> infer per-cell CNV -> Leiden-cluster cells in CNV space -> label malignant clones (epithelial-dominated, CNV signal > 1.5x reference) -> pseudobulk each clone.

## Result 1 — clones and subclonal structure are real

8 of 10 patients showed clear whole-tumor malignant CNV signal (epithelial/reference CNV ratio 1.58-2.78x). Two tumors sat below the 1.5x whole-tumor threshold: P3 (ratio 1.29, and no malignant clone called) contributed just 223 epithelial cells, and P1 (ratio 1.20) had a diffuse epithelial CNV profile in which only a single high-CNV subcluster passed the per-clone threshold. Both P1 and P3 therefore contributed 0 or 1 clone and did not enter the multi-clone subclonal analyses (Results 2-5), which rest on the 7-8 tumors with >=2 malignant clones. Those informative tumors resolved into 2-4 malignant CNV subclones each (25 subclones total across 9 patients, 25,142 malignant cells labeled). The WES table independently confirms canonical cSCC driver architecture in this cohort (TP53 12, NOTCH1 13, NOTCH2 15, CDKN2A 6, plus KMT2C/D, FAT1, HRAS, ARID2), so the malignant biology is what it claims to be.

## Result 2 — strong pseudobulk imprint (effect size)

In the top-variance pseudobulk gene space, subclones of the same tumor are separated by a median 1.10 SD per gene, nearly as much as clones from different tumors (1.29 SD). The top signal genes shift by ~3.5 SD between subclones of a patient. On the power curve calibrated earlier (weak=0.6, moderate=1.0, strong=1.5 SD on the signal genes), this sits at or above the strong end. In the earlier simulation, a strong imprint was recoverable at 3-class balanced accuracy ~0.73 even at N=25 and ~0.97 by N=67.

## Result 3 — not a CNV-clustering artifact

Because clones were defined from expression-derived CNV, the separation could be circular. Two guards argue it is not.

Guard A (genomic diffuseness): between-subclone expression variance is only ~10% concentrated on the single most-variable chromosome (median 0.104 across patients), against a 4% uniform baseline and a >30% threshold for CNV-arm-dominated (circular) signal. The imprint is genome-wide, not confined to the aberrant arms that seeded the clustering.

Guard B (independent readout + permutation null): re-clustering malignant cells in highly-variable-gene expression space (PCA), independent of the CNV features, the clone partition scored a silhouette 6.4 to 63.5 SD above a within-patient random-label null (median z ~24) in every one of the 7 multi-clone patients. Subclones carry a reproducible transcriptional identity that survives on features they were not defined by.

## Result 4 — two recurrent classes, immune-coupled

Clustering the 25 subclones on hallmark transcriptional programs yields a stable two-class axis that recurs across patients: a differentiated class (16 clones / 9 patients, high epithelial-differentiation program) and a basal/stem-like class (9 clones / 6 patients, high basal/stem program). Both span multiple patients, so they are transferable classes rather than per-patient idiosyncrasies. The axis couples to antigen presentation: the basal/stem-like class carries significantly higher interferon/MHC program activity (r=+0.49, p=0.012), and differentiation anti-correlates with it (r=-0.41, p=0.044). That coupling is the mechanistic hook the project thesis predicts, a clonal-state axis wired to the antigen-presentation machinery that governs checkpoint-blockade response.

## Result 5 — the imprint survives realistic tumor-purity dilution

The pseudobulk in Results 2-4 is malignant-cell-only, an optimistic proxy for real bulk. To test the true endpoint, each patient's malignant clones were remixed with that patient's own non-malignant (immune/stromal) cells at controlled purity, and class separability was re-measured. The signal degrades gracefully rather than collapsing: mean clone-class silhouette runs 0.78 (pure) -> 0.61 (60% purity) -> 0.42 (30%) -> 0.19 (10%). At 60% purity, typical for solid-tumor bulk, ~78% of the pure-compartment separability is retained; the aggregate signal stays clearly above the ~0 noise floor even at 10% purity (mean 0.19), with only the single weakest tumor near indistinguishability. So the clone-class imprint is recoverable from a realistically-admixed bulk sample, not just from sorted malignant cells. This is the single most important robustness check for the bulk-classifier endpoint, and it passes.

## Scope and limitations

Establishes: clone-defined trajectory classes imprint whole-tumor expression strongly and genome-wide, the imprint is not a CNV-clustering artifact, the classes recur across patients, the recurrent axis is coupled to an ICB-relevant immune program, and the imprint survives physiological stromal/immune dilution. On the calibrated power curve this places the project in the favorable regime where the bulk classifier is trainable at the N the public matched cohorts provide. The Phase 0 gate is passed.

Does not establish: (1) any link to actual ICB outcome (the immune correlate is a program score, not measured response). (2) Transfer across cancer types: this is cSCC only, and the classes may be tissue-specific. (3) SNV-level clone identity: CNV clones are a coarser partition than a full mutational phylogeny. (The stroma-dilution concern from an earlier draft is now addressed by Result 5.)

## Next steps

1. DONE — Dilution test (Result 5): the imprint survives physiological purity. No longer a gating unknown.
2. Cross-cohort recurrence: repeat the pipeline on a second solid-tumor matched cohort (DDLPS GSE221492/3) and test whether the differentiation/stem axis or an analogous recurrent axis appears, probing the tissue-specificity limitation.
3. Bulk classifier + ICB test: with the class labels now in hand, train the pseudobulk-to-class classifier and apply it to the bulk ICB cohorts (Riaz, Hugo, Auslander) to test association with response, the Phase 2-3 core.
