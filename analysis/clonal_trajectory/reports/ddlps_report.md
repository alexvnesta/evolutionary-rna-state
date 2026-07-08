# Cross-cohort test and matched-bulk reality check (DDLPS liposarcoma)

## Why this cohort

cSCC Phase 0 showed clone-defined trajectory classes imprint pseudobulk expression, survive dilution, and recur across patients. Two things it could not test: whether the result is tissue-specific, and whether the imprint transfers to REAL bulk RNA (as opposed to pseudobulk of the same cells). Dedifferentiated liposarcoma (GSE221493 scRNA + GSE221492 bulk) answers both: it is a mesenchymal tumor mechanistically unrelated to squamous cSCC, and it ships matched bulk RNA-seq for the same patients. This is the first test on real bulk.

## Pipeline

102,041 cells x 28,026 genes across 31 patients, no author cell-type labels. Compartment identity was called de novo: per-patient Leiden clustering, then immune/stromal reference clusters identified by canonical markers (PTPRC, CD3D, LYZ, CD68, MS4A1, PECAM1, VWF) and used as the diploid reference for expression-based CNV clone inference. Same infercnvpy pipeline as cSCC. 29 patients had enough cells; 19 resolved >=2 malignant clones.

## Result 1 — the clone-imprint result replicates in a different tumor type

17 of 19 multi-clone tumors show subclones separated well above a within-patient permutation null in independent expression space (silhouette z = 5 to 84, median ~15). Two are genuine negatives (P18 z=0.01, P4 z=-4.35) where CNV clones carry no independent transcriptional identity, expected in a label-free de-novo pipeline. The core Phase 0 finding is not specific to cSCC.

## Result 2 — recurrent classes replicate, with an immune axis

The 55 malignant subclones (19 patients) fall into two recurrent classes: an interferon/MHC-high class (20 clones / 11 patients) and an MDM2/CDK4-amplicon + EMT-high class (35 clones / 18 patients). Both span many patients. As in cSCC, the recurrent axis couples to antigen presentation, and the coupling is sharper here and directly genetic: interferon/MHC is anti-correlated with the MDM2/CDK4 12q-amplicon program (r=-0.63, p<0.001), the defining oncogenic lesion of DDLPS. The immune-cold subclones are the ones carrying the amplified driver. In both cohorts a clonal state governs antigen presentation, but the immune-hot pole is tissue-specific (basal/stem in cSCC, adipocytic/differentiated in DDLPS) while the immune-cold pole tracks each tumor's own oncogenic program.

## Result 3 — the critical negative: naive scoring does NOT transfer to real bulk

Scoring the same programs directly in the 53 matched bulk samples breaks the malignant-compartment relationships. The scRNA clone-level anti-correlation (MDM2 vs IFN/MHC, r=-0.63) becomes a weak POSITIVE in bulk (r=+0.28); EMT and adipocytic correlations with IFN/MHC also flip sign. The cause is compartment mixing (a Simpson's paradox): bulk IFN/MHC is almost entirely explained by an immune-infiltration proxy (r=+0.94, ~88% of variance, p=3.5e-25), so the bulk immune signal reads out infiltrating lymphocytes and myeloid cells, not malignant antigen-presentation state. Adjusting IFN/MHC for infiltration removes the confounded positive (r drops to +0.16, n.s.) but does not recover the -0.63 malignant-intrinsic value: simple linear deconvolution is insufficient.

Patient-level transfer of malignant-intrinsic programs from scRNA clones to matched bulk is also weak: MDM2-amplicon r=+0.08, adipocytic r=+0.33, EMT r=+0.26, none significant at n=19. Two contributors are plausible and not separable here: scRNA and bulk are likely from different tumor pieces (spatial sampling mismatch), and collapsing each patient's clones to a mean discards the compositional signal a real classifier would exploit.

## Result 4 — supervised bulk recovery is untestable at this n

The project's actual proposed method is supervised: train a classifier on bulk with the scRNA-derived clone class as the label, letting the model find compartment-robust features rather than relying on hand-picked programs. Tested directly on the 19 matched bulk samples (dominant clone class per patient, 13:6 balance, top-3000-variance features, leave-one-out CV, L2-regularized logistic regression): balanced accuracy 0.55, AUC 0.55, permutation p=0.274. Not significant.

This is not evidence the method fails, it is evidence the test is underpowered. n=19 with 3000 features is a severe small-n/large-p regime, and it sits below the usable floor from the earlier power analysis (moderate effect needs N~40 for balanced accuracy >0.62, N~67 for >0.82). On top of that, the bulk and scRNA are likely different tumor pieces, adding label noise. The supervised approach cannot be adjudicated at this sample size; it needs the larger matched cohorts or same-tissue sampling.

## What this changes

The pseudobulk and dilution successes in cSCC established the signal is recoverable in principle when the readout is the same cells. This matched-bulk test tempers that: with a naive pathway-score readout on real, independently-sampled bulk, the malignant-compartment class signal is masked by immune admixture and sampling. The implication for the project design is concrete and useful:

1. The bulk classifier cannot use raw pathway scores of immune-related programs; those are dominated by infiltration. It must either (a) train on malignant-compartment-specific features recovered by expression deconvolution (CIBERSORTx / BayesPrism), or (b) learn the class directly as a supervised label from bulk with the scRNA-derived label as ground truth, letting the model find compartment-robust features rather than hand-chosen programs.
2. Genetic-amplicon signals (MDM2/CDK4 in DDLPS), being compartment-specific and dosage-driven, are the most promising bulk-transferable features and should be prioritized over immune programs for the malignant-state readout.
3. Matched same-tissue scRNA+bulk (not different pieces) is needed to cleanly separate sampling mismatch from method failure; the cSCC cohort's lack of bulk and DDLPS's likely piece-mismatch both limit this.

## Bottom line

The biology replicates across two unrelated tumor types, which strengthens the thesis that a recurrent clonal state governs antigenicity. The engineering is harder than the pseudobulk proxy implied: real-bulk transfer needs deconvolution or supervised label-recovery, not naive program scoring. This is the reality check that should precede any scale-up, and it reshapes Phase 2 from pathway-scoring to deconvolution-based or supervised classification.
