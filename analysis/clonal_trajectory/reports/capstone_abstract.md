# Clonal-evolutionary immune state imprints the tumor transcriptome, but reading it from bulk is compartment-limited

## Abstract

**Premise.** The project thesis holds that a tumor's clonal-evolutionary trajectory sets an
underlying RNA state that governs antigenicity and immune-checkpoint-blockade (ICB)
response. We tested a specific, tractable version: that clone-defined transcriptional states
imprint the bulk-readable transcriptome strongly enough to be recovered by a supervised
classifier and, ultimately, to predict ICB response. The design is supervised label
recovery (treat the single-cell-derived trajectory class as a training label), not de-novo
phylogeny reconstruction from bulk, which is underdetermined.

**Clone-imprint replicates across four cancers.** Expression-based CNV clone inference
(infercnvpy) on single-cell cohorts from four unrelated tumor lineages — cutaneous squamous
carcinoma (cSCC), dedifferentiated liposarcoma (DDLPS), acute myeloid leukemia (AML), and
lung adenocarcinoma (LUAD) — recovered malignant subclones whose partition, validated in an
independent highly-variable-gene expression space against a within-patient permutation null,
scored well above null in the large majority of multi-clone tumors (8/8 cSCC, 17/19 DDLPS,
7/7 AML, 15/16 LUAD; median silhouette z 10-24). Within-patient subclone separation is a
substantial fraction of between-patient separation. The imprint is real, genome-wide,
non-circular, and recurrent.

**The immune-class axis is tumor-type-dependent.** In cSCC, DDLPS, and AML the clones
organize into two discrete recurrent classes with an immune-hot pole (about 36% of clones)
coupled to antigen presentation. In LUAD the interferon/MHC axis is continuous rather than a
discrete hot/cold split. This is a boundary condition on the model, not a failure of it.

**Supervised recovery is powered, with a structural ceiling.** A harmonized immune-state
label, pooled across the three discrete-axis cohorts (N=99 clone-level pseudobulk), is
recovered at grouped-by-patient cross-validated balanced accuracy 0.91 (AUC 0.94,
permutation p=0.005), directly confirming that an earlier n=19 failure was sample size, not
concept. Adding LUAD (N=145) keeps the result significant but lowers it to 0.79, because the
graded axis forced into a median-split label adds noise. The binding constraint past the
power floor is how discretely a tissue's immune axis is structured, not sample count.
Leave-one-tissue-out transfer stays near chance except for AML, so the method is
within-tissue, not zero-shot.

**Reading the malignant-intrinsic state from real bulk is the blocker.** On real DDLPS
bulk, naive program scoring sign-flips (Simpson's paradox: bulk interferon/MHC is about 88%
explained by immune-cell infiltration, r=+0.94), and NNLS composition deconvolution does not
recover the within-malignant axis. Restricting to malignant-compartment-specific features
lifts real-bulk balanced accuracy from 0.55 to 0.64 (not significant at n=19) and improves
low-purity separability in a controlled dilution test, confirming the mechanism while
remaining under-powered. At the ICB endpoint (Riaz melanoma bulk, 49 pre-treatment samples,
10 responders), neither the malignant-intrinsic interferon/MHC signature (AUC 0.56) nor an
immune-infiltration proxy (AUC 0.62) predicts response, and the two correlate at r=0.84: a
coherent negative that localizes the blocker to the bulk-transfer step, not the biology.

**Conclusion.** The clonal-evolutionary immune state imprints the tumor transcriptome
robustly and generalizably, and is recoverable by a properly-powered supervised classifier
within a tumor type. Translating that to bulk-only ICB prediction requires compartment
resolution that bulk alone does not provide. The two experiments that would close the gap —
a matched single-cell-plus-bulk-plus-ICB cohort at adequate N, or a compartment-specific
feature set certified at N about 40 in one tissue — are gated on data availability, not
method.

## Figure

`figures/capstone_synthesis.png` : (A) replication across four cancers, (B) the powered
classifier and its LUAD-induced ceiling, (C) leave-one-tissue-out transfer, (D) the
bulk-transfer barrier through to the ICB endpoint.
