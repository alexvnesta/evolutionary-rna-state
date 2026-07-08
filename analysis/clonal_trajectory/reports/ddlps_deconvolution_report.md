# Deconvolution remedy for the DDLPS matched-bulk sign-flip

## Question

In prior work on this project, DDLPS liposarcoma single-cell RNA-seq (GSE221493)
clone inference found that malignant subclones split into two classes: an
**interferon/MHC-high** class (class 1) and an **MDM2/CDK4-amplicon + EMT-high**
class (class 2). At the clone level, the MDM2-amplicon program is strongly
**anti-correlated** with the IFN/MHC program (**r = −0.63**, n = 55 clones).

On the 53 matched **real bulk** samples (GSE221492), naive whole-bulk program
scoring **flips the sign to r = +0.28** (p = 0.041). The diagnosis was a
compartment-mixing Simpson's paradox: bulk IFN/MHC is dominated by immune
infiltration rather than malignant-intrinsic biology, so the two programs
co-vary positively across samples of differing purity even though they
anti-vary *within* the malignant compartment.

**Remedy hypothesis tested here:** deconvolve the bulk into its constituent
cell populations, recover the malignant-compartment expression, and test whether
the malignant-intrinsic MDM2-vs-IFN negative sign returns.

## Methods

1. **Per-cell labels.** Re-ran the deterministic (`np.random.seed(0)`) DDLPS
   clone-inference pipeline: per patient (≥1500 cells, 29 patients), Leiden
   clustering, marker-based identification of the diploid immune/stromal
   reference (PTPRC, CD3D, LYZ, CD68, MS4A1, PECAM1, VWF), inferCNV
   (`infercnvpy`, window 100), Leiden on the CNV-PCA embedding to call malignant
   clones. This labeled **100,760 cells** (35,298 malignant). Malignant clones
   were joined to their `class_k2` label from the prior classified clone table
   (all 55 prior clones matched); **30,586 malignant cells** carried a class
   (class 1 IFN/MHC-high: 9,297; class 2 MDM2-amplicon: 21,289). Non-malignant
   reference cells were assigned a broad type (T cell, Myeloid, B cell,
   Endothelial, Fibroblast) by marker score.

2. **Signature matrix S.** Mean CPM profiles for 7 populations
   (malignant_class1, malignant_class2, T, Myeloid, B, Endo, Fibroblast) over
   the 19,221 genes shared between scRNA and bulk. Marker selection kept the top
   ~80 population-specific genes each (log-fold-change vs. the best competing
   population, expressed >5 CPM), yielding **560 signature genes**. Condition
   number **50.9** (well-conditioned).

3. **NNLS deconvolution.** For each of the 53 CPM-normalized bulk samples,
   estimated non-negative population coefficients (`scipy.optimize.nnls`)
   against S, normalized to fractions (CIBERSORT-style support).

4. **Sign-recovery readouts (all reported honestly, side by side):**
   - **Reference-based residual recovery.** Subtracted the fixed non-malignant
     reference contribution (Σ frac × mean profile) from each bulk sample and
     divided by the estimated malignant fraction to recover malignant-compartment
     CPM; re-scored MDM2 and IFN/MHC programs on 44 samples with malignant
     fraction ≥ 0.05.
   - **Partial correlation** of the two bulk program scores controlling for all
     five non-malignant fractions.
   - **High-purity subset:** naive correlation restricted to the top-50%
     malignant-fraction samples (n = 27).

## Results

### Fraction estimates are valid

The deconvolution recovers a **low-purity, immune-rich mixture**: mean malignant
fraction **0.37**, mean total immune fraction **0.44** (myeloid-dominated, mean
0.36). Estimates are sane against an orthogonal marker — bulk CD45/PTPRC
expression tracks the estimated immune fraction (**r = +0.62**, p = 7×10⁻⁷) and
anti-tracks the estimated malignant fraction (**r = −0.39**, p = 0.004). The
confound is confirmed directly: the estimated immune fraction correlates with the
naive bulk IFN/MHC program (r = +0.46).

### Deconvolution does NOT rescue the malignant-intrinsic sign

| Readout | r | p | n |
|---|---|---|---|
| scRNA clone-level (ground-truth target) | **−0.630** | — | 55 |
| Naive bulk (the Simpson flip) | +0.282 | 0.041 | 53 |
| **Deconvolution-recovered malignant expression** | **+0.124** | 0.421 | 44 |
| Partial correlation (control non-malig fractions) | +0.315 | 0.021 | 53 |
| High-purity subset (top-50% malignant) | +0.331 | 0.092 | 27 |

**The honest verdict: deconvolution did not recover the negative sign.**
Reference-based residual recovery moved the correlation from **+0.28 → +0.12**
— i.e. it **collapsed the spurious positive toward zero** (the immune confound
is largely removed, and the result is no longer significant), but it did **not**
cross into the malignant-intrinsic negative regime (−0.63). The two
compartment-controlled alternatives (partial correlation, high-purity subset)
stay clearly positive (+0.32, +0.33). No approach approached −0.63.

For comparison, the earlier simple linear immune-adjustment moved r from +0.28
only to +0.16 (n.s.). Full multi-population NNLS deconvolution does modestly
better at *neutralizing* the confound (+0.12, closer to zero) but is no better
at *reversing* it.

## Why the sign does not return — interpretation

1. **Class collinearity / class-1 is nearly invisible in bulk.** The two
   malignant-class signatures are correlated (r = +0.55 in log-CPM), and NNLS
   assigns the malignant compartment almost entirely to class 2 (mean class-2
   share of malignant = 0.95, i.e. class-1 share of malignant ~0.05; class-1
   fraction of total bulk mean 0.02). The bulk simply does
   not contain enough resolvable class-1 (IFN/MHC-high) malignant signal for a
   population-fraction contrast to encode the within-compartment MDM2-vs-IFN
   axis. Consistently, the estimated class-2 fraction does **not** track the bulk
   MDM2 program (r = −0.11).

2. **The −0.63 relationship is a *within-malignant, between-clone* axis, not a
   between-sample mixing axis.** Deconvolution corrects *composition* (how much
   malignant vs. immune), but the clone-level anti-correlation lives in the
   *expression state* of malignant cells. Mean-profile deconvolution assumes one
   fixed malignant expression profile per class; it cannot reconstruct
   sample-specific malignant expression heterogeneity well enough — residual
   recovery amplifies noise (division by a small, uncertain malignant fraction)
   rather than cleanly isolating the intrinsic axis.

3. **Bulk IFN/MHC remains partly immune-driven even after subtraction**, because
   the fixed reference profile cannot capture sample-to-sample variation in
   immune activation state.

## What this means for the project design

- **Composition confounding is real and only partly removable.** Multi-population
  deconvolution successfully *identifies and neutralizes* the immune-infiltration
  confound (positive correlation collapses to non-significant ~0), which is the
  correct, defensible outcome — but it does **not** manufacture the
  malignant-intrinsic signal that whole-bulk mixing destroyed.
- **The malignant-intrinsic MDM2-vs-IFN anti-correlation is a single-cell / clone-level
  observation that does not transfer to these 53 matched bulk samples**, even
  with deconvolution. Any project claim about this axis should be scoped to the
  single-cell compartment and explicitly flagged as **not recovered in matched
  bulk** — reporting +0.12 (deconvolved, n.s.) vs. −0.63 (scRNA) is the honest
  framing.
- **Design implication:** to test malignant-intrinsic RNA-state relationships in
  bulk cohorts, either (a) use purity-stratified / high-malignant-content samples
  with paired estimates, or (b) prefer physical enrichment (LCM, sorted malignant
  cells) over in-silico deconvolution when the signal of interest is a
  within-compartment expression-state axis rather than a composition axis.

## Files

- `ddlps_deconvolution.png` — 3-panel figure (fractions, sign-recovery test, validation)
- `ddlps_deconv_fractions.csv` — per-sample population fractions + program scores + validation markers
- `ddlps_deconv_signature.csv` — 560-gene × 7-population signature matrix (CPM)
- `ddlps_deconv_correlations.csv` — naive vs. deconvolution-adjusted correlation table
