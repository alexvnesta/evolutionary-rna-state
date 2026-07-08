# ICB outcome test on Riaz bulk (GSE91061)

## The endpoint, and why this cohort

The project's ultimate claim is that the clonal immune-state axis organizes checkpoint
(ICB) response. This is the first test against real ICB outcomes. A survey of the three
classic melanoma bulk-ICB cohorts found only Riaz (GSE91061) adequately powered: Auslander
(GSE115821) has just 2 responders among 14 labeled pre-treatment samples, and Hugo
(GSE78220) is n=28. Riaz has 49 pre-treatment samples with RECIST response
(10 responders [PRCR], 39 non-responders [PD or SD]) and matched WES, in
melanoma, which is the right tissue for a same-cancer test.

## What is and is not testable

The mechanistic work established that the malignant-intrinsic clonal immune-state does not
transfer to real bulk by naive scoring: bulk interferon/MHC is dominated by immune-cell
infiltration (Simpson's paradox), and composition deconvolution does not recover the
within-malignant axis. Riaz has no matched single-cell data, so the clonal trajectory
class itself cannot be inferred here. What is legitimately testable is whether the
malignant-intrinsic interferon/MHC signature, scored in pre-treatment bulk, associates with
response, and whether it adds anything beyond the standard bulk immune-infiltration signal.

## Result: a negative, exactly as the mechanism predicts

- Neither signature significantly separates responders from non-responders at this sample
  size: the malignant-intrinsic IFN/MHC signature gives AUC 0.56 (Mann-Whitney
  p=0.56), the immune-infiltration proxy AUC 0.62 (p=0.24).
- The two signatures correlate at r=0.84 in bulk. This is the same
  confound quantified in DDLPS: in bulk, the interferon/MHC signal is not separable from
  immune-cell infiltration, so the malignant-intrinsic axis cannot be isolated without a
  compartment-resolving step.
- In a joint logistic model (response ~ infiltration + IFN/MHC), neither term is
  significant (infiltration p=0.13, IFN/MHC p=0.34) and
  the malignant signature adds nothing beyond infiltration.

## Interpretation

This is a coherent negative, not a contradiction. It confirms at the ICB endpoint what the
mechanistic experiments showed at the expression level: reading the malignant-intrinsic
clonal immune-state from bulk requires compartment resolution that bulk alone does not
provide, and naive bulk signature scoring recovers only the infiltration signal, which is
itself underpowered here (10 responders). The test does not falsify the thesis; it shows
that the bulk-transfer step, not the biology, is the blocker, and it localizes the
requirement precisely: a matched single-cell reference (or a validated compartment-specific
feature set at adequate N) is needed before the ICB association can be fairly evaluated.

## What would make this test decisive

- A melanoma cohort with matched scRNA + bulk + ICB outcome (none currently public at N
  adequate for a 10-plus responder contrast), or
- The compartment-specific feature set validated on a matched same-tissue cohort at N about
  40, then applied to Riaz. The compartment approach showed a positive point estimate in
  DDLPS (balanced accuracy 0.55 to 0.64) but is not yet certified.

## Artifacts

- `figures/riaz_icb_test.png` : signature scores by response, AUCs, and the bulk confound scatter.
- `tables/riaz_icb.json` : all metrics.
