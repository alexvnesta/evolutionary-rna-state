# Pre-registered evaluation protocol — non-reference RNA features vs the immune floor

**Locked 2026-07-10 by session `23cf8106`, BEFORE any result was computed.** No decision below may be
changed after seeing an AUROC. Deviations, if forced, are logged in a dated "Deviations" section with
reason. This is the anti-tuning contract.

## Question
Do the hypothesis's **non-reference RNA-processing features** (RNA editing, intron retention, splicing,
TE/ERV family activity) carry ICB-response signal **beyond** what the immune floor (an IFN-γ / T-cell
expression signature) already provides?

This is the first test ever run on the actual non-reference feature block — every prior "DEAD" verdict
used expression, expression-PCA, or a DNA/WES antigen-quantity proxy (`AUDIT_SYNTHESIS.md`).

## Data & blocks (fixed)
- **Label:** `y` (responder=1/non-responder=0), pre-treatment only, from `reconciled_frame_n106.parquet`.
- **Immune floor block (5 feat):** `gep_tcell_inflamed, ifng_score, teff, tgfb, teff_tgfb_balance`
  (`immune_floor_block.parquet`). This is the established positive control (prior LOCO 0.767, perm p=0.001).
- **Non-reference block (≤66 feat):** editing (AEI %, S/N, A→G, A-cov), intron retention (mean/median/
  frac>0.1/n_eval), splicing (junction burden), TE/ERV family counts (~57 families). Built by the unified
  pipeline / `results/nonref_run/` callers, assembled into one parquet keyed on `run_accession`.
- **Sample set:** whatever samples have BOTH blocks. Recorded honestly with an n and per-cohort breakdown.

## Primary analysis (fixed)
1. **Two blocks, three models per split:** (A) floor only, (B) non-ref only, (C) floor + non-ref.
   Classifier: L2-penalized logistic regression (same family as the established floor harness).
2. **Cross-validation frame — decided by cohort availability, declared in advance:**
   - **If ≥2 cohorts each with ≥10 samples and both classes:** **leave-one-cohort-out (LOCO)** is the
     primary. This is the real cross-cohort transfer test.
   - **If effectively one cohort (the honest current state: ~34 Gide, ~2 Riaz, 0 Hugo):** primary is
     **grouped 5-fold CV within Gide, grouped by `patient_id`** (never by sample). LOCO is reported as
     "pending Hugo/Riaz alignment" and is the headline ONLY when the ≥2-cohort condition is met.
   - Hugo alignment is running in the background specifically to move this from the second case to the first.
3. **Metric:** ROC AUROC, out-of-fold. Report mean + Hanley-McNeil 95% CI.
4. **Incremental test:** does block C beat block A? Report ΔAUROC (C − A) and a DeLong / bootstrap CI on it.
5. **Null:** 5000-permutation label shuffle (within-cohort exchangeability), report perm p for each block.
6. **Success criterion (stated up front, not tunable):** the non-reference layer is judged to carry
   independent signal ONLY if block B beats its permutation null (perm p < 0.05) AND block C's ΔAUROC over
   block A has a CI excluding 0. Anything weaker is reported as "no independent non-ref signal at this n" —
   a legitimate, publishable negative, given the honest prior (ADAR editing collapsed to 0.535 at n=106).

## Leakage discipline (fixed — signature-rigour test 3)
- Every data-dependent transform (imputation, `StandardScaler`, any feature selection) is fit **inside the
  training fold only**. No global standardization before CV.
- Folds grouped by **`patient_id`**, not `run_accession` — a patient with multiple runs never straddles
  the split.
- TE-family features are high-dimensional (~57) at low n: L2 regularization + in-fold selection only; no
  peeking at the full matrix for feature choice.

## Circularity check (fixed — test 1)
- The non-ref block shares **no genes** with the immune floor by construction (editing indices, intron/exon
  read ratios, TE-family counts vs IFN-γ pathway genes) — report the explicit zero-overlap.
- Size-matched **random-feature permutation null**: draw random same-size feature blocks from the expression
  matrix; if they match the non-ref block's incremental AUROC, the non-ref block is not special. Report.

## Composition confounding check (fixed — test 2)
- TE/ERV and editing signals can rise from infiltrate rather than tumor. Re-fit block C with **tumor purity**
  (InstaPrism malignant fraction, already computed n=106, `mediation_frame_purity_n106.parquet`) as a
  covariate. A non-ref signal that vanishes under purity adjustment is flagged as composition-driven, not
  tumor-intrinsic.

## Provenance check (fixed — test 4)
- Cohorts resolved by accession (Gide ERR2208*, Hugo SRR318*, Riaz SRR508*) via the committed crosswalks,
  not title guesses. Caller tool versions (samtools, subread/featureCounts, regtools, hisat2) recorded from
  the `nonref-callers`/`editing` conda envs into the result manifest.

## What a NEGATIVE result means (fixed)
If the non-ref block shows no independent signal: this is the **first honest test** of the hypothesis's own
features, not a proxy — a real, reportable negative at the achieved n. It does NOT retroactively validate the
earlier proxy-based "DEAD" verdicts, and it does NOT test the learned-representation half (scoped separately).

## Deviations (logged, dated)
- **2026-07-10, after first run:** the permutation null was initially coded to shuffle labels within the
  CV grouping variable (`patient_id`). In the within-Gide frame patients are singletons, so within-group
  shuffling is a no-op → degenerate `perm_p=1.0` for every block. **Fix:** the permutation's exchangeable
  unit is now the **cohort** (`perm_block=cohort`), decoupled from the CV split grouping. This is a
  correctness fix to a degenerate null, not a change to the models, features, CV frame, or success
  criterion — the point AUROCs and CIs are unaffected; only the perm-p values change from the degenerate
  1.0 to their true values. Logged here per the pre-registration contract.

## Outputs (fixed)
`results/eval/nonref_vs_floor_loco.json` (or `_gidecv.json`), `two_block_result.png`,
and a row appended to the result brief. Analysis script: `analysis/two_block_eval.py` (committed with this
protocol, run only AFTER the matrix exists).
