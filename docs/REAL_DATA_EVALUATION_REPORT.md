# Real-Data Evaluation Report — Evolutionary RNA-State & ICB Response

*Everything in this report is computed on REAL pilot RNA-seq, not proxies or
synthetic fixtures. Cerebral Valley hackathon · toward July 13, 2026.*

## Data provenance

The raw-read pipeline has delivered gene-level quantification for **52 samples
across 3 public melanoma ICB cohorts** — Gide 2019 (n=30), Hugo 2016 (n=12),
Riaz 2017 (n=10), all pre-treatment. Response labels come from the harmonized
iAtlas frame. The pipeline keys samples on ENA/SRA run accessions; the labels key
on study sample IDs. We built and **independently validated** the crosswalk: our
ENA-derived mapping agrees 40/40 with the crosswalk a separate pipeline session
produced, and the `run_catalog.csv` path extends it to all 3 cohorts (52/52 mapped,
52/52 labels attached). No sample's response label is guessed or imputed.

## What is computable NOW vs. what is pending

| Layer | Feature | Real-data status |
|---|---|---|
| Baseline floor | T-cell-inflamed GEP, IFN-γ, Teff-TGFβ | **computed on n=52** |
| Differentiated (expression proxy) | splicing-factor / broad-RBP / ADAR-editing activity | **computed on n=52** |
| Differentiated (de-novo antigen) | splice / TE / intron-retention / editing / fusion neoantigen burden | **pending** — matrices not yet emitted |

The de-novo antigen burdens — the actual substrate of the hypothesis — require the
pipeline's splicing-junction, TE-locus, intron-retention, editing-site, fusion and
variant matrices, which have not landed. `pilot_ingest.py` fills those reserved
slots and re-runs evaluation automatically when they do.

## Result 1 — Baseline GEP behaves as a good biomarker (and exposes why LOCO matters)

Real T-cell-inflamed GEP over the runnable floor (TIDE; TMB unavailable for these
RNA cohorts):

| Evaluation | ΔAUROC | 95% CI | boot p |
|---|---|---|---|
| Pooled LOPO | +0.190 | [+0.024, +0.356] | 0.014 |
| LOCO −Gide | +0.259 | [+0.134, +0.384] | 0.000 |
| LOCO −Hugo | −0.306 | [−0.694, +0.083] | 0.933 |
| LOCO −Riaz | +0.080 | [−0.140, +0.280] | 0.231 |
| LOCO mean | +0.094 | — | — |

Within-cohort AUROC: Gide **0.91**, Hugo 0.56, Riaz 0.52. GEP adds signal pooled,
but LOCO shows it **does not transport uniformly** — strong in Gide, near-chance or
negative elsewhere. Even a proven biomarker fails to generalize cleanly across
cohorts here. This is the single most important methodological point: **pooled CV
overstates; LOCO is the honest bar.**

Two distinct things are being measured and must not be conflated: (i) *feature-value
batch-robustness* — is the GEP score itself distorted by which platform/cohort a
sample came from? — and (ii) *predictive transport* — does the score's association
with response generalize to a held-out cohort? On the real n=52 data GEP is robust in
sense (i): after within-cohort harmonization its cohort-explained variance is η²≈0.0
and pooled-vs-within-batch AUROC are 0.751 vs 0.754 (`batch_robustness_real_n52.csv`,
verdict ROBUST) — the score is comparable across batches. But it is **not** uniform in
sense (ii): the within-cohort AUROCs above (0.91/0.56/0.52) and the negative
leave-Hugo-out ΔAUROC show the score→response relationship does not transport. A
feature can be batch-robust as a *measurement* yet fail to *predict* out-of-cohort;
GEP is exactly that case, and it is why we report LOCO, not just batch-robustness.

## Result 2 — Expression-derived RNA-processing activity does not beat the floor

Named regulator-activity features (rbp-activity-scorer: 36 splicing factors, 17
broad RBPs, 3 ADAR enzymes), scored on all 52 samples, tested over the TIDE floor:

| Feature | Pooled ΔAUROC | 95% CI | boot p |
|---|---|---|---|
| Splicing-factor activity | −0.044 | [−0.164, +0.076] | 0.77 |
| Broad-RBP activity | −0.033 | [−0.105, +0.030] | 0.84 |
| ADAR-editing activity | +0.076 | [−0.057, +0.214] | 0.14 |

None beats the floor; every CI crosses zero. ADAR-editing activity has the only
positive point estimate but is not significant. Consistent with the project's prior
finding that summed WES-proxy antigen burdens showed no marginal signal over TMB.

## Result 3 — The RNA-processing machinery co-varies, but that axis does not organize response

Testing the latent-state hypothesis at the level of regulator expression (n=52,
within-cohort standardized):

- **Two coherent co-variation blocks.** Immune signatures cluster (GEP-IFNγ ρ=0.96,
  GEP-Teff ρ=0.57). Separately, **splicing-factor and broad-RBP activity co-vary
  (ρ=0.73)** — a coherent RNA-processing axis. ADAR-editing activity is a largely
  independent third axis.
- **The RNA-processing block is distinct from immune signal** — weakly *anti*-
  correlated with the immune block (ρ = −0.07 to −0.20), so it is not a restatement
  of infiltration.
- **But it does not predict response.** RNA-processing PC1 (57% of axis variance,
  loadings 0.71/0.70 on splicing/RBP) vs. response: direction-agnostic AUROC 0.519,
  **permutation p = 0.82.** Only the immune GEP block separates responders (0.75).

**Reading:** there is measurable, coordinated RNA-processing-machinery activity —
partial support for a shared latent RNA state — but at the level of bulk regulator
expression it carries no ICB-response signal beyond the immune floor, at n=52.
Whether the *de-novo antigen burden* those regulators produce carries signal is the
open question, and the one this pilot cannot yet answer.

## Honest bottom line

On real 3-cohort data we have: a validated crosswalk, a working
crosswalk→score→LOCO machine, a proven-biomarker positive control (GEP) that behaves
correctly and demonstrates the LOCO-vs-pooled lesson, and two honest interim nulls
for expression-derived RNA-processing features. The differentiated de-novo antigen
features — the point of the project — remain pending their pipeline matrices, wired
and ready. Nothing here is fabricated; every number traces to a saved artifact.

---
*Artifacts: results_real_gep_loco_n52.json · results_regulator_activity_loco_n52.json
· results_covariation_real_n52.json · fig_real_gep_n52.png · fig_regulator_activity_n52.png
· fig_covariation_real_n52.png · pilot_n52_crosswalk.parquet · pilot_ingest.py (orchestrator).*
