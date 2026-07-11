# Why the non-reference sequence layer never adds signal — the aberrancy×dosage insight

**Date:** 2026-07-11  **Status:** COMPLETE  **Compute:** local CPU only (no new Modal spend)

## The question this answers
Every non-reference RNA layer we tested — editing, intron retention, TE/ERV, splicing burden, and the
Evo2 novel-junction *aberrancy* score — came back null against the immune floor. This analysis asks the
mechanistically-motivated question the hypothesis actually rests on, and finds *why* they all fail.

## The hypothesis's real claim
Not "aberrant junctions predict response," but that a junction which is BOTH sequence-surprising
(high Evo2 aberrancy) AND abundant (high read support) is the aberrant-and-EXPRESSED event that could
generate neoantigens. That is an **interaction**, not a marginal aberrancy score. We built it explicitly:
per sample, aberrancy(−delta) × log10(reads) summaries over the top-200 novel junctions — the
"aberrant-expressed burden."

## Result — the interaction is also null, but informatively so
Two-block test (20-seed 5-fold CV, fold-contained residualization, cohort-internal permutation):

| | Gide n=32 | Hugo n=22 |
|---|---|---|
| Immune floor | 0.792 | 0.593 |
| Interaction alone | **0.612** | 0.383 |
| Floor + interaction | 0.763 | 0.535 |
| Interaction residualized on floor | 0.421 | 0.528 |
| permutation p (residual) | **0.92** | 0.219 |

The interaction alone beats chance in Gide (0.612 > aberrancy-alone's 0.513), yet collapses to 0.421 once
residualized on the floor. That specific pattern — predictive alone, nothing left after removing the floor —
is the signature of a **proxy**, not an independent signal.

## The mechanism (the actual insight)
The aberrant-expressed junction burden is a **correlate of the inflamed-tumor state**, not an independent
evolutionary readout:
- In Gide it correlates positively with every immune-infiltration axis (GEP +0.36, IFNγ +0.29, Teff +0.29)
  and trends with response (frac-aberrant-abundant vs y: ρ=+0.33, p=0.064).
- **Its coupling to inflammation FLIPS sign across cohorts** — positive in Gide, negative in Hugo. Same-feature
  GEP correlation: the aberrant-expressed *burden* (ix_sum_aberr×log-reads) goes +0.36 (Gide) → −0.15 (Hugo);
  the aberrant-*and-abundant fraction* goes +0.31 (Gide) → −0.37 (Hugo). The flip is systematic across all five
  inflammation axes for both features (see figure panels B and C).

So the measurable non-reference "aberrancy burden" is entangled with immune composition and inherits the
**exact same cross-cohort direction flip** that makes the immune floor itself fail to transfer (Phase 2).
This is a single, unifying explanation for why NONE of the non-reference layers — interpretable callers or
Evo2 sequence scores — ever add independent, portable ICB signal at this scale: they are downstream
readouts of the tumor-immune state, measured through a lens that reorients between cohorts.

## What would be needed to salvage the hypothesis's independent-signal claim
- A feature that is aberrancy-driven but **decoupled** from bulk immune infiltration (e.g. per-neojunction
  HLA-presented peptide load via the presentation layer, not a bulk burden) — the ENCODER_PHASE_PROTOCOL
  presentation layer.
- Cohorts where the immune floor itself transfers (so residual signal is interpretable), i.e. a single
  large cohort with internal splits, not 3 small flipping cohorts.

Figure: `evo2_interaction_insight.png`. Data: `results/eval/evo2_interaction_test.json`.
