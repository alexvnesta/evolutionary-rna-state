# Hypothesis kill-tests (analysis side)

Tests of the core hypothesis — "bulk RNA-seq reconstructs a latent evolutionary state that
determines antigenicity and ICB response" — run on the iAtlas melanoma cohorts. Top-level
read: `docs/FINAL_VERDICT.md`. Master figure: `docs/figures/fig_master_verdict.png`.

## Files
- `move2_autorun.py` (v7) — idempotent auto-runner: detect grown expression matrix → score
  regulator activity → build/validate per-cohort crosswalk → run all test families → save.
  Fires four families: T1/T2/T3 activity, reframing P1-P3, learned-representation head-to-head,
  and per-cohort ID-mapping diagnostics. Wired for liu2019 (passthrough maps sampleId OR
  patientId → 122/122). Reproduce in any session: `run_if_changed(host)`.
- `headtohead_repr.py` — the untested-core experiment: unsupervised RNA embedding (PCA/kPCA/AE,
  fit inside CV) vs immune floor, under random-5fold + leave-one-cohort-out with permutation null.
- `reframing_tests.py` — P1 coordination / P2 immune-cold coupling / P3 IFN-incremental.
- `run_activity_response_tests.py` — T1 mechanistic / T2 shared-state / T3 LOCO-AUROC.
- `rbp_activity_scorer.py`, `move2_watch_live.py` — regulator scorer + in-session poll loop.
- `outputs/` — all JSON/CSV outputs (40-sample 2-cohort and 66-sample 3-cohort).

## Result (current data: 66 samples, gide+hugo+riaz)
Learned RNA representation is at CHANCE across cohorts (LOCO AUROC 0.49, perm p=0.41); immune
floor real-but-weak (0.59, p=0.044); combining degrades the floor. Tumor-intrinsic RNA adds
nothing over immune composition (residual p=0.38). The hypothesis's antigen-quantity core is
dead (416 samples) and its learned-representation method fails cross-cohort transfer. What
survives: immune-composition dominance + a DNA-side heterogeneity effect.

liu2019 expression (n~122, has TMB) is the powered adjudication — 4th held-out cohort, ~n180.
The auto-runner fires it automatically when liu lands.
