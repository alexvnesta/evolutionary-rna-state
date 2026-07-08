"""Evaluation harness — the bar every RNA feature must clear.

- ``batch_robustness`` — per-feature cross-batch reproducibility verdicts
  (cohort-explained variance eta^2 before/after within-cohort harmonization,
  pooled vs. within-batch AUROC). Separates a feature being a robust *measurement*
  from it *predicting* out-of-cohort.
- ``loco_lopo`` — COMPASS-style leave-one-cohort-out / leave-one-patient-out
  incremental AUROC over the dMMR/TMB/PD-L1+GEP floor, with bootstrap CIs.
"""
