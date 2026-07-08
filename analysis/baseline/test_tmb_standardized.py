#!/usr/bin/env python
"""
test_tmb_standardized.py — unit validation for the standardized-TMB module.

Two kinds of check:
  1. Numerical properties of the transform on a small controlled frame:
       * per-Mb rate reproduces a known callable size;
       * harmonization collapses a synthetic per-batch offset;
       * the transform is rank-preserving WITHIN each batch (Spearman ~ 1.0);
       * a fitted alignment applied to a held-out sample lands near the reference.
  2. The banked-data assertion (only if analysis_frame.parquet is present):
       harmonization reduces between-cohort eta^2 and keeps within-cohort rho ~1.

Run:  python analysis/baseline/test_tmb_standardized.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tmb_standardized as ts  # noqa: E402


def _synthetic_frame(seed: int = 0) -> pd.DataFrame:
    """Two cohorts with the SAME latent TMB distribution but a 3x assay offset
    (cohort B calls 3x more mutations on identical biology)."""
    rng = np.random.default_rng(seed)
    latent = rng.lognormal(mean=1.5, sigma=0.8, size=120)   # muts/Mb, shared truth
    mb = 30.0
    rows = []
    for i, lat in enumerate(latent[:60]):
        rows.append(dict(run_accession=f"A{i}", cohort="A",
                         TMB_NONSYNONYMOUS=lat, MUTATION_COUNT=lat * mb))
    for i, lat in enumerate(latent[60:]):
        infl = lat * 3.0                                     # assay inflation
        rows.append(dict(run_accession=f"B{i}", cohort="B",
                         TMB_NONSYNONYMOUS=infl, MUTATION_COUNT=infl * mb))
    return pd.DataFrame(rows)


def test_rate_and_callable_inference():
    n = np.array([300.0, 600.0])
    mb = np.array([30.0, 30.0])
    assert np.allclose(ts.per_mb_rate(n, mb), [10.0, 20.0])
    # infer callable Mb back from (count, rate)
    inferred = ts.infer_callable_mb(n, np.array([10.0, 20.0]))
    assert np.allclose(inferred, [30.0, 30.0])
    # zero callable -> NaN rate
    assert np.isnan(ts.per_mb_rate([100.0], [0.0])[0])
    print("  rate + callable inference: OK")


def test_harmonization_collapses_offset():
    f = _synthetic_frame()
    out, fit = ts.standardize_tmb(f, batch_col="cohort", method="center")
    med = out.groupby("cohort")["tmb_standardized"].median()
    # the 3x offset should be largely removed: medians within ~15%
    ratio = med["B"] / med["A"]
    print(f"  post-harmonization median ratio B/A = {ratio:.3f} (raw was ~3.0)")
    assert 0.85 <= ratio <= 1.18, f"offset not collapsed (ratio {ratio:.3f})"

    # between-batch eta^2 must drop sharply
    log_raw = np.log10(out["tmb_rate"] + ts.LOG_PSEUDOCOUNT)
    log_std = np.log10(out["tmb_standardized"] + ts.LOG_PSEUDOCOUNT)
    e_before = ts.eta_squared(log_raw, out["cohort"])
    e_after = ts.eta_squared(log_std, out["cohort"])
    print(f"  eta^2: {e_before:.3f} -> {e_after:.3f}")
    assert e_after < e_before / 3, "harmonization did not reduce batch separation"


def test_rank_preserving_within_batch():
    f = _synthetic_frame()
    out, _ = ts.standardize_tmb(f, batch_col="cohort", method="center")
    rho = ts.rank_preservation(out["tmb_rate"], out["tmb_standardized"], out["cohort"])
    print(f"  within-batch Spearman rho = {rho:.4f}")
    assert rho > 0.999, "location shift must preserve within-batch ranking"


def test_fit_apply_holdout():
    f = _synthetic_frame()
    _, fit = ts.standardize_tmb(f, batch_col="cohort", method="center")
    # apply the FITTED alignment to a new cohort-B sample at B's median
    med_B_rate = f.loc[f.cohort == "B", "TMB_NONSYNONYMOUS"].median()
    val = fit.apply([med_B_rate], ["B"])[0]
    ref = 10 ** fit.ref_loc - fit.pseudocount
    print(f"  held-out B median rate {med_B_rate:.2f} -> standardized {val:.2f} "
          f"(ref {ref:.2f})")
    assert abs(np.log10(val + 1e-3) - fit.ref_loc) < 0.05


def test_center_scale_equalizes_spread():
    f = _synthetic_frame()
    out, _ = ts.standardize_tmb(f, batch_col="cohort", method="center_scale")
    log_std = np.log10(out["tmb_standardized"] + ts.LOG_PSEUDOCOUNT)
    mad = out.assign(l=log_std).groupby("cohort")["l"].apply(
        lambda x: np.median(np.abs(x - np.median(x))))
    print(f"  per-cohort MAD after center_scale: {mad.round(3).to_dict()}")
    assert abs(mad["A"] - mad["B"]) < 0.15, "center_scale did not equalize spread"


def test_banked_data_if_present():
    candidates = [
        Path(__file__).resolve().parents[2] / "results" / "features" / "_analysis_frame_tmp.parquet",
    ]
    frame_path = next((p for p in candidates if p.exists()), None)
    if frame_path is None:
        print("  banked-data check: SKIPPED (analysis_frame not on disk)")
        return
    frame = pd.read_parquet(frame_path)
    frame = frame[frame["TMB_NONSYNONYMOUS"].notna()].copy()
    out, _ = ts.standardize_tmb(frame, batch_col="cohort", method="center")
    log_raw = np.log10(out["tmb_rate"] + ts.LOG_PSEUDOCOUNT)
    log_std = np.log10(out["tmb_standardized"] + ts.LOG_PSEUDOCOUNT)
    e_before = ts.eta_squared(log_raw, out["cohort"])
    e_after = ts.eta_squared(log_std, out["cohort"])
    rho = ts.rank_preservation(out["tmb_rate"], out["tmb_standardized"], out["cohort"])
    print(f"  banked: eta^2 {e_before:.4f} -> {e_after:.4f}, within-cohort rho {rho:.4f}")
    assert e_after < e_before, "no cross-cohort improvement on banked data"
    assert rho > 0.999


if __name__ == "__main__":
    print("Standardized-TMB module validation:")
    test_rate_and_callable_inference()
    test_harmonization_collapses_offset()
    test_rank_preserving_within_batch()
    test_fit_apply_holdout()
    test_center_scale_equalizes_spread()
    test_banked_data_if_present()
    print("\nALL TMB STANDARDIZATION TESTS PASSED")
