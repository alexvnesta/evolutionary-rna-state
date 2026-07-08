#!/usr/bin/env python
"""
analysis/baseline/tmb_standardized.py

BASELINE FLOOR — Standardized tumor mutational burden (`tmb_standardized`).

WHY THIS EXISTS
---------------
Raw TMB (nonsynonymous mutations per Mb) is the most-cited ICB biomarker, but
its cross-study comparability collapsed once people pooled panels/platforms:
different assays call different numbers of mutations on the *same* tumour
because of panel size, callable-territory definition, VAF/depth cutoffs, and
germline-filtering strategy. The Friends-of-Cancer-Research TMB Harmonization
Project (Merino et al., J Immunother Cancer 2020;8:e000147) showed that naive
per-Mb TMB is NOT interchangeable across assays and prescribed *aligning each
assay's TMB distribution to a reference* before any cross-cohort use.

This module reproduces that floor with an explicit, documented two-stage
transform, and — critically — a **batch/panel covariate hook** so the
harmonization target (`batch`) can be a cohort, a sequencing platform, or a
capture panel, whichever is the real source of the assay effect.

    stage 1 (rate):        muts / callable_Mb                  -> per-Mb rate
    stage 2 (harmonize):   align per-batch log10(rate) to a    -> tmb_standardized
                           common reference (location, or
                           location+scale)

Optionally a stage-3 within-batch z-score (`tmb_standardized_z`) is emitted for
models that pool cohorts, exactly as the FEATURE_CONTRACT_v2 batch-robustness
note prescribes ("z-score within cohort before pooling; never pool raw counts").

WHAT MAKES IT BATCH-ROBUST
--------------------------
* Denominator is *callable* Mb, not a fixed 38 Mb exome guess — supplied per
  sample (WES BED callable size) or inferred from a paired raw count + per-Mb
  rate when only banked TMB is available.
* Stage-2 harmonization removes the per-assay offset (and optionally scale)
  that the FoCR project identified as the dominant non-biological variance,
  while a monotone (rank-preserving) transform keeps within-cohort ordering
  intact — so a sample that was high-TMB in its own cohort stays high.
* The `batch` argument is the covariate hook: pass `cohort` today; pass
  `platform`/`panel` when that metadata lands, no code change.

VALIDATED ON BANKED DATA
------------------------
`main()` runs on results analysis_frame.parquet (TMB_NONSYNONYMOUS +
MUTATION_COUNT, 4 melanoma ICB cohorts) and reports between-cohort variance
(eta^2 of log-TMB) before vs after harmonization and within-cohort Spearman
rho (must stay ~1.0). This is runnable now — it is the one baseline feature we
can reproduce on banked data without new sequencing.

Feature column produced: `tmb_standardized` (float, muts/Mb, harmonized).
Also emits `tmb_standardized_z` (float, within-batch z) and `callable_mb`
(float, QC) when the inputs allow.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
LOG_PSEUDOCOUNT = 1e-3     # added before log10 so rate==0 is finite & monotone
DEFAULT_EXOME_MB = 38.0    # fallback callable size (WES capture, if nothing else)
KEY_COLS = ["run_accession", "cohort"]

# Feature column names this module writes (see FEATURE_CONTRACT_v2.md)
FEATURE_COL = "tmb_standardized"
FEATURE_COL_Z = "tmb_standardized_z"


# ---------------------------------------------------------------------------
# Stage 1 — per-Mb rate over callable territory
# ---------------------------------------------------------------------------
def per_mb_rate(
    n_nonsyn: Sequence[float],
    callable_mb: Sequence[float],
) -> np.ndarray:
    """Nonsynonymous coding mutations per Mb of *callable* territory.

    callable_mb is the size (Mb) of the region over which mutations could be
    called for that sample — the FoCR-critical denominator. Zero/NaN callable
    sizes yield NaN (cannot standardize what was not callable).
    """
    n = np.asarray(n_nonsyn, float)
    mb = np.asarray(callable_mb, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        rate = np.where((mb > 0) & np.isfinite(mb), n / mb, np.nan)
    return rate


def infer_callable_mb(
    n_nonsyn: Sequence[float],
    reported_rate: Sequence[float],
    fallback_mb: float = DEFAULT_EXOME_MB,
) -> np.ndarray:
    """Recover per-sample callable Mb from a banked (count, per-Mb rate) pair.

    When a cohort banked both a raw nonsynonymous count and a per-Mb TMB, the
    implied callable size is count / rate. This lets us standardize banked data
    whose BED files we do not have. Falls back to `fallback_mb` when the rate is
    missing or zero.
    """
    n = np.asarray(n_nonsyn, float)
    r = np.asarray(reported_rate, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        mb = np.where((r > 0) & np.isfinite(r), n / r, np.nan)
    mb = np.where(np.isfinite(mb) & (mb > 0), mb, fallback_mb)
    return mb


# ---------------------------------------------------------------------------
# Stage 2 — cross-batch harmonization (FoCR / Merino 2020 style)
# ---------------------------------------------------------------------------
@dataclass
class HarmonizationFit:
    """The learned per-batch alignment — persist it to apply to future samples.

    Applying a *fitted* alignment (rather than re-deriving offsets from a new
    batch's own distribution) is what makes a held-out sample's standardized
    TMB comparable to the training cohorts.
    """
    method: str                       # "center" | "center_scale"
    ref_loc: float
    ref_scale: float
    batch_loc: dict                   # batch -> median(log10 rate)
    batch_scale: dict                 # batch -> MAD(log10 rate)*1.4826
    pseudocount: float

    def apply(self, rate: Sequence[float], batch: Sequence[str]) -> np.ndarray:
        rate = np.asarray(rate, float)
        batch = np.asarray(batch)
        x = np.log10(rate + self.pseudocount)
        out = x.copy()
        for b in np.unique(batch):
            m = batch == b
            loc = self.batch_loc.get(b, self.ref_loc)
            if self.method == "center":
                out[m] = x[m] - loc + self.ref_loc
            elif self.method == "center_scale":
                scale = self.batch_scale.get(b, self.ref_scale) or 1.0
                out[m] = (x[m] - loc) / scale * self.ref_scale + self.ref_loc
        return np.where(np.isfinite(x), 10 ** out - self.pseudocount, np.nan)


def fit_harmonization(
    rate: Sequence[float],
    batch: Sequence[str],
    method: str = "center",
    ref_batch: str | None = None,
    pseudocount: float = LOG_PSEUDOCOUNT,
) -> HarmonizationFit:
    """Learn the per-batch alignment of log10(per-Mb rate) to a reference.

    method="center"        : subtract each batch's median log-rate offset (the
                             FoCR-observed dominant assay effect); rank-preserving.
    method="center_scale"  : additionally equalize each batch's spread (MAD) to
                             the reference — use when assays differ in dynamic
                             range, not just offset.
    ref_batch=None         : reference = pooled median/MAD across all samples;
                             otherwise the named batch is the anchor.

    Robust statistics (median, MAD) are used so a few hypermutators do not move
    the alignment. NaN rates are ignored in fitting and stay NaN on apply.
    """
    if method not in ("center", "center_scale"):
        raise ValueError(f"unknown method {method!r}")
    rate = np.asarray(rate, float)
    batch = np.asarray(batch)
    x = np.log10(rate + pseudocount)
    ok = np.isfinite(x)

    def _loc_scale(vals):
        loc = float(np.median(vals))
        scale = float(np.median(np.abs(vals - loc)) * 1.4826)
        return loc, (scale if scale > 1e-9 else 1.0)

    if ref_batch is None:
        ref_loc, ref_scale = _loc_scale(x[ok])
    else:
        sel = ok & (batch == ref_batch)
        if not sel.any():
            raise ValueError(f"ref_batch {ref_batch!r} has no finite rates")
        ref_loc, ref_scale = _loc_scale(x[sel])

    batch_loc, batch_scale = {}, {}
    for b in np.unique(batch):
        sel = ok & (batch == b)
        if sel.any():
            loc, scale = _loc_scale(x[sel])
        else:
            loc, scale = ref_loc, ref_scale
        batch_loc[b], batch_scale[b] = loc, scale

    return HarmonizationFit(method, ref_loc, ref_scale,
                            batch_loc, batch_scale, pseudocount)


# ---------------------------------------------------------------------------
# Stage 3 — within-batch z-score (optional pooling helper)
# ---------------------------------------------------------------------------
def within_batch_z(values: Sequence[float], batch: Sequence[str]) -> np.ndarray:
    """Within-batch z-score of log10 values (robust: median/MAD)."""
    v = np.asarray(values, float)
    batch = np.asarray(batch)
    x = np.log10(v + LOG_PSEUDOCOUNT)
    out = np.full_like(x, np.nan, dtype=float)
    for b in np.unique(batch):
        m = batch == b
        xb = x[m]
        ok = np.isfinite(xb)
        if ok.sum() < 2:
            continue
        loc = np.median(xb[ok])
        scale = np.median(np.abs(xb[ok] - loc)) * 1.4826 or 1.0
        out[m] = (xb - loc) / scale
    return out


# ---------------------------------------------------------------------------
# Top-level entry — build the tmb_standardized feature column
# ---------------------------------------------------------------------------
def standardize_tmb(
    frame: pd.DataFrame,
    *,
    batch_col: str = "cohort",
    nonsyn_col: str = "TMB_NONSYNONYMOUS",
    count_col: str | None = "MUTATION_COUNT",
    callable_mb_col: str | None = None,
    method: str = "center",
    ref_batch: str | None = None,
    emit_z: bool = True,
) -> tuple[pd.DataFrame, HarmonizationFit]:
    """Compute `tmb_standardized` (and `tmb_standardized_z`) for a sample frame.

    Parameters
    ----------
    frame        : per-sample table (keyed on run_accession + cohort ideally).
    batch_col    : THE COVARIATE HOOK. 'cohort' today; set to 'platform' or
                   'panel' when that column exists — the harmonization then
                   aligns the true assay batches.
    nonsyn_col   : per-Mb nonsynonymous TMB, OR (if callable_mb given) the raw
                   nonsynonymous *count*. See count/callable logic below.
    count_col    : raw nonsynonymous mutation count (used with a banked per-Mb
                   rate to infer callable Mb). Optional.
    callable_mb_col : per-sample callable size in Mb (WES BED). If given, the
                   rate is count/callable_mb directly; otherwise inferred.
    method       : 'center' (offset only) or 'center_scale'.
    ref_batch    : anchor batch, or None for the pooled reference.

    Returns
    -------
    (out, fit) : `out` has columns tmb_standardized [+ tmb_standardized_z,
                 callable_mb, tmb_rate]; `fit` is the persistable alignment.
    """
    f = frame.copy()

    # --- resolve the per-Mb rate + callable Mb --------------------------------
    if callable_mb_col is not None and callable_mb_col in f:
        # count / callable_mb is the cleanest path (real BED sizes)
        n = f[count_col] if (count_col and count_col in f) else f[nonsyn_col]
        callable_mb = f[callable_mb_col].to_numpy(float)
        rate = per_mb_rate(n, callable_mb)
    elif count_col is not None and count_col in f:
        # banked (count, per-Mb rate) -> infer callable Mb, then rate=count/Mb
        callable_mb = infer_callable_mb(f[count_col], f[nonsyn_col])
        rate = per_mb_rate(f[count_col], callable_mb)
    else:
        # only a per-Mb rate available: standardize it directly on default Mb
        rate = f[nonsyn_col].to_numpy(float)
        callable_mb = np.full(len(f), DEFAULT_EXOME_MB)

    f["callable_mb"] = callable_mb
    f["tmb_rate"] = rate

    # --- harmonize ------------------------------------------------------------
    fit = fit_harmonization(rate, f[batch_col], method=method, ref_batch=ref_batch)
    f[FEATURE_COL] = fit.apply(rate, f[batch_col])

    if emit_z:
        f[FEATURE_COL_Z] = within_batch_z(f[FEATURE_COL], f[batch_col])

    return f, fit


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------
def eta_squared(values: Sequence[float], groups: Sequence[str]) -> float:
    """Between-group variance fraction (0..1) — cross-batch separation metric."""
    v = np.asarray(values, float)
    g = np.asarray(groups)
    ok = np.isfinite(v)
    v, g = v[ok], g[ok]
    if len(v) == 0:
        return np.nan
    grand = v.mean()
    ss_tot = ((v - grand) ** 2).sum()
    if ss_tot == 0:
        return 0.0
    ss_bet = sum(len(v[g == u]) * (v[g == u].mean() - grand) ** 2
                 for u in np.unique(g))
    return float(ss_bet / ss_tot)


def rank_preservation(rate_before, values_after, batch) -> float:
    """Mean within-batch Spearman rho between raw rate and standardized value.

    ~1.0 confirms harmonization is monotone within each batch (no re-ordering
    of a cohort's own samples), which is the safety property of a location/
    location-scale transform.
    """
    from scipy.stats import spearmanr
    rate_before = np.asarray(rate_before, float)
    values_after = np.asarray(values_after, float)
    batch = np.asarray(batch)
    rhos = []
    for b in np.unique(batch):
        m = batch == b
        ok = np.isfinite(rate_before[m]) & np.isfinite(values_after[m])
        if ok.sum() >= 3:
            rhos.append(spearmanr(rate_before[m][ok], values_after[m][ok]).statistic)
    return float(np.nanmean(rhos)) if rhos else np.nan


# ---------------------------------------------------------------------------
# Runnable demo on banked data
# ---------------------------------------------------------------------------
def main(argv: Sequence[str] | None = None) -> int:
    import argparse
    p = argparse.ArgumentParser(description="Standardize TMB on banked cohort data.")
    p.add_argument("--frame", required=True, help="analysis_frame.parquet")
    p.add_argument("--batch-col", default="cohort")
    p.add_argument("--method", default="center", choices=["center", "center_scale"])
    p.add_argument("--out", default=None, help="optional output parquet path")
    args = p.parse_args(argv)

    frame = pd.read_parquet(args.frame)
    frame = frame[frame["TMB_NONSYNONYMOUS"].notna()].copy()

    out, fit = standardize_tmb(frame, batch_col=args.batch_col, method=args.method)

    log_raw = np.log10(out["tmb_rate"] + LOG_PSEUDOCOUNT)
    log_std = np.log10(out[FEATURE_COL] + LOG_PSEUDOCOUNT)
    eta_before = eta_squared(log_raw, out[args.batch_col])
    eta_after = eta_squared(log_std, out[args.batch_col])
    rho = rank_preservation(out["tmb_rate"], out[FEATURE_COL], out[args.batch_col])

    print("Standardized-TMB harmonization report")
    print("=" * 50)
    med_before = out.groupby(args.batch_col)["tmb_rate"].median().round(2)
    med_after = out.groupby(args.batch_col)[FEATURE_COL].median().round(2)
    print("per-batch median muts/Mb  (raw -> standardized):")
    for b in med_before.index:
        print(f"  {b:12s} {med_before[b]:8.2f} -> {med_after[b]:8.2f}")
    print(f"\nbetween-batch eta^2 (log): {eta_before:.4f} -> {eta_after:.4f} "
          f"({eta_before/eta_after:.1f}x reduction)")
    print(f"within-batch Spearman rho (raw vs standardized): {rho:.4f}")
    print("\nfeature columns:", FEATURE_COL, ",", FEATURE_COL_Z)

    if args.out:
        keep = [c for c in KEY_COLS if c in out] + [
            "callable_mb", "tmb_rate", FEATURE_COL, FEATURE_COL_Z]
        out[keep].to_parquet(args.out, index=False)
        print("wrote", args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
