"""batch_robustness.py

Shared cross-platform / cross-batch reproducibility harness for the
evolutionary-RNA-state feature set.

Motivation
----------
Naive tumor mutational burden (TMB) is *not* comparable across cohorts because
each cohort used a different assay / callable-region size: pooling raw counts
mixes a real biological signal with a technical batch axis. That harmonization
failure is the reference failure mode this harness is built to catch for EVERY
feature, before any feature is allowed into the response model.

Given any per-sample feature matrix with a batch label (cohort / platform /
panel), for each feature the harness:

  1. Quantifies the batch effect BEFORE any correction:
       - eta^2  : fraction of total variance explained by batch
                  (one-way ANOVA sum-of-squares decomposition; the effect size).
       - Kruskal-Wallis H test (non-parametric): p-value that >=1 batch differs.
  2. Applies a harmonization transform and RE-quantifies the residual batch
     effect:
       - "zscore"  : within-batch z-score (center+scale per batch)  [default]
       - "rank"    : within-batch rank -> uniform (scale-free, monotone-robust)
       - "combat"  : parametric empirical-Bayes location/scale adjustment
                     (Johnson, Li & Rabinovic 2007, Biostatistics 8:118-127),
                     self-contained reference implementation below. Use when a
                     biological covariate must be PRESERVED while batch is
                     removed (pass covariate=).
  3. Emits a per-feature verdict. A pure location/scale harmonizer always drives
     between-batch eta^2 toward zero, so residual eta^2 alone cannot separate a
     safe feature from a confounded one. The decisive question is whether the
     feature's association WITH THE OUTCOME survives harmonization. When an
     outcome column is supplied, the harness compares the pooled (raw) feature
     -> outcome AUROC against the within-batch (harmonized) AUROC:

       - ROBUST              : batch effect already negligible (eta^2 < 0.06);
                               pool the raw feature as-is.
       - NEEDS-HARMONIZATION : material batch effect raw (eta^2 >= 0.06), a
                               transform removes it (eta^2 -> ~0), AND the
                               outcome association is preserved after
                               harmonization (within-batch AUROC stays on the
                               same side of 0.5 and within ~0.05 of pooled) ->
                               safe to pool AFTER harmonization.
       - BATCH-CONFOUNDED    : material batch effect raw, and the pooled outcome
                               association COLLAPSES toward the null once the
                               batch axis is removed (|within-batch AUROC-0.5|
                               drops by >0.05, or flips sign). The apparent
                               signal was the cohort axis, not the biology
                               (Simpson's-paradox risk) -- exactly the failure
                               that sank naive pooled TMB. Do NOT pool; report
                               per clinical context only.

     With no outcome column the harness falls back to a variance-only verdict
     (ROBUST if eta^2_raw < bar, else NEEDS-HARMONIZATION if the transform clears
     the bar) and flags the feature as OUTCOME-UNTESTED in the note.

The eta^2 = 0.06 bar is Cohen's "medium effect" threshold; it is also where the
banked TMB harmonization crossed (raw 0.062 -> 0.006). Configurable via
`ETA2_BAR`; the AUROC-collapse tolerance is `AUROC_COLLAPSE`.

The harness is feature-agnostic: point it at any DataFrame + batch column +
list of numeric feature columns. `run_batch_robustness()` writes the report CSV
and the before/after variance figure.

Author: batch-robustness / eval session, evolutionary-RNA-state project.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ETA2_BAR = 0.06          # Cohen medium-effect; the TMB harmonization crossing point
RHO_PRESERVE = 0.80      # min within-batch Spearman rho(raw, harmonized) [diagnostic]
AUROC_COLLAPSE = 0.05    # max tolerated drop in |AUROC-0.5| from pooled to within-batch
RANDOM_STATE = 0


# --------------------------------------------------------------------------- #
# Batch-effect quantification
# --------------------------------------------------------------------------- #
def eta_squared(values: np.ndarray, batch: np.ndarray) -> float:
    """Fraction of total variance explained by batch (one-way ANOVA eta^2).

    eta^2 = SS_between / SS_total, in [0, 1]. 0 = batch explains nothing;
    1 = batch explains everything. This is the batch *effect size*.
    """
    v = np.asarray(values, float)
    b = np.asarray(batch)
    ok = ~np.isnan(v)
    v, b = v[ok], b[ok]
    if len(v) < 3 or len(np.unique(b)) < 2:
        return np.nan
    grand = v.mean()
    ss_total = np.sum((v - grand) ** 2)
    if ss_total == 0:
        return 0.0
    ss_between = 0.0
    for g in np.unique(b):
        vg = v[b == g]
        ss_between += len(vg) * (vg.mean() - grand) ** 2
    return float(ss_between / ss_total)


def kruskal_p(values: np.ndarray, batch: np.ndarray) -> float:
    """Kruskal-Wallis p-value that at least one batch has a different location."""
    v = np.asarray(values, float)
    b = np.asarray(batch)
    ok = ~np.isnan(v)
    v, b = v[ok], b[ok]
    groups = [v[b == g] for g in np.unique(b) if (b == g).sum() > 0]
    groups = [g for g in groups if len(g) > 0]
    if len(groups) < 2:
        return np.nan
    try:
        return float(stats.kruskal(*groups).pvalue)
    except ValueError:
        return np.nan


# --------------------------------------------------------------------------- #
# Harmonization transforms
# --------------------------------------------------------------------------- #
def harmonize_zscore(values: np.ndarray, batch: np.ndarray) -> np.ndarray:
    """Within-batch z-score: (x - mean_b) / sd_b for each batch b."""
    v = np.asarray(values, float).copy()
    b = np.asarray(batch)
    out = np.full_like(v, np.nan)
    for g in np.unique(b):
        m = b == g
        vg = v[m]
        ok = ~np.isnan(vg)
        if ok.sum() < 2:
            out[m] = 0.0
            continue
        mu, sd = np.nanmean(vg), np.nanstd(vg)
        out[m] = (vg - mu) / sd if sd > 0 else 0.0
    return out


def harmonize_rank(values: np.ndarray, batch: np.ndarray) -> np.ndarray:
    """Within-batch rank mapped to (0,1): scale-free, monotone-preserving."""
    v = np.asarray(values, float).copy()
    b = np.asarray(batch)
    out = np.full_like(v, np.nan)
    for g in np.unique(b):
        m = b == g
        vg = v[m]
        ok = ~np.isnan(vg)
        if ok.sum() == 0:
            continue
        r = stats.rankdata(vg[ok])
        vv = np.full(len(vg), np.nan)
        vv[ok] = (r - 0.5) / ok.sum()
        out[m] = vv
    return out


def combat(values: np.ndarray, batch: np.ndarray,
           covariate: np.ndarray | None = None) -> np.ndarray:
    """Parametric empirical-Bayes ComBat batch correction for a single feature.

    Reference implementation of Johnson, Li & Rabinovic (2007), Biostatistics
    8(1):118-127, specialized to one feature (1 x n). Removes additive (gamma)
    and multiplicative (delta) batch effects via an L/S model with empirical-
    Bayes shrinkage, optionally preserving a biological covariate design.

    Model:  y_ij = alpha + X beta + gamma_i + delta_i * eps_ij
    where i indexes batch, j sample; gamma/delta are shrunk toward the
    across-batch hyper-priors before the data are re-standardized.
    """
    y = np.asarray(values, float)
    b = np.asarray(batch)
    ok = ~np.isnan(y)
    if ok.sum() < 4 or len(np.unique(b[ok])) < 2:
        return y.copy()

    yv, bv = y[ok], b[ok]
    batches = np.unique(bv)
    n = len(yv)

    # design: intercept + optional covariate (one-hot) preserved through correction
    cols = [np.ones(n)]
    if covariate is not None:
        cov = np.asarray(covariate)[ok]
        for lv in np.unique(cov):
            cols.append((cov == lv).astype(float))
        # drop last covariate col to avoid collinearity with intercept
        cols = cols[:1] + cols[1:-1] if len(cols) > 2 else cols
    # batch design (one-hot, no intercept collinearity handled via grand mean)
    Bd = np.column_stack([(bv == g).astype(float) for g in batches])
    X = np.column_stack(cols + [Bd])

    # grand-mean standardization (pooled variance)
    beta, *_ = np.linalg.lstsq(X, yv, rcond=None)
    n_cov = len(cols)
    # batch sizes for weighted grand mean of batch terms
    ns = Bd.sum(axis=0)
    grand_batch = (beta[n_cov:] * ns).sum() / ns.sum()
    alpha = beta[:n_cov] @ np.array([1.0] + [0.0] * (n_cov - 1)) + grand_batch \
        if n_cov > 0 else grand_batch
    stand_mean = np.full(n, alpha)
    if n_cov > 1:
        stand_mean = X[:, :n_cov] @ beta[:n_cov] + grand_batch
    resid = yv - X @ beta
    var_pooled = (resid ** 2).mean()
    var_pooled = var_pooled if var_pooled > 0 else 1e-8
    Z = (yv - stand_mean) / np.sqrt(var_pooled)

    # per-batch L/S estimates
    gamma_hat, delta_hat = {}, {}
    for g in batches:
        zg = Z[bv == g]
        gamma_hat[g] = zg.mean()
        delta_hat[g] = zg.var(ddof=1) if len(zg) > 1 else 1.0

    # empirical-Bayes hyper-priors
    g_bar = np.mean(list(gamma_hat.values()))
    t2 = np.var(list(gamma_hat.values()), ddof=1) if len(batches) > 1 else 1.0
    d = np.array(list(delta_hat.values()))
    d_bar, d_var = d.mean(), d.var(ddof=1) if len(d) > 1 else (d.mean(), 1.0)
    a_prior = (2 * d_var + d_bar ** 2) / d_var if d_var > 0 else 1.0
    b_prior = (d_bar * d_var + d_bar ** 3) / d_var if d_var > 0 else 1.0

    # closed-form EB posterior (parametric shrinkage)
    gamma_star, delta_star = {}, {}
    for g in batches:
        ng = (bv == g).sum()
        gamma_star[g] = (ng * t2 * gamma_hat[g] + var_pooled * g_bar) / \
                        (ng * t2 + var_pooled)
        num = 0.5 * np.sum((Z[bv == g] - gamma_star[g]) ** 2) + b_prior
        den = ng / 2.0 + a_prior - 1.0
        delta_star[g] = num / den if den > 0 else delta_hat[g]

    # adjust
    Zc = Z.copy()
    for g in batches:
        m = bv == g
        Zc[m] = (Z[m] - gamma_star[g]) / np.sqrt(max(delta_star[g], 1e-8))
    corrected = Zc * np.sqrt(var_pooled) + stand_mean

    out = y.copy()
    out[ok] = corrected
    return out


HARMONIZERS = {"zscore": harmonize_zscore, "rank": harmonize_rank, "combat": combat}


# --------------------------------------------------------------------------- #
# Signal-preservation check
# --------------------------------------------------------------------------- #
def within_batch_preservation(raw, harmonized, batch) -> float:
    """Min across batches of Spearman rho(raw, harmonized) within that batch.

    A good harmonization removes the between-batch offset but preserves the
    WITHIN-batch ordering of samples (the biology). rho ~ 1 => order kept.
    """
    raw = np.asarray(raw, float)
    har = np.asarray(harmonized, float)
    b = np.asarray(batch)
    rhos = []
    for g in np.unique(b):
        m = b == g
        r, h = raw[m], har[m]
        ok = ~np.isnan(r) & ~np.isnan(h)
        if ok.sum() < 3 or np.nanstd(r[ok]) == 0 or np.nanstd(h[ok]) == 0:
            continue
        rho, _ = stats.spearmanr(r[ok], h[ok])
        if not np.isnan(rho):
            rhos.append(abs(rho))
    return float(min(rhos)) if rhos else np.nan


# --------------------------------------------------------------------------- #
# Outcome-association (the decisive test for confounding)
# --------------------------------------------------------------------------- #
from sklearn.metrics import roc_auc_score  # noqa: E402


def _pooled_auroc(values, y):
    v = np.asarray(values, float)
    y = np.asarray(y, float)
    ok = ~np.isnan(v) & ~np.isnan(y)
    v, y = v[ok], y[ok].astype(int)
    if len(np.unique(y)) < 2:
        return np.nan
    return float(roc_auc_score(y, v))


def _within_batch_auroc(values, y, batch):
    """Sample-weighted mean of within-batch AUROC (batch axis removed).

    Each batch with >=1 event and >=1 non-event contributes its own AUROC of the
    feature vs the outcome; averaging within batch is what removes the between-
    batch axis from the association (the pooled-vs-within contrast IS the
    Simpson's-paradox test).
    """
    v = np.asarray(values, float)
    y = np.asarray(y, float)
    b = np.asarray(batch)
    aucs, ws = [], []
    for g in np.unique(b):
        m = b == g
        vg, yg = v[m], y[m]
        ok = ~np.isnan(vg) & ~np.isnan(yg)
        vg, yg = vg[ok], yg[ok].astype(int)
        if len(np.unique(yg)) < 2 or len(yg) < 8:
            continue
        aucs.append(roc_auc_score(yg, vg))
        ws.append(len(yg))
    if not aucs:
        return np.nan, 0
    return float(np.average(aucs, weights=ws)), len(aucs)


# --------------------------------------------------------------------------- #
# Per-feature evaluation + verdict
# --------------------------------------------------------------------------- #
def evaluate_feature(df, feature, batch_col="cohort", method="zscore",
                     covariate_col=None, outcome_col=None) -> dict:
    keep = [feature, batch_col] + ([covariate_col] if covariate_col else []) \
        + ([outcome_col] if outcome_col else [])
    sub = df[keep].copy()
    sub = sub[sub[feature].notna()]
    base = dict(feature=feature, method=method, n=len(sub),
                n_batches=int(sub[batch_col].nunique()) if len(sub) else 0)
    if sub[batch_col].nunique() < 2 or len(sub) < 6:
        return dict(**base, eta2_raw=np.nan, kruskal_p_raw=np.nan,
                    eta2_harmonized=np.nan, within_batch_rho=np.nan,
                    auroc_pooled=np.nan, auroc_within_batch=np.nan,
                    verdict="INSUFFICIENT-DATA", note="need >=2 batches and >=6 samples")
    v = sub[feature].to_numpy(float)
    b = sub[batch_col].to_numpy()
    cov = sub[covariate_col].to_numpy() if covariate_col else None

    eta_raw = eta_squared(v, b)
    kp_raw = kruskal_p(v, b)
    fn = HARMONIZERS[method]
    vh = fn(v, b, cov) if method == "combat" else fn(v, b)
    eta_h = eta_squared(vh, b)
    rho = within_batch_preservation(v, vh, b)

    # outcome-association test
    auroc_pool = auroc_wb = np.nan
    n_batch_auc = 0
    note = ""
    if outcome_col and outcome_col in sub.columns:
        y = sub[outcome_col].to_numpy(float)
        auroc_pool = _pooled_auroc(v, y)
        auroc_wb, n_batch_auc = _within_batch_auroc(v, y, b)

    if eta_raw < ETA2_BAR:
        verdict = "ROBUST"
    elif outcome_col and not np.isnan(auroc_pool) and not np.isnan(auroc_wb):
        # decisive test: does the outcome association survive batch removal?
        drop = abs(auroc_pool - 0.5) - abs(auroc_wb - 0.5)
        flipped = np.sign(auroc_pool - 0.5) != np.sign(auroc_wb - 0.5)
        if drop > AUROC_COLLAPSE or flipped:
            verdict = "BATCH-CONFOUNDED"
            note = f"outcome assoc collapses pooled->within (dAUROC-from-0.5={drop:+.3f}"
            note += ", sign-flip" if flipped else ""
            note += f"); {n_batch_auc} batches testable"
        else:
            verdict = "NEEDS-HARMONIZATION"
            note = f"outcome assoc preserved after harmonization ({n_batch_auc} batches testable)"
    else:
        # no outcome available -> variance-only fallback
        verdict = "NEEDS-HARMONIZATION" if eta_h < ETA2_BAR else "BATCH-CONFOUNDED"
        note = "OUTCOME-UNTESTED: variance-only verdict (no outcome column)"

    return dict(**base, eta2_raw=round(eta_raw, 4), kruskal_p_raw=kp_raw,
                eta2_harmonized=round(eta_h, 4),
                within_batch_rho=round(rho, 3) if not np.isnan(rho) else np.nan,
                variance_reduction=round(eta_raw - eta_h, 4),
                auroc_pooled=round(auroc_pool, 3) if not np.isnan(auroc_pool) else np.nan,
                auroc_within_batch=round(auroc_wb, 3) if not np.isnan(auroc_wb) else np.nan,
                verdict=verdict, note=note)


def run_batch_robustness(df, features, batch_col="cohort", method="zscore",
                         covariate_col=None, outcome_col=None, outdir="results/eval",
                         report_name="batch_robustness_report.csv"):
    """Evaluate a list of features; write the report CSV. Returns the DataFrame."""
    rows = [evaluate_feature(df, f, batch_col, method, covariate_col, outcome_col)
            for f in features if f in df.columns]
    rep = pd.DataFrame(rows)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rep.to_csv(outdir / report_name, index=False)
    return rep


# --------------------------------------------------------------------------- #
# Figure
# --------------------------------------------------------------------------- #
def plot_batch_robustness(rep, path, title="Batch effect before vs after harmonization"):
    import matplotlib.pyplot as plt
    try:
        apply_figure_style(sizes=(8, 7, 6))  # noqa: F821
    except (NameError, Exception):
        import matplotlib as mpl
        mpl.rcParams.update({"figure.dpi": 120, "savefig.dpi": 300,
                             "savefig.bbox": "tight", "axes.spines.top": False,
                             "axes.spines.right": False})
    d = rep[rep["eta2_raw"].notna()].copy().sort_values("eta2_raw", ascending=True)
    y = np.arange(len(d))
    vcolor = {"ROBUST": "#1b7837", "NEEDS-HARMONIZATION": "#2166ac",
              "BATCH-CONFOUNDED": "#b2182b", "INSUFFICIENT-DATA": "0.6"}
    fig, ax = plt.subplots(figsize=(6.6, 0.44 * len(d) + 1.4))
    for yi, (_, r) in zip(y, d.iterrows()):
        c = vcolor.get(r["verdict"], "0.6")
        ax.plot([r["eta2_harmonized"], r["eta2_raw"]], [yi, yi],
                color="0.75", lw=1.2, zorder=1)
        ax.scatter(r["eta2_raw"], yi, marker="o", facecolors="none",
                   edgecolors=c, s=46, zorder=3)
        ax.scatter(r["eta2_harmonized"], yi, marker="o", color=c, s=46, zorder=3)
    ax.axvline(ETA2_BAR, ls="--", lw=0.9, color="0.4")
    ax.annotate(f"harmonization bar\n(eta$^2$={ETA2_BAR})", (ETA2_BAR, len(d) - 0.5),
                fontsize=6, color="0.4", ha="left", va="top",
                xytext=(4, 0), textcoords="offset points")
    ax.set_yticks(y)
    ax.set_yticklabels(d["feature"])
    ax.set_xlabel("Batch variance explained  (eta$^2$)   |   open = raw, filled = harmonized")
    ax.set_title(title)
    ax.margins(x=0.04)
    handles = [plt.Line2D([], [], marker="o", ls="", mfc=vcolor[k], mec=vcolor[k],
                          label=k) for k in ["ROBUST", "NEEDS-HARMONIZATION",
                                             "BATCH-CONFOUNDED"]]
    ax.legend(handles=handles, loc="lower right", fontsize=6)
    fig.savefig(path)
    plt.close(fig)
    return fig


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", default="results/analysis_frame.parquet")
    ap.add_argument("--outdir", default="results/eval")
    ap.add_argument("--method", default="zscore", choices=list(HARMONIZERS))
    args = ap.parse_args()
    frame = pd.read_parquet(args.parquet)
    feats = [c for c in frame.columns if frame[c].dtype.kind in "fi"]
    rep = run_batch_robustness(frame, feats, method=args.method, outdir=args.outdir)
    print(rep.to_string(index=False))
