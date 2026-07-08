"""baseline_benchmarks.py

Baseline benchmark harness for immune-checkpoint-blockade (ICB) response
prediction on the banked melanoma analysis frame (results/analysis_frame.parquet;
416 pretreatment samples, 5 cohorts).

Purpose
-------
Establish the *field-standard* baselines against which the latent evolutionary
RNA-state S will later be scored, and make that comparison a one-line addition:

    from baseline_benchmarks import run_benchmarks, MODELS, ModelSpec
    df_s = df.merge(S_scores, on="sampleId")          # S lands as a column
    MODELS.append(ModelSpec("Clinical + S", ["TMB_NONSYNONYMOUS","MUTATION_COUNT",
                                             "ICI_TARGET","S"], _feat_clinical_S))
    run_benchmarks(df_s, outdir="results/baselines")  # re-benchmarks with S

Baselines (aligned with the ICB modeling literature review):
  (a) TMB alone                       -- the clinically validated, weakly
                                         predictive single marker.
  (b) TIDE_RESPONDER alone            -- best-in-melanoma legacy transcriptomic
                                         classifier (Jiang 2018).
  (c) Clinical/burden panel           -- TMB + MUTATION_COUNT + ICI_TARGET, the
                                         "compact tabular" idea (Chowell 2021 /
                                         LORIS-style deployable model).
  (d) S / RNA-state                   -- PENDING placeholder slot.

Design guarantees
-----------------
* CV-leakage guard: StandardScaler + feature handling sit INSIDE a Pipeline that
  is fit only on each training fold; the outer train/test split is above all
  scaling and hyper-parameter selection (nested CV).
* AUROC with stratified bootstrap 95% CIs; PR-AUC (average precision) reported
  alongside. Pooled and per-cohort.
* Cohort-held-out external validation: train {gide2019, liu2019} -> test
  {riaz2017, hugo2016}. Missingness is handled explicitly and reported (dfci2019
  has no RESPONDER/TIDE; gide2019 has no TMB, so TMB-based models cannot use it).

Author: baseline-harness session, evolutionary-RNA-state project.
"""
from __future__ import annotations

import json
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Sequence

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve

RANDOM_STATE = 0
C_GRID = [0.01, 0.1, 1.0, 10.0]
N_BOOT = 2000

EXTERNAL_TRAIN = ["gide2019", "liu2019"]
EXTERNAL_TEST = ["riaz2017", "hugo2016"]

# --------------------------------------------------------------------------- #
# Feature builders. Each returns an (n, k) float array from the frame.
# Burden features are log1p-transformed (heavy right skew); TIDE is already
# binary; ICI_TARGET is encoded as a single combo indicator (CTLA4+PD1 vs PD1).
# --------------------------------------------------------------------------- #
def _log_tmb(df):
    return np.log1p(df["TMB_NONSYNONYMOUS"].to_numpy(float))


def _log_mut(df):
    return np.log1p(df["MUTATION_COUNT"].to_numpy(float))


def _ici_combo(df):
    return (df["ICI_TARGET"].astype("string") == "CTLA4 PD1").to_numpy(float)


def _feat_tmb(df):
    return _log_tmb(df).reshape(-1, 1)


def _feat_tide(df):
    return df["TIDE_RESPONDER"].to_numpy(float).reshape(-1, 1)


def _feat_clinical(df):
    return np.column_stack([_log_tmb(df), _log_mut(df), _ici_combo(df)])


@dataclass
class ModelSpec:
    name: str
    required: Sequence[str]          # columns that must be non-null to include a row
    feature_fn: Callable            # df -> (n, k) array
    pending: bool = False           # True = reserved slot, no data yet


# Registry. Append here (or pass extra_models=) to add S later.
MODELS: list[ModelSpec] = [
    ModelSpec("TMB", ["TMB_NONSYNONYMOUS"], _feat_tmb),
    ModelSpec("TIDE", ["TIDE_RESPONDER"], _feat_tide),
    ModelSpec("Clinical panel (TMB+MUT+ICI)",
              ["TMB_NONSYNONYMOUS", "MUTATION_COUNT", "ICI_TARGET"], _feat_clinical),
    ModelSpec("S / RNA-state (pending)", ["S"], lambda df: np.empty((len(df), 1)),
              pending=True),
]

# --------------------------------------------------------------------------- #
# Core estimator + CV helpers
# --------------------------------------------------------------------------- #
def _pipeline():
    return Pipeline([
        ("scale", StandardScaler()),
        ("lr", LogisticRegression(max_iter=2000, solver="liblinear")),
    ])


def _safe_k(y, kmax=5):
    """Largest feasible number of stratified folds given the minority class."""
    minc = int(np.min(np.bincount(y.astype(int))))
    return max(2, min(kmax, minc))


def nested_cv_oof(X, y, kmax=5, seed=RANDOM_STATE):
    """Out-of-fold predicted probabilities from NESTED stratified CV.

    Outer StratifiedKFold defines train/test (above all scaling & tuning); an
    inner StratifiedKFold GridSearch selects the L2 penalty C on each outer-train
    fold only. Guarantees no scaling/selection leakage across the split.
    """
    y = np.asarray(y).astype(int)
    k = _safe_k(y, kmax)
    outer = StratifiedKFold(n_splits=k, shuffle=True, random_state=seed)
    oof = np.full(len(y), np.nan)
    for tr, te in outer.split(X, y):
        ki = _safe_k(y[tr], kmax)
        inner = StratifiedKFold(n_splits=ki, shuffle=True, random_state=seed)
        gs = GridSearchCV(_pipeline(), {"lr__C": C_GRID}, cv=inner,
                          scoring="roc_auc", n_jobs=1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gs.fit(X[tr], y[tr])
            oof[te] = gs.predict_proba(X[te])[:, 1]
    return oof


def fit_select(X, y, kmax=5, seed=RANDOM_STATE):
    """Fit pipeline with inner-CV-selected C on ALL of (X,y). For external holdout."""
    y = np.asarray(y).astype(int)
    ki = _safe_k(y, kmax)
    inner = StratifiedKFold(n_splits=ki, shuffle=True, random_state=seed)
    gs = GridSearchCV(_pipeline(), {"lr__C": C_GRID}, cv=inner,
                      scoring="roc_auc", n_jobs=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gs.fit(X, y)
    return gs.best_estimator_


# --------------------------------------------------------------------------- #
# Metrics + stratified bootstrap CIs
# --------------------------------------------------------------------------- #
def _boot_ci(y, s, metric, B=N_BOOT, seed=RANDOM_STATE):
    y = np.asarray(y).astype(int)
    s = np.asarray(s, float)
    rng = np.random.default_rng(seed)
    i0, i1 = np.where(y == 0)[0], np.where(y == 1)[0]
    if len(i0) == 0 or len(i1) == 0:
        return (np.nan, np.nan)
    vals = np.empty(B)
    for b in range(B):
        bi = np.concatenate([rng.choice(i0, len(i0), True),
                             rng.choice(i1, len(i1), True)])
        vals[b] = metric(y[bi], s[bi])
    return tuple(np.nanpercentile(vals, [2.5, 97.5]))


def score_predictions(y, s, seed=RANDOM_STATE):
    y = np.asarray(y).astype(int)
    s = np.asarray(s, float)
    ok = ~np.isnan(s)
    y, s = y[ok], s[ok]
    if len(np.unique(y)) < 2:
        return dict(auroc=np.nan, auroc_lo=np.nan, auroc_hi=np.nan,
                    prauc=np.nan, prauc_lo=np.nan, prauc_hi=np.nan, n=len(y),
                    n_pos=int(y.sum()))
    auroc = roc_auc_score(y, s)
    prauc = average_precision_score(y, s)
    alo, ahi = _boot_ci(y, s, roc_auc_score, seed=seed)
    plo, phi = _boot_ci(y, s, average_precision_score, seed=seed)
    return dict(auroc=auroc, auroc_lo=alo, auroc_hi=ahi,
                prauc=prauc, prauc_lo=plo, prauc_hi=phi,
                n=len(y), n_pos=int(y.sum()))


# --------------------------------------------------------------------------- #
# Evaluation drivers
# --------------------------------------------------------------------------- #
def _eligible(df, spec, target="RESPONDER"):
    m = df[target].notna()
    for c in spec.required:
        m &= df[c].notna() if c in df.columns else False
    return df[m].copy()


def eval_model(df, spec, target="RESPONDER"):
    """Return a list of metric rows for one model: pooled_cv, each cohort, external."""
    rows = []
    if spec.pending:
        rows.append(dict(model=spec.name, evaluation="pooled_cv", **_nan_metrics(),
                         note="PENDING: awaiting latent RNA-state S"))
        return rows, None

    sub = _eligible(df, spec, target)
    y = sub[target].astype(int).to_numpy()
    X = spec.feature_fn(sub)

    # (1) pooled CV
    oof = nested_cv_oof(X, y)
    m = score_predictions(y, oof)
    rows.append(dict(model=spec.name, evaluation="pooled_cv", **m,
                     cohorts_used="|".join(sorted(sub["cohort"].unique()))))
    roc_data = _roc_xy(y, oof)

    # (2) per-cohort CV
    for coh, g in sub.groupby("cohort"):
        yc = g[target].astype(int).to_numpy()
        Xc = spec.feature_fn(g)
        if len(np.unique(yc)) < 2 or int(np.min(np.bincount(yc))) < 3:
            rows.append(dict(model=spec.name, evaluation=f"cohort:{coh}",
                             **_nan_metrics(n=len(yc), n_pos=int(yc.sum())),
                             note="too few events for within-cohort CV"))
            continue
        oofc = nested_cv_oof(Xc, yc, kmax=5)
        rows.append(dict(model=spec.name, evaluation=f"cohort:{coh}",
                         **score_predictions(yc, oofc)))

    # (3) external cohort holdout
    rows.append(_external_eval(df, spec, target))
    return rows, roc_data


def _external_eval(df, spec, target="RESPONDER"):
    tr = _eligible(df[df["cohort"].isin(EXTERNAL_TRAIN)], spec, target)
    te = _eligible(df[df["cohort"].isin(EXTERNAL_TEST)], spec, target)
    note_bits = [f"train={'+'.join(sorted(tr['cohort'].unique()))}(n={len(tr)})",
                 f"test={'+'.join(sorted(te['cohort'].unique()))}(n={len(te)})"]
    dropped = sorted(set(EXTERNAL_TRAIN) - set(tr["cohort"].unique()))
    if dropped:
        note_bits.append(f"DROPPED_from_train(missing feature): {','.join(dropped)}")
    ytr = tr[target].astype(int).to_numpy()
    if len(tr) < 10 or len(np.unique(ytr)) < 2 or len(te) < 5:
        return dict(model=spec.name, evaluation="external_holdout",
                    **_nan_metrics(n=len(te), n_pos=int(te[target].sum()) if len(te) else 0),
                    note="; ".join(note_bits + ["insufficient for external fit"]))
    est = fit_select(spec.feature_fn(tr), ytr)
    yte = te[target].astype(int).to_numpy()
    ste = est.predict_proba(spec.feature_fn(te))[:, 1]
    m = score_predictions(yte, ste)
    return dict(model=spec.name, evaluation="external_holdout", **m,
                note="; ".join(note_bits))


def _nan_metrics(n=0, n_pos=0):
    return dict(auroc=np.nan, auroc_lo=np.nan, auroc_hi=np.nan,
                prauc=np.nan, prauc_lo=np.nan, prauc_hi=np.nan, n=n, n_pos=n_pos)


def _roc_xy(y, s):
    ok = ~np.isnan(s)
    y, s = np.asarray(y)[ok], np.asarray(s)[ok]
    if len(np.unique(y)) < 2:
        return None
    fpr, tpr, _ = roc_curve(y, s)
    return dict(fpr=fpr.tolist(), tpr=tpr.tolist(), auroc=float(roc_auc_score(y, s)),
                n=int(len(y)))


# --------------------------------------------------------------------------- #
# Plotting (self-contained publication style; no external skill dependency)
# --------------------------------------------------------------------------- #
def _setup_style():
    """Apply the figure-style skill's apply_figure_style() when available
    (preferred: brings the skill's correctness defaults). Fall back to an
    equivalent self-contained rcParams block so the module stays importable
    and runnable without the skill kernel plugin."""
    try:
        apply_figure_style(sizes=(8, 7, 6))  # noqa: F821 (injected by skill)
        return
    except (NameError, Exception):
        pass
    import matplotlib as mpl
    mpl.rcParams.update({
        "figure.dpi": 120, "savefig.dpi": 300, "savefig.bbox": "tight",
        "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
        "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 7,
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.direction": "out", "ytick.direction": "out",
        "legend.frameon": False, "axes.linewidth": 0.8,
        "font.family": "sans-serif",
    })


def plot_roc(roc_by_model, path):
    import matplotlib.pyplot as plt
    _setup_style()
    fig, ax = plt.subplots(figsize=(4.2, 4.0))
    colors = ["#2166ac", "#b2182b", "#1b7837", "#762a83"]
    for (name, rd), c in zip(roc_by_model.items(), colors):
        if rd is None:
            continue
        ax.plot(rd["fpr"], rd["tpr"], color=c, lw=1.8,
                label=f"{name} (AUROC {rd['auroc']:.2f}, n={rd.get('n','?')})")
    ax.plot([0, 1], [0, 1], ls="--", lw=0.9, color="0.6")
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("Pooled cross-validated ROC, ICB response baselines")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal")
    ax.legend(loc="lower right")
    fig.savefig(path)
    plt.close(fig)


def plot_forest(metrics_df, path):
    import matplotlib.pyplot as plt
    _setup_style()
    order_eval = ["pooled_cv", "external_holdout"]
    d = metrics_df[metrics_df["evaluation"].isin(order_eval)].copy()
    # keep model order as in registry, evaluation grouped
    model_order = [m.name for m in MODELS]
    d["mrank"] = d["model"].map({m: i for i, m in enumerate(model_order)})
    d["erank"] = d["evaluation"].map({e: i for i, e in enumerate(order_eval)})
    d = d.sort_values(["mrank", "erank"]).reset_index(drop=True)

    labels, y_auroc, y_lo, y_hi, is_pending = [], [], [], [], []
    for _, r in d.iterrows():
        labels.append(f"{r['model']}  \u2013  {r['evaluation']}")
        y_auroc.append(r["auroc"])
        y_lo.append(r["auroc_lo"])
        y_hi.append(r["auroc_hi"])
        is_pending.append(bool(pd.isna(r["auroc"])))

    ypos = np.arange(len(labels))[::-1]
    fig, ax = plt.subplots(figsize=(6.4, 0.42 * len(labels) + 1.2))
    for yp, a, lo, hi, pend in zip(ypos, y_auroc, y_lo, y_hi, is_pending):
        if pend:
            ax.scatter([0.5], [yp], marker="s", facecolors="none",
                       edgecolors="0.5", s=42, zorder=3)
            ax.annotate("slot reserved \u2013 fill when S lands", (0.505, yp),
                        va="center", ha="left", fontsize=6, color="0.45")
        else:
            ax.plot([lo, hi], [yp, yp], color="0.35", lw=1.4, zorder=2)
            ax.scatter([a], [yp], color="#2166ac", s=30, zorder=3)
    ax.axvline(0.5, ls="--", lw=0.9, color="0.6")
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels)
    ax.set_xlim(0.3, 1.0)
    ax.set_xlabel("AUROC (point + bootstrap 95% CI); dashed line = chance (0.5)")
    ax.set_title("Baseline AUROC by model and evaluation\n(RNA-state S slot pending)")
    fig.savefig(path)
    plt.close(fig)


def plot_km_by_tmb(df, path):
    import matplotlib.pyplot as plt
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import multivariate_logrank_test
    _setup_style()
    s = df.dropna(subset=["TMB_NONSYNONYMOUS", "OS_MONTHS", "OS_STATUS"]).copy()
    s["event"] = s["OS_STATUS"].astype(str).str.startswith("1").astype(int)
    s["tert"] = pd.qcut(s["TMB_NONSYNONYMOUS"], 3, labels=["Low", "Mid", "High"])
    fig, ax = plt.subplots(figsize=(4.6, 4.0))
    colors = {"Low": "#2166ac", "Mid": "#999999", "High": "#b2182b"}
    kmf = KaplanMeierFitter()
    for lvl in ["Low", "Mid", "High"]:
        g = s[s["tert"] == lvl]
        kmf.fit(g["OS_MONTHS"], g["event"], label=f"{lvl} TMB (n={len(g)})")
        kmf.plot_survival_function(ax=ax, color=colors[lvl], lw=1.8, ci_show=False)
    lr = multivariate_logrank_test(s["OS_MONTHS"], s["tert"], s["event"])
    p = lr.p_value
    ax.set_xlabel("Overall survival (months)")
    ax.set_ylabel("Survival probability")
    ax.set_ylim(0, 1.02)
    ptxt = "p < 0.001" if p < 1e-3 else f"log-rank p = {p:.3f}"
    ax.set_title("Overall survival by TMB tertile")
    ax.annotate(ptxt, (0.96, 0.96), xycoords="axes fraction", ha="right",
                va="top", fontsize=7)
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.90))
    fig.savefig(path)
    plt.close(fig)
    return float(p)


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def run_benchmarks(df, outdir="results/baselines", extra_models=None, target="RESPONDER"):
    """Run all baselines, write metrics CSV + figures. Returns the metrics DataFrame.

    To benchmark S: merge its column into df and append a ModelSpec (see header),
    then call this again.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    models = list(MODELS) + list(extra_models or [])

    all_rows, roc_by_model = [], {}
    for spec in models:
        rows, roc = eval_model(df, spec, target=target)
        all_rows.extend(rows)
        if roc is not None:
            roc_by_model[spec.name] = roc

    cols = ["model", "evaluation", "auroc", "auroc_lo", "auroc_hi",
            "prauc", "prauc_lo", "prauc_hi", "n", "n_pos", "cohorts_used", "note"]
    mdf = pd.DataFrame(all_rows)
    for c in cols:
        if c not in mdf.columns:
            mdf[c] = np.nan
    mdf = mdf[cols]
    mdf.to_csv(outdir / "baseline_metrics.csv", index=False)

    plot_roc(roc_by_model, outdir / "fig_baseline_roc.png")
    plot_forest(mdf, outdir / "fig_baseline_forest.png")
    km_p = plot_km_by_tmb(df, outdir / "km_by_tmb_tertile.png")

    (outdir / "km_logrank_p.json").write_text(json.dumps({"logrank_p": km_p}))
    return mdf


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", default="results/analysis_frame.parquet")
    ap.add_argument("--outdir", default="results/baselines")
    args = ap.parse_args()
    frame = pd.read_parquet(args.parquet)
    metrics = run_benchmarks(frame, outdir=args.outdir)
    with pd.option_context("display.width", 160, "display.max_columns", 20):
        print(metrics.to_string(index=False))
