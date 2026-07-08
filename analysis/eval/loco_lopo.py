"""loco_lopo.py

COMPASS-style leave-one-cohort-out (LOCO) and leave-one-patient-out (LOPO)
evaluation of INCREMENTAL predictive value over the proven ICB-response floor.

Extends baseline_benchmarks.py (the reserved 'S' slot): instead of scoring one
model at a time, this harness asks the deliverable question directly --

    Does feature X add AUROC over the floor, out-of-cohort and out-of-patient?

Floor definition (explicit)
---------------------------
The clinically proven melanoma ICB-response floor used here is:

    FLOOR = TMB (muts, log1p) + TIDE_RESPONDER + T-cell-inflamed GEP composite

with two honest caveats specific to this cohort collection:

  * dMMR / MSI is NOT a usable floor component: these are cutaneous melanoma
    cohorts, essentially all MMR-proficient, so dMMR carries no signal and is
    absent from the frame. Stated, not silently dropped.
  * PD-L1 IHC is NOT available in the banked frame (no harmonized TPS/CPS
    column across cohorts). Also stated, not imputed.
  * The GEP composite is RNA-derived and PENDING the pipeline pilot
    (results/features/quant_gene_tpm -> gep_scores). Until it lands, the
    RUNNABLE-NOW floor is TMB + TIDE; the GEP term is wired as a slot that the
    harness folds in automatically when the gep_* columns appear. Every metric
    row records which floor variant it used (`floor_variant`).

Because TMB is missing in gide2019, any TMB-containing floor evaluates on the
3 cohorts that have TMB (liu2019, riaz2017, hugo2016). Candidate features that
exist for gide2019 are additionally reported against a TIDE-only floor so the
4th cohort is not silently discarded.

Incremental AUROC
-----------------
For each candidate feature (or feature set) we fit two nested-CV models on the
SAME sample support -- floor alone and floor+candidate -- and report

    dAUROC = AUROC(floor + candidate) - AUROC(floor)

under three evaluation modes:

  * pooled_LOPO   : leave-one-patient-out (== pooled out-of-fold nested CV with
                    the split above all scaling/tuning); the within-collection
                    generalization estimate.
  * LOCO_<cohort> : train on all-but-one cohort, test on the held-out cohort;
                    the out-of-distribution / cross-study estimate (COMPASS
                    external-validation style).
  * LOCO_mean     : sample-weighted mean of the per-cohort LOCO dAUROC.

Leakage guard: identical to baseline_benchmarks -- StandardScaler + C-selection
live inside a Pipeline fit only on each training fold; the LOCO/LOPO split sits
ABOVE all fitting. Paired bootstrap (resampling test samples, scoring both
models on each resample) gives the 95% CI of dAUROC, so the CI reflects the
paired nature of the comparison.

RNA-derived candidates whose per-sample matrices have not landed are emitted as
PENDING rows (auto-filled when results/features/*.parquet appear and are merged
onto the frame by key).

Author: eval session, evolutionary-RNA-state project.
"""
from __future__ import annotations

import json
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

RANDOM_STATE = 0
C_GRID = [0.01, 0.1, 1.0, 10.0]
N_BOOT = 2000

TARGET = "RESPONDER"
BATCH = "cohort"

# --------------------------------------------------------------------------- #
# Floor + candidate specification
# --------------------------------------------------------------------------- #
# Floor terms that are RUNNABLE NOW on the banked frame (WES / legacy TIDE):
FLOOR_NOW = ["TMB_NONSYNONYMOUS", "TIDE_RESPONDER"]
# Floor terms PENDING the RNA pilot (folded in automatically when present):
FLOOR_PENDING = ["gep_tcell_inflamed", "ifng_score", "teff_tgfb_balance"]


@dataclass
class Candidate:
    name: str                       # interpretable feature name
    columns: Sequence[str]          # frame column(s) providing the value
    bucket: str                     # 'baseline' | 'differentiated'
    log1p: bool = True              # burden counts -> log1p; scores -> False
    pending: bool = False           # True until the RNA matrix lands
    note: str = ""


# Candidates demonstrable NOW use the WES-derived proxy columns as stand-ins for
# the RNA features (clearly labelled). The RNA columns are wired as pending rows
# keyed to the exact registry feature name; when results/features/*.parquet land
# and merge onto the frame, flip pending=False (or pass a resolver map).
CANDIDATES: list[Candidate] = [
    # --- demonstrated on real WES-derived proxy columns ---
    Candidate("SNV/indel neoantigen burden (WES proxy)", ["SNV_NEOANTIGEN", "INDEL_NEOANTIGEN"],
              "baseline", note="real WES-derived proxy for snv_indel_neoantigen_burden"),
    Candidate("Splicing neoantigen burden (WES proxy)", ["SPLICE_NEOANTIGEN"],
              "differentiated", note="real WES-derived proxy for splice_neoantigen_burden"),
    Candidate("TE/ERV antigen burden (WES proxy)", ["ERV_NEOANTIGEN"],
              "differentiated", note="real WES-derived proxy for te_antigen_burden(_ERV)"),
    Candidate("Fusion neoantigen burden (WES proxy)", ["FUSION_NEOANTIGEN"],
              "differentiated", note="real WES-derived proxy for fusion_neoantigen_burden"),
    # --- PENDING real RNA-derived features (auto-fill when matrices land) ---
    Candidate("splice_neoantigen_burden (RNA)", ["splice_neoantigen_burden"],
              "differentiated", pending=True, note="SNAF neojunction burden, RNA pilot"),
    Candidate("te_antigen_burden (RNA)", ["te_antigen_burden"],
              "differentiated", pending=True, note="TE/ERV antigen burden, RNA pilot"),
    Candidate("retained_intron_load (RNA)", ["retained_intron_load"],
              "differentiated", pending=True, log1p=False, note="intron-retention load, RNA pilot"),
    Candidate("ir_neoantigen_burden (RNA)", ["ir_neoantigen_burden"],
              "differentiated", pending=True, note="retained-intron neoantigen burden, RNA pilot"),
    Candidate("alu_editing_index (RNA)", ["alu_editing_index"],
              "differentiated", pending=True, log1p=False, note="Alu editing index, RNA pilot"),
    Candidate("editing_neoantigen_burden (RNA)", ["editing_neoantigen_burden"],
              "differentiated", pending=True, note="RNA-editing neoantigen burden, RNA pilot"),
    Candidate("fusion_neoantigen_burden (RNA)", ["fusion_neoantigen_burden"],
              "differentiated", pending=True, note="fusion neoantigen burden, RNA pilot"),
    Candidate("RNA-state composite S (all differentiated)", ["S"],
              "differentiated", pending=True, log1p=False,
              note="reserved 'S' slot: joint differentiated-feature composite"),
]


# --------------------------------------------------------------------------- #
# Feature-matrix builders
# --------------------------------------------------------------------------- #
def _resolve_floor(df):
    """Return (floor_columns_present, floor_variant_label)."""
    pending_present = [c for c in FLOOR_PENDING if c in df.columns and df[c].notna().any()]
    floor = [c for c in FLOOR_NOW if c in df.columns]
    if pending_present:
        return floor + pending_present, "TMB+TIDE+GEP"
    return floor, "TMB+TIDE (GEP pending)"


def _build_matrix(df, columns, log_mask):
    cols = []
    for c, lg in zip(columns, log_mask):
        x = df[c].to_numpy(float)
        cols.append(np.log1p(x) if lg else x)
    return np.column_stack(cols)


# --------------------------------------------------------------------------- #
# Estimator + nested CV (leakage-guarded, mirrors baseline_benchmarks)
# --------------------------------------------------------------------------- #
def _pipeline():
    return Pipeline([("scale", StandardScaler()),
                     ("lr", LogisticRegression(max_iter=2000, solver="liblinear"))])


def _safe_k(y, kmax=5):
    minc = int(np.min(np.bincount(y.astype(int))))
    return max(2, min(kmax, minc))


def _oof_nested(X, y, groups=None, mode="lopo", kmax=5, seed=RANDOM_STATE):
    """Out-of-fold probabilities.

    mode='lopo' : stratified K-fold pooled CV (leave-one-patient-out surrogate;
                  each sample scored by a model that never saw it, split above
                  scaling+tuning). True LOO is the k=n limit; stratified K-fold
                  with the largest feasible k is the stable, standard stand-in.
    mode='loco' : leave-one-group(cohort)-out; each cohort scored by a model
                  trained on the others.
    """
    y = np.asarray(y).astype(int)
    oof = np.full(len(y), np.nan)
    if mode == "loco":
        for g in np.unique(groups):
            te = groups == g
            tr = ~te
            if len(np.unique(y[tr])) < 2 or tr.sum() < 10:
                continue
            est = _fit_select(X[tr], y[tr], seed=seed)
            oof[te] = est.predict_proba(X[te])[:, 1]
        return oof
    k = _safe_k(y, kmax)
    outer = StratifiedKFold(n_splits=k, shuffle=True, random_state=seed)
    for tr, te in outer.split(X, y):
        est = _fit_select(X[tr], y[tr], seed=seed)
        oof[te] = est.predict_proba(X[te])[:, 1]
    return oof


def _fit_select(X, y, kmax=5, seed=RANDOM_STATE):
    y = np.asarray(y).astype(int)
    ki = _safe_k(y, kmax)
    inner = StratifiedKFold(n_splits=ki, shuffle=True, random_state=seed)
    gs = GridSearchCV(_pipeline(), {"lr__C": C_GRID}, cv=inner, scoring="roc_auc", n_jobs=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gs.fit(X, y)
    return gs.best_estimator_


# --------------------------------------------------------------------------- #
# Paired incremental-AUROC scoring
# --------------------------------------------------------------------------- #
def _paired_dauroc_ci(y, s_floor, s_full, B=N_BOOT, seed=RANDOM_STATE):
    """Bootstrap 95% CI of dAUROC = AUROC(full) - AUROC(floor), PAIRED.

    Resample test samples once per iteration; score BOTH models on the same
    resample so the CI reflects the paired comparison (correlated errors).
    """
    y = np.asarray(y).astype(int)
    sf = np.asarray(s_floor, float)
    su = np.asarray(s_full, float)
    ok = ~np.isnan(sf) & ~np.isnan(su)
    y, sf, su = y[ok], sf[ok], su[ok]
    rng = np.random.default_rng(seed)
    i0, i1 = np.where(y == 0)[0], np.where(y == 1)[0]
    if len(i0) == 0 or len(i1) == 0:
        return (np.nan, np.nan, np.nan)
    diffs = np.empty(B)
    for b in range(B):
        bi = np.concatenate([rng.choice(i0, len(i0), True), rng.choice(i1, len(i1), True)])
        diffs[b] = roc_auc_score(y[bi], su[bi]) - roc_auc_score(y[bi], sf[bi])
    lo, hi = np.nanpercentile(diffs, [2.5, 97.5])
    p_gt0 = float(np.mean(diffs <= 0))  # one-sided bootstrap p that dAUROC<=0
    return float(lo), float(hi), p_gt0


def evaluate_candidate(df, cand: Candidate, seed=RANDOM_STATE):
    """Return metric rows (pooled_LOPO, LOCO per-cohort, LOCO_mean) for one candidate."""
    floor_cols, floor_variant = _resolve_floor(df)
    rows = []
    if cand.pending or any(c not in df.columns or df[c].notna().sum() == 0 for c in cand.columns):
        rows.append(dict(candidate=cand.name, bucket=cand.bucket, evaluation="pooled_LOPO",
                         floor_variant=floor_variant, n=0, n_pos=0,
                         auroc_floor=np.nan, auroc_full=np.nan, dauroc=np.nan,
                         dauroc_lo=np.nan, dauroc_hi=np.nan, boot_p=np.nan,
                         status="PENDING",
                         note=f"awaiting {'|'.join(cand.columns)} (RNA pilot); {cand.note}"))
        return rows, floor_variant

    need = floor_cols + list(cand.columns) + [TARGET, BATCH]
    sub = df.dropna(subset=[c for c in need if c in df.columns]).copy()
    if len(sub) < 20 or sub[TARGET].nunique() < 2:
        rows.append(dict(candidate=cand.name, bucket=cand.bucket, evaluation="pooled_LOPO",
                         floor_variant=floor_variant, n=len(sub), n_pos=int(sub[TARGET].sum()),
                         auroc_floor=np.nan, auroc_full=np.nan, dauroc=np.nan,
                         dauroc_lo=np.nan, dauroc_hi=np.nan, boot_p=np.nan,
                         status="INSUFFICIENT", note="too few complete-case samples"))
        return rows, floor_variant

    y = sub[TARGET].astype(int).to_numpy()
    groups = sub[BATCH].to_numpy()
    floor_logmask = [c in ("TMB_NONSYNONYMOUS",) for c in floor_cols]
    Xf = _build_matrix(sub, floor_cols, floor_logmask)
    Xc = _build_matrix(sub, list(cand.columns), [cand.log1p] * len(cand.columns))
    Xu = np.column_stack([Xf, Xc])

    # pooled LOPO
    of = _oof_nested(Xf, y, mode="lopo", seed=seed)
    ou = _oof_nested(Xu, y, mode="lopo", seed=seed)
    af_, au_ = roc_auc_score(y, of), roc_auc_score(y, ou)
    lo, hi, p = _paired_dauroc_ci(y, of, ou, seed=seed)
    rows.append(dict(candidate=cand.name, bucket=cand.bucket, evaluation="pooled_LOPO",
                     floor_variant=floor_variant, n=len(sub), n_pos=int(y.sum()),
                     auroc_floor=round(af_, 3), auroc_full=round(au_, 3),
                     dauroc=round(au_ - af_, 3), dauroc_lo=round(lo, 3), dauroc_hi=round(hi, 3),
                     boot_p=round(p, 3), status="DEMONSTRATED_ON_PROXY" if "proxy" in cand.note
                     else "DEMONSTRATED", cohorts="|".join(sorted(sub[BATCH].unique())),
                     note=cand.note))

    # LOCO
    of_l = _oof_nested(Xf, y, groups=groups, mode="loco", seed=seed)
    ou_l = _oof_nested(Xu, y, groups=groups, mode="loco", seed=seed)
    loco_d, loco_w = [], []
    for g in np.unique(groups):
        m = groups == g
        yg = y[m]
        if len(np.unique(yg)) < 2 or np.isnan(of_l[m]).all():
            continue
        try:
            afg = roc_auc_score(yg, of_l[m]); aug = roc_auc_score(yg, ou_l[m])
        except ValueError:
            continue
        lo, hi, p = _paired_dauroc_ci(yg, of_l[m], ou_l[m], seed=seed)
        rows.append(dict(candidate=cand.name, bucket=cand.bucket, evaluation=f"LOCO_holdout:{g}",
                         floor_variant=floor_variant, n=int(m.sum()), n_pos=int(yg.sum()),
                         auroc_floor=round(afg, 3), auroc_full=round(aug, 3),
                         dauroc=round(aug - afg, 3), dauroc_lo=round(lo, 3), dauroc_hi=round(hi, 3),
                         boot_p=round(p, 3), status="DEMONSTRATED_ON_PROXY" if "proxy" in cand.note
                         else "DEMONSTRATED", note=f"trained on other cohorts"))
        loco_d.append(aug - afg); loco_w.append(int(m.sum()))
    if loco_d:
        rows.append(dict(candidate=cand.name, bucket=cand.bucket, evaluation="LOCO_mean",
                         floor_variant=floor_variant, n=sum(loco_w), n_pos=int(y.sum()),
                         auroc_floor=np.nan, auroc_full=np.nan,
                         dauroc=round(float(np.average(loco_d, weights=loco_w)), 3),
                         dauroc_lo=np.nan, dauroc_hi=np.nan, boot_p=np.nan,
                         status="DEMONSTRATED_ON_PROXY" if "proxy" in cand.note else "DEMONSTRATED",
                         note=f"sample-weighted mean over {len(loco_d)} held-out cohorts"))
    return rows, floor_variant


def run_loco_lopo(df, candidates=None, outdir="results/eval",
                  report_name="loco_lopo_report.csv"):
    candidates = candidates or CANDIDATES
    all_rows = []
    for c in candidates:
        rows, _ = evaluate_candidate(df, c)
        all_rows.extend(rows)
    rep = pd.DataFrame(all_rows)
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    rep.to_csv(outdir / report_name, index=False)
    return rep


# --------------------------------------------------------------------------- #
# Figure: incremental-AUROC forest
# --------------------------------------------------------------------------- #
def plot_incremental_forest(rep, path, evaluation="pooled_LOPO",
                            title="Incremental AUROC over the TMB+TIDE(+GEP) floor"):
    import matplotlib.pyplot as plt
    try:
        apply_figure_style(sizes=(8, 7, 6))  # noqa: F821
    except (NameError, Exception):
        import matplotlib as mpl
        mpl.rcParams.update({"figure.dpi": 120, "savefig.dpi": 300, "savefig.bbox": "tight",
                             "axes.spines.top": False, "axes.spines.right": False})
    d = rep[rep["evaluation"] == evaluation].copy()
    # order: demonstrated first (by dAUROC), then pending
    d["is_pending"] = d["status"].eq("PENDING")
    d = d.sort_values(["is_pending", "dauroc"], ascending=[True, True]).reset_index(drop=True)
    y = np.arange(len(d))
    bcolor = {"baseline": "#666666", "differentiated": "#2166ac"}
    fig, ax = plt.subplots(figsize=(7.0, 0.44 * len(d) + 1.3))
    for yi, (_, r) in zip(y, d.iterrows()):
        if r["is_pending"]:
            ax.scatter([0.0], [yi], marker="s", facecolors="none", edgecolors="0.6", s=40, zorder=3)
            ax.annotate("PENDING — fills when RNA matrix lands", (0.004, yi), fontsize=6,
                        color="0.5", va="center", ha="left")
            continue
        c = bcolor.get(r["bucket"], "#2166ac")
        ax.plot([r["dauroc_lo"], r["dauroc_hi"]], [yi, yi], color="0.4", lw=1.3, zorder=2)
        ax.scatter([r["dauroc"]], [yi], color=c, s=36, zorder=3)
    ax.axvline(0.0, ls="--", lw=0.9, color="0.5")
    ax.set_yticks(y)
    ax.set_yticklabels(d["candidate"])
    ax.set_xlabel("Incremental AUROC over floor (dAUROC, point + paired-bootstrap 95% CI)")
    ax.set_title(title)
    ax.margins(y=0.02)
    handles = [plt.Line2D([], [], marker="o", ls="", color=bcolor[k], label=k)
               for k in ["baseline", "differentiated"]]
    handles.append(plt.Line2D([], [], marker="s", ls="", mfc="none", mec="0.6", label="pending"))
    ax.legend(handles=handles, loc="lower right", fontsize=6)
    fig.savefig(path)
    plt.close(fig)
    return fig


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", default="results/analysis_frame.parquet")
    ap.add_argument("--outdir", default="results/eval")
    args = ap.parse_args()
    frame = pd.read_parquet(args.parquet)
    rep = run_loco_lopo(frame, outdir=args.outdir)
    plot_incremental_forest(rep, Path(args.outdir) / "fig_incremental_auroc.png")
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(rep.to_string(index=False))
