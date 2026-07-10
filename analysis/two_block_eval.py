#!/usr/bin/env python3
"""Pre-registered two-block evaluation: non-reference RNA features vs the immune floor.

Implements docs/EVAL_PROTOCOL.md EXACTLY. Run only AFTER the non-ref matrix exists.
Three models per split: (A) floor, (B) nonref, (C) floor+nonref.
CV frame auto-selected per the protocol: LOCO if >=2 cohorts w/ >=10 & both classes, else
grouped-by-patient 5-fold within the dominant cohort. Fold-internal standardization only.

Usage:
  python analysis/two_block_eval.py --nonref <matrix.parquet> --out results/eval
"""
import argparse, json, os, sys
import numpy as np, pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score

ERS = "/Users/alex/OrchestratedBiosciences/evolutionary-rna-state"
FLOOR_COLS = ["gep_tcell_inflamed", "ifng_score", "teff", "tgfb", "teff_tgfb_balance"]

def _oof_auroc(X, y, groups, splits, seed=0):
    """Out-of-fold AUROC, fold-internal scaling + L2 logit. splits: list of (tr_idx,te_idx)."""
    oof = np.full(len(y), np.nan)
    for tr, te in splits:
        if len(np.unique(y[tr])) < 2:
            continue
        scl = StandardScaler().fit(X[tr])
        Xtr, Xte = scl.transform(X[tr]), scl.transform(X[te])
        oof[te] = LogisticRegression(max_iter=2000, C=0.5).fit(Xtr, y[tr]).predict_proba(Xte)[:, 1]
    m = ~np.isnan(oof)
    return roc_auc_score(y[m], oof[m]) if len(np.unique(y[m])) > 1 else np.nan, oof

def _loco_splits(cohort):
    idx = np.arange(len(cohort))
    return [(idx[cohort != c], idx[cohort == c]) for c in pd.unique(cohort)]

def _grouped_splits(groups, y, seed=0, k=5):
    gkf = GroupKFold(n_splits=min(k, len(np.unique(groups))))
    return list(gkf.split(np.zeros(len(y)), y, groups))

def _hanley_mcneil_ci(auc, n_pos, n_neg, z=1.96):
    if np.isnan(auc) or n_pos == 0 or n_neg == 0:
        return (np.nan, np.nan)
    q1 = auc / (2 - auc); q2 = 2 * auc**2 / (1 + auc)
    se = np.sqrt((auc*(1-auc) + (n_pos-1)*(q1-auc**2) + (n_neg-1)*(q2-auc**2)) / (n_pos*n_neg))
    return (max(0, auc - z*se), min(1, auc + z*se))

def _perm_p(build_X, y, groups, splits, obs, n_perm=5000, seed=0):
    rng = np.random.default_rng(seed)
    # permute WITHIN cohort/group blocks (exchangeability)
    ge = 0
    for _ in range(n_perm):
        yp = y.copy()
        for g in np.unique(groups):
            mask = groups == g
            yp[mask] = rng.permutation(y[mask])
        a, _ = _oof_auroc(build_X, yp, groups, splits, seed)
        if not np.isnan(a) and a >= obs:
            ge += 1
    return (ge + 1) / (n_perm + 1)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nonref", required=True, help="non-ref feature matrix parquet (index or col run_accession)")
    ap.add_argument("--frame", default=None, help="label frame; default reconciled_frame_n106 artifact path")
    ap.add_argument("--floor", default=f"{ERS}/results/predictor/immune_floor_block.parquet")
    ap.add_argument("--out", default=f"{ERS}/results/eval")
    ap.add_argument("--n_perm", type=int, default=5000)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    # ---- assemble ----
    nb = pd.read_parquet(args.nonref)
    if "run_accession" not in nb.columns:
        nb = nb.reset_index().rename(columns={nb.index.name or "index": "run_accession"})
    fl = pd.read_parquet(args.floor)
    frame = pd.read_parquet(args.frame) if args.frame else pd.read_parquet(
        os.environ.get("N106_FRAME", f"{ERS}/results/predictor/frozen_analysis_set.parquet"))
    keep = ["run_accession", "cohort", "resp" if "resp" in frame.columns else "y"]
    if "patient_id" in frame.columns: keep.append("patient_id")
    lab = frame[keep].rename(columns={"resp": "y"})
    df = lab.merge(fl, on="run_accession").merge(nb, on="run_accession")
    df = df.dropna(subset=["y"]).reset_index(drop=True)
    df["y"] = df["y"].astype(int)
    if "patient_id" not in df.columns: df["patient_id"] = df["run_accession"]

    nonref_cols = [c for c in nb.columns if c != "run_accession"]
    y = df["y"].values; cohort = df["cohort"].values; patient = df["patient_id"].values
    Xf = df[FLOOR_COLS].apply(pd.to_numeric, errors="coerce").fillna(0).values
    Xn = df[nonref_cols].apply(pd.to_numeric, errors="coerce").fillna(0).values
    Xc = np.hstack([Xf, Xn])

    # ---- CV frame decision (per protocol, BEFORE looking at AUROCs) ----
    coh_ok = [c for c in pd.unique(cohort)
              if (cohort == c).sum() >= 10 and len(np.unique(y[cohort == c])) == 2]
    if len(coh_ok) >= 2:
        frame_kind = "LOCO"; splits = _loco_splits(cohort); groups = cohort
    else:
        dom = pd.Series(cohort).value_counts().idxmax()
        m = cohort == dom
        frame_kind = f"grouped5fold_within_{dom}"
        y, cohort, patient = y[m], cohort[m], patient[m]
        Xf, Xn, Xc = Xf[m], Xn[m], Xc[m]
        splits = _grouped_splits(patient, y); groups = patient

    npos, nneg = int((y == 1).sum()), int((y == 0).sum())
    res = {"frame": frame_kind, "n": int(len(y)), "n_pos": npos, "n_neg": nneg,
           "cohorts": {c: int((cohort == c).sum()) for c in pd.unique(cohort)},
           "nonref_n_features": len(nonref_cols), "blocks": {}}
    for name, X in [("A_floor", Xf), ("B_nonref", Xn), ("C_floor_plus_nonref", Xc)]:
        auc, _ = _oof_auroc(X, y, groups, splits)
        ci = _hanley_mcneil_ci(auc, npos, nneg)
        p = _perm_p(X, y, groups, splits, auc, n_perm=args.n_perm) if not np.isnan(auc) else np.nan
        res["blocks"][name] = {"auroc": None if np.isnan(auc) else round(float(auc), 4),
                               "ci95": [None if np.isnan(x) else round(float(x), 4) for x in ci],
                               "perm_p": None if np.isnan(p) else round(float(p), 4)}
    a = res["blocks"]["A_floor"]["auroc"]; c = res["blocks"]["C_floor_plus_nonref"]["auroc"]
    res["delta_C_minus_A"] = None if (a is None or c is None) else round(c - a, 4)
    outp = os.path.join(args.out, f"nonref_vs_floor_{frame_kind.split('_')[0].lower()}.json")
    json.dump(res, open(outp, "w"), indent=2)
    print(json.dumps(res, indent=2)); print("wrote", outp)

if __name__ == "__main__":
    main()
