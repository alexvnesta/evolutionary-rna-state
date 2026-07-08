"""
analysis/pilot/deepen_analysis.py — turn the n=12 pilot into a defensible
result on the expanded de-novo Salmon feature set.

Three additions the pilot flagged as needed:

  1. SCALE — quantify the full Gide 2019 PRE set (not just 12) so the
     antigen-presentation axis is estimated on adequate n.
  2. TUMOR-INTRINSIC vs INFILTRATION — the immune axis could be a proxy for
     'more immune cells'. We estimate an immune-infiltration score from the
     SAME de-novo data (leukocyte/stromal marker panel) and test whether the
     antigen-presentation axis still separates response AFTER residualizing
     infiltration out. If it does, the signal is not merely infiltration.
  3. HELD-OUT VALIDATION — fit the response model on Gide, evaluate on Riaz
     (a completely separate cohort), which is the honest external test.

Consumes: results/{pilot_salmon,expand_salmon}/<run>/quant.sf
          refs/tx2gene.tsv, data/manifests/selection_manifest.csv
Produces: results/deepen_features.csv, results/deepen_results.json,
          results/fig_deepen.png
"""
from __future__ import annotations

import glob
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).resolve().parents[2]

# curated Ensembl gene panels (version-stripped ids)
ANTIGEN = {"B2M": "ENSG00000166710", "HLA_A": "ENSG00000206503",
           "TAP1": "ENSG00000168394", "TAP2": "ENSG00000204267",
           "HLA_B": "ENSG00000234745", "PSMB9": "ENSG00000240065"}
CYTOLYTIC = {"GZMA": "ENSG00000145649", "PRF1": "ENSG00000180644"}
IFN = {"IFIT1": "ENSG00000185745", "OAS1": "ENSG00000089127",
       "MX1": "ENSG00000157601", "ISG15": "ENSG00000187608",
       "DDX58": "ENSG00000107201", "IFIH1": "ENSG00000115267",
       "STAT1": "ENSG00000115415", "CXCL10": "ENSG00000169245"}
# pan-leukocyte / infiltration markers (the confound to residualize)
INFILTRATION = {"PTPRC": "ENSG00000081237", "CD3D": "ENSG00000167286",
                "CD3E": "ENSG00000198851", "CD8A": "ENSG00000153563",
                "CD4": "ENSG00000010610", "CD2": "ENSG00000116824",
                "CD68": "ENSG00000129226", "ITGAM": "ENSG00000169896"}


def load_gene_tpm(dirs: list[Path], tx2gene: Path) -> pd.DataFrame:
    """Load per-sample transcript TPM from all quant dirs, collapse to gene."""
    m = pd.read_csv(tx2gene, sep="\t", header=None,
                    names=["tx", "gene"]).set_index("tx")["gene"]
    cols = {}
    for d in dirs:
        for q in glob.glob(str(d / "*" / "quant.sf")):
            run = os.path.basename(os.path.dirname(q))
            cols[run] = pd.read_csv(q, sep="\t").set_index("Name")["TPM"]
    tpm = pd.DataFrame(cols)
    g = tpm.join(m).groupby("gene").sum(numeric_only=True)
    g.index = g.index.str.split(".").str[0]
    return g


def axis_score(g: pd.DataFrame, panel: dict) -> pd.Series:
    """Mean z of log-TPM across a gene panel (per sample)."""
    log_g = np.log1p(g)
    present = [eid for eid in panel.values() if eid in g.index]
    sub = log_g.loc[present].T  # sample x gene
    z = StandardScaler().fit_transform(sub)
    return pd.Series(z.mean(axis=1), index=g.columns)


def residualize(target: pd.Series, covar: pd.Series) -> pd.Series:
    """Return target with the linear effect of covar removed (OLS residuals)."""
    X = np.column_stack([np.ones(len(covar)), covar.values])
    beta, *_ = np.linalg.lstsq(X, target.values, rcond=None)
    return pd.Series(target.values - X @ beta, index=target.index)


def loo_auroc(x: np.ndarray, y: np.ndarray) -> float:
    """LOO on a PRECOMPUTED scalar axis (logistic map only). Use
    foldcontained_loo when the axis involves standardization to avoid leakage."""
    x = np.asarray(x).reshape(-1, 1)
    oof = np.full(len(y), np.nan)
    for tr, te in LeaveOneOut().split(x):
        lr = LogisticRegression(max_iter=1000).fit(x[tr], y[tr])
        oof[te] = lr.predict_proba(x[te])[:, 1]
    return roc_auc_score(y, oof)


def foldcontained_loo(log_panel: pd.DataFrame, y, panel_ids,
                      covar_log_ids=None):
    """Leakage-free LOO for a z-scored gene-panel axis.

    log_panel: sample x gene log-TPM (raw, unscaled). Within each fold, the
    panel mean/std are estimated on TRAIN samples only, the held-out sample is
    projected with those train statistics, then a logistic map is fit on train.
    If covar_log_ids is given, the axis is residualized on that (train-fit)
    covariate axis inside the fold too (tumor-intrinsic test, leakage-free).
    Returns (auroc, oof_probs).
    """
    ids = [i for i in panel_ids if i in log_panel.columns]
    X = log_panel[ids].values
    y = np.asarray(y); n = len(y)
    C = log_panel[[i for i in (covar_log_ids or []) if i in log_panel.columns]].values \
        if covar_log_ids else None
    oof = np.full(n, np.nan)
    for tr, te in LeaveOneOut().split(X):
        mu, sd = X[tr].mean(0), X[tr].std(0); sd[sd == 0] = 1
        a_tr = ((X[tr] - mu) / sd).mean(1)
        a_te = ((X[te] - mu) / sd).mean(1)
        if C is not None:
            cmu, csd = C[tr].mean(0), C[tr].std(0); csd[csd == 0] = 1
            c_tr = ((C[tr] - cmu) / csd).mean(1)
            c_te = ((C[te] - cmu) / csd).mean(1)
            b = np.polyfit(c_tr, a_tr, 1)          # train-fit residualization
            a_tr = a_tr - np.polyval(b, c_tr)
            a_te = a_te - np.polyval(b, c_te)
        lr = LogisticRegression(max_iter=1000).fit(a_tr.reshape(-1, 1), y[tr])
        oof[te] = lr.predict_proba(a_te.reshape(-1, 1))[:, 1]
    return roc_auc_score(y, oof), oof


def perm_p(x: np.ndarray, y: np.ndarray, n_perm=5000, seed=0) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    obs = roc_auc_score(y, x); obs = max(obs, 1 - obs)
    null = np.array([max(a, 1 - a) for a in
                     (roc_auc_score(rng.permutation(y), x) for _ in range(n_perm))])
    return obs, float((1 + (null >= obs).sum()) / (n_perm + 1))


def build(dirs, manifest_paths):
    tx2gene = REPO / "refs" / "tx2gene.tsv"
    g = load_gene_tpm([Path(d) for d in dirs], tx2gene)
    man = pd.concat([pd.read_csv(p) for p in manifest_paths]).drop_duplicates("run_accession")
    man = man.set_index("run_accession")
    runs = [r for r in g.columns if r in man.index]
    g = g[runs]
    meta = man.loc[runs, ["cohort", "resp_NR", "treatment_arm", "recist"]].copy()
    meta["y"] = (meta.resp_NR.str.upper() == "R").astype(int)
    F = pd.DataFrame({
        "antigen_presentation": axis_score(g, ANTIGEN),
        "cytolytic": axis_score(g, CYTOLYTIC),
        "viral_mimicry_IFN": axis_score(g, IFN),
        "infiltration": axis_score(g, INFILTRATION),
    }, index=g.columns).loc[runs]
    # infiltration-residualized antigen presentation (tumor-intrinsic test)
    F["antigen_presentation_resid"] = residualize(F["antigen_presentation"], F["infiltration"])
    return g, meta, F
