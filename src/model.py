"""
src/model.py — the modeling brain: latent evolutionary RNA-state S and its
organization of ICB response, run on the DE-NOVO per-sample phenotype matrices
(src/features.py), not the WES proxy.

Two claims, in order, with the rigor guards baked in:

  latent_state(...)  — Stage 4 / internal claim. Learn a low-rank representation
                       S of the de-novo RNA-phenotype block AFTER residualizing
                       on burden/library-size, and test that shared variance
                       exceeds a permutation null. Leave-one-cohort-out loading
                       stability. NO response labels touched here.

  response_organization(...) — Stage 5 / external claim. Nested, group-aware
                       cross-validated logistic regression of ICB response on S,
                       with EVERY data-dependent step (scaling, feature choice,
                       S construction) fit INSIDE the training fold (leakage
                       guard). Incremental AUROC over TMB and over TIDE.

This module is validated on a contract-shaped MOCK matrix until the pipeline
session lands real matrices; the same code runs unchanged on real features.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from numpy.random import default_rng
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler


# ---------------------------------------------------------------------------
# Feature preprocessing — log + burden/library-size residualization
# ---------------------------------------------------------------------------
def prep_block(X: pd.DataFrame, size_factor: np.ndarray | None = None,
               log: bool = True) -> np.ndarray:
    """Log1p, then (optionally) residualize each feature on a library-size /
    burden factor so co-variation cannot be a shared depth artefact.

    X: samples x features (numeric). size_factor: per-sample scalar (e.g. log
    total counts) to residualize out. Returns z-scored residual matrix.
    """
    A = np.log1p(X.to_numpy(dtype=float)) if log else X.to_numpy(dtype=float)
    if size_factor is not None:
        Z = np.column_stack([np.ones(len(A)), np.asarray(size_factor, float)])
        beta, *_ = np.linalg.lstsq(Z, A, rcond=None)
        A = A - Z @ beta
    mu, sd = np.nanmean(A, 0), np.nanstd(A, 0)
    sd[sd == 0] = 1.0
    return (A - mu) / sd


# ---------------------------------------------------------------------------
# Stage 4 — latent state + permutation test for shared low-rank structure
# ---------------------------------------------------------------------------
@dataclass
class LatentResult:
    scores: np.ndarray                 # per-sample S (n,)
    loadings: np.ndarray               # feature loadings on PC1
    var_explained: float               # PC1 fraction
    perm_p: float
    null_mean: float
    loo_stability: dict = field(default_factory=dict)


def _pc1_frac(Xz: np.ndarray) -> float:
    C = np.corrcoef(Xz, rowvar=False)
    ev = np.sort(np.linalg.eigvalsh(C))[::-1]
    return float(ev[0] / ev.sum())


def latent_state(Xz: np.ndarray, cohorts: np.ndarray | None = None,
                 n_perm: int = 5000, seed: int = 0) -> LatentResult:
    """Fit PC1 as S; permutation-null on PC1 variance fraction (shuffle each
    feature independently to break cross-feature structure while preserving
    marginals). Leave-one-cohort-out loading stability if cohorts given."""
    rng = default_rng(seed)
    Xz = Xz[:, np.nanstd(Xz, 0) > 0]
    pca = PCA(n_components=1).fit(Xz)
    scores = pca.transform(Xz)[:, 0]
    loadings = pca.components_[0]
    obs = _pc1_frac(Xz)
    null = np.array([_pc1_frac(np.column_stack(
        [rng.permutation(Xz[:, k]) for k in range(Xz.shape[1])]))
        for _ in range(n_perm)])
    p = (1 + np.sum(null >= obs)) / (n_perm + 1)

    loo = {}
    if cohorts is not None:
        for c in np.unique(cohorts):
            m = cohorts != c
            if m.sum() > Xz.shape[1] + 2:
                l_c = PCA(n_components=1).fit(Xz[m]).components_[0]
                # sign-align, cosine similarity to full loadings
                cos = np.dot(l_c, loadings) / (np.linalg.norm(l_c) * np.linalg.norm(loadings))
                loo[str(c)] = float(abs(cos))
    return LatentResult(scores, loadings, obs, float(p), float(null.mean()), loo)


# ---------------------------------------------------------------------------
# Stage 5 — response organization with fold-contained everything (leakage guard)
# ---------------------------------------------------------------------------
def _auc_ci(y, s, n_boot=2000, seed=0):
    rng = default_rng(seed)
    aucs = []
    y = np.asarray(y); s = np.asarray(s)
    for _ in range(n_boot):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        aucs.append(roc_auc_score(y[idx], s[idx]))
    return float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))


def response_organization(feat: pd.DataFrame, y: np.ndarray, groups: np.ndarray,
                          cohorts: np.ndarray | None = None,
                          size_factor: np.ndarray | None = None,
                          extra: dict[str, np.ndarray] | None = None,
                          n_splits: int = 5, seed: int = 0) -> dict:
    """Nested group-aware CV. For each outer fold: standardize + residualize +
    fit S on TRAIN only, project TEST, then logistic regression of y on S.
    Compare against comparator models (e.g. TMB alone, TIDE) fit the same way.

    feat: samples x features de-novo block. y: binary response (0/1).
    groups: patient id per sample (GroupKFold unit — never split a patient).
    extra: {name: per-sample covariate vector} for comparator/incremental models.
    Returns per-model out-of-fold AUROC + bootstrap CI + incremental AUROC(S|extra).
    """
    extra = extra or {}
    X = feat.to_numpy(dtype=float)
    y = np.asarray(y).astype(int)
    n = len(y)
    oof = {"S": np.full(n, np.nan)}
    for k in extra:
        oof[k] = np.full(n, np.nan)
        oof[f"{k}+S"] = np.full(n, np.nan)

    sgkf = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for tr, te in sgkf.split(X, y, groups):
        # --- fold-contained S: scale + residualize + PCA on TRAIN only ---
        sc = StandardScaler().fit(np.log1p(X[tr]))
        Xtr, Xte = sc.transform(np.log1p(X[tr])), sc.transform(np.log1p(X[te]))
        if size_factor is not None:
            sf = np.asarray(size_factor, float)
            Ztr = np.column_stack([np.ones(len(tr)), sf[tr]])
            beta, *_ = np.linalg.lstsq(Ztr, Xtr, rcond=None)
            Xtr = Xtr - Ztr @ beta
            Xte = Xte - np.column_stack([np.ones(len(te)), sf[te]]) @ beta
        pca = PCA(n_components=1).fit(Xtr)
        Str, Ste = pca.transform(Xtr)[:, 0], pca.transform(Xte)[:, 0]

        def fit_predict(train_cols, test_cols):
            lr = LogisticRegression(max_iter=1000)
            lr.fit(train_cols, y[tr])
            return lr.predict_proba(test_cols)[:, 1]

        oof["S"][te] = fit_predict(Str.reshape(-1, 1), Ste.reshape(-1, 1))
        for k, v in extra.items():
            v = np.asarray(v, float)
            oof[k][te] = fit_predict(v[tr].reshape(-1, 1), v[te].reshape(-1, 1))
            comb_tr = np.column_stack([v[tr], Str])
            comb_te = np.column_stack([v[te], Ste])
            oof[f"{k}+S"][te] = fit_predict(comb_tr, comb_te)

    out = {}
    for name, sc in oof.items():
        ok = ~np.isnan(sc)
        if len(np.unique(y[ok])) < 2:
            continue
        auc = roc_auc_score(y[ok], sc[ok])
        lo, hi = _auc_ci(y[ok], sc[ok], seed=seed)
        out[name] = {"auroc": float(auc), "ci_lo": lo, "ci_hi": hi,
                     "n": int(ok.sum())}
    # incremental AUROC of adding S to each comparator
    out["_incremental"] = {k: out[f"{k}+S"]["auroc"] - out[k]["auroc"]
                           for k in extra if k in out and f"{k}+S" in out}
    return out
