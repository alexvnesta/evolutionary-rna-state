"""
analysis/01_covariation.py — Internal claim (falsifiable claim #1).

Do the mechanism-resolved RNA-phenotype neoantigen loads *co-vary*, and does
that shared structure survive removal of mutational burden? Tested WITHOUT any
reference to response labels.

Design honesty (set by the data, see results/coverage_matrix.csv):
  * The iAtlas proxy lights up only SPLICE (99% >0) and FUSION (68% >0) among
    the three thesis RNA-phenotypes; ERV/TE is degenerate (1.9% >0, max=2) and
    cannot carry a covariation signal from WES-derived calls. This is the
    ceiling of the proxy stage and the motivation for the raw-read arm.
  * Burden covariates (TMB, mutation count) exist only in Hugo/Liu/Riaz (192
    samples); Gide (72) has no WES. So the burden-adjusted test runs on the
    192-sample burden set; the raw association is shown on all 264 for context.

Proxy-circularity guard (harness test 1): every neoantigen category scales with
mutational/expression burden, so a raw SPLICE-FUSION correlation could be pure
burden. The claim requires the association to SURVIVE residualizing on burden
(TMB + mutation count) and on cohort (batch guard).

Outputs (results/):
  covariation_stats.csv   — raw & partial Spearman, per-cohort + pooled, with CIs
  covariation_perm.json   — permutation-null test of the low-rank component
  fig_covariation.png     — raw vs partial correlation heatmaps + scree/null
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results"

RNA = ["SPLICE_NEOANTIGEN", "ERV_NEOANTIGEN", "FUSION_NEOANTIGEN"]
MUT = ["INDEL_NEOANTIGEN", "SNV_NEOANTIGEN"]
BURDEN = ["TMB_NONSYNONYMOUS", "MUTATION_COUNT"]
LIVE = ["SPLICE_NEOANTIGEN", "FUSION_NEOANTIGEN"]   # non-degenerate RNA proxies
BURDEN_COHORTS = ["hugo2016", "liu2019", "riaz2017"]

RNG = np.random.default_rng(0)


# ---------------------------------------------------------------------------
# Partial Spearman: rank-transform, residualize on covariates (incl. intercept),
# correlate residuals. Covariate matrix Z may include burden + cohort dummies.
# ---------------------------------------------------------------------------
def _resid(y: np.ndarray, Z: np.ndarray) -> np.ndarray:
    Z1 = np.column_stack([np.ones(len(y)), Z]) if Z.size else np.ones((len(y), 1))
    beta, *_ = np.linalg.lstsq(Z1, y, rcond=None)
    return y - Z1 @ beta


def partial_spearman(df, a, b, covars, boot=2000, seed=0):
    """Partial Spearman between a and b controlling for covars (list of cols).

    Ranks all variables (Spearman = Pearson on ranks), residualizes a and b on
    ranked covariates, correlates residuals. Bootstrap CI + permutation p.
    """
    cols = [a, b] + covars
    d = df[cols].dropna()
    n = len(d)
    R = d.rank().to_numpy()
    ya, yb = R[:, 0], R[:, 1]
    Z = R[:, 2:] if len(covars) else np.empty((n, 0))
    ra, rb = _resid(ya, Z), _resid(yb, Z)
    rho = stats.pearsonr(ra, rb)[0]
    # bootstrap CI
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(boot):
        idx = rng.integers(0, n, n)
        Rb = d.iloc[idx].rank().to_numpy()
        Zb = Rb[:, 2:] if len(covars) else np.empty((n, 0))
        boots.append(stats.pearsonr(_resid(Rb[:, 0], Zb),
                                    _resid(Rb[:, 1], Zb))[0])
    lo, hi = np.percentile(boots, [2.5, 97.5])
    # permutation p (shuffle b's residuals)
    perm = []
    for _ in range(boot):
        perm.append(abs(stats.pearsonr(ra, rng.permutation(rb))[0]))
    p = (1 + np.sum(np.array(perm) >= abs(rho))) / (boot + 1)
    return dict(a=a, b=b, covars=";".join(covars) or "none",
                n=int(n), rho=float(rho), ci_lo=float(lo), ci_hi=float(hi),
                perm_p=float(p))


def cohort_dummies(df):
    return pd.get_dummies(df["cohort"], prefix="coh", drop_first=True).astype(float)


def run():
    RESULTS.mkdir(exist_ok=True)
    pre = pd.read_parquet(RESULTS / "analysis_frame.parquet")
    ph = pre[pre["in_phenotype_set"]].copy()
    for c in BURDEN:
        ph[c + "_log"] = np.log1p(ph[c])
    burden_log = [c + "_log" for c in BURDEN]
    ph_burden = ph[ph[BURDEN[0]].notna()].copy()          # 192, burden known

    rows = []

    # (1) raw SPLICE-FUSION, all phenotype cohorts (no burden control), pooled
    rows.append({"analysis": "raw_all_cohorts", "cohort": "pooled",
                 **partial_spearman(ph, *LIVE, [])})
    # per-cohort raw (replication guard)
    for coh, sub in ph.groupby("cohort"):
        if sub[LIVE].dropna().shape[0] >= 15:
            rows.append({"analysis": "raw_percohort", "cohort": coh,
                         **partial_spearman(sub, *LIVE, [])})

    # (2) burden-adjusted (proxy-circularity guard), burden cohorts only
    #     pooled: control burden + cohort dummies (batch guard)
    phb = ph_burden.copy()
    dums = cohort_dummies(phb)
    phb = pd.concat([phb, dums], axis=1)
    rows.append({"analysis": "burden_adj_pooled", "cohort": "pooled(H/L/R)",
                 **partial_spearman(phb, *LIVE, burden_log + list(dums.columns))})
    # per-cohort burden-adjusted
    for coh in BURDEN_COHORTS:
        sub = ph_burden[ph_burden["cohort"] == coh]
        if sub[LIVE].dropna().shape[0] >= 15:
            rows.append({"analysis": "burden_adj_percohort", "cohort": coh,
                         **partial_spearman(sub, *LIVE, burden_log)})

    stats_df = pd.DataFrame(rows)
    stats_df.to_csv(RESULTS / "covariation_stats.csv", index=False)

    # (3) full raw + partial correlation matrices (for the heatmap figure)
    feats = RNA + MUT
    raw_mat = ph[feats].corr(method="spearman")
    # partial matrix on burden set: each pair residualized on burden+cohort
    dums_b = cohort_dummies(ph_burden)
    phb2 = pd.concat([ph_burden, dums_b], axis=1)
    part_mat = pd.DataFrame(np.eye(len(feats)), index=feats, columns=feats)
    cov = burden_log + list(dums_b.columns)
    for i, fi in enumerate(feats):
        for fj in feats[i + 1:]:
            r = partial_spearman(phb2, fi, fj, cov, boot=500)["rho"]
            part_mat.loc[fi, fj] = part_mat.loc[fj, fi] = r
    raw_mat.to_csv(RESULTS / "covariation_raw_matrix.csv")
    part_mat.to_csv(RESULTS / "covariation_partial_matrix.csv")

    # (4) permutation-null low-rank test on burden-residualized LIVE block
    #     residualize each live proxy on burden+cohort, PCA, PC1 variance vs null
    Rr = ph_burden[LIVE + BURDEN].rank()
    dums_r = cohort_dummies(ph_burden)
    Zc = np.column_stack([Rr[BURDEN].to_numpy(), dums_r.to_numpy()])
    resid = np.column_stack([_resid(Rr[c].to_numpy(), Zc) for c in LIVE])
    resid = (resid - resid.mean(0)) / resid.std(0)

    def pc1_frac(X):
        C = np.corrcoef(X, rowvar=False)
        ev = np.sort(np.linalg.eigvalsh(C))[::-1]
        return ev[0] / ev.sum()

    obs = pc1_frac(resid)
    null = []
    for _ in range(5000):
        Xp = np.column_stack([RNG.permutation(resid[:, k])
                              for k in range(resid.shape[1])])
        null.append(pc1_frac(Xp))
    null = np.array(null)
    perm = {"metric": "PC1_variance_fraction_of_burden_residualized_LIVE_block",
            "features": LIVE, "n_samples": int(len(resid)),
            "obs_pc1_frac": float(obs), "null_mean": float(null.mean()),
            "null_p95": float(np.percentile(null, 95)),
            "perm_p": float((1 + np.sum(null >= obs)) / (len(null) + 1)),
            "n_perm": int(len(null))}
    (RESULTS / "covariation_perm.json").write_text(json.dumps(perm, indent=2))

    return stats_df, raw_mat, part_mat, perm


if __name__ == "__main__":
    s, raw, part, perm = run()
    pd.set_option("display.width", 200, "display.max_columns", 20)
    print(s.to_string(index=False))
    print("\nPERM:", json.dumps(perm, indent=2))
