"""
run_activity_response_tests.py
==============================================================================
Move-2 harness (see HANDOFF_rna_state_next.md): run the regulator-activity tests
once the full 5-cohort expression matrix lands. Reuses the validated gide2019
crosswalk (gide2019_id_crosswalk.csv) to join expression-keyed activity scores
(run_accession) to the iAtlas burden/response matrix (analysis_frame.parquet).

Three tests, faithful to io_rna_antigen_analysis.md / io_rbp_and_clonality.md:

  T1  MECHANISTIC   regulator activity vs matched-source antigen burden (Spearman)
                    SPLICING_FACTOR~SPLICE, ADAR_EDITING~ERV, RBP_BROAD~FUSION,
                    plus off-target and *~TMB controls.
  T2  SHARED-STATE  PC1 %-variance of the burden-adjusted activity block vs a
                    permutation null (residualize on log-TMB, log-mutation-count,
                    cohort; 5000 column-permutations). Mirrors the burden-feature
                    covariation test, now run on ACTIVITY scores.
  T3  PREDICTION    leave-one-cohort-out (LOCO) logistic-regression AUROC:
                    TMB alone / activity alone / TMB+activity. The bar is the
                    TMB floor (~0.62); activity must add over it.

Run now (pilot smoke test, n=12, 1 cohort -> T2/T3 degrade gracefully):
    python run_activity_response_tests.py \
        --activity rbp_activity_pilot_scores.csv \
        --crosswalk gide2019_id_crosswalk.csv \
        --frame analysis_frame.parquet --out results_pilot

Run at scale (all 5 cohorts): pass the full activity matrix. If you have a TPM
matrix + HGNC->Ensembl map instead of scores, call activity_from_tpm() first.

Activity input contract: a CSV/parquet with a run_accession column (or index)
and one column per regulator set (SPLICING_FACTOR, RBP_BROAD, ADAR_EDITING, ...).
"""
from __future__ import annotations
import argparse, json, sys
import numpy as np
import pandas as pd
from scipy import stats

# matched regulator -> antigen source for the mechanistic test
MATCHED = [
    ("SPLICING_FACTOR", "SPLICE_NEOANTIGEN"),
    ("ADAR_EDITING",    "ERV_NEOANTIGEN"),
    ("RBP_BROAD",       "FUSION_NEOANTIGEN"),
]
# off-target controls (regulator vs a source it should NOT drive)
OFFTARGET = [
    ("SPLICING_FACTOR", "FUSION_NEOANTIGEN"),
    ("SPLICING_FACTOR", "ERV_NEOANTIGEN"),
]
BURDEN_COVARS = ["TMB_NONSYNONYMOUS", "MUTATION_COUNT"]


# ---------------------------------------------------------------------------
# optional: activity from a TPM matrix (wraps the project's validated scorer)
# ---------------------------------------------------------------------------
def activity_from_tpm(tpm_df, gene_id_map, sets=None, min_genes=3):
    """Thin wrapper over rbp_activity_scorer.score_regulator_activity.
    tpm_df: rows=samples(run_accession), cols=Ensembl gene IDs, linear TPM.
    gene_id_map: {HGNC_symbol: ENSG...}. Returns samples x set-scores."""
    from rbp_activity_scorer import score_regulator_activity
    return score_regulator_activity(tpm_df, gene_id_map, sets=sets, min_genes=min_genes)


# ---------------------------------------------------------------------------
# data assembly
# ---------------------------------------------------------------------------
def _read_any(path):
    return pd.read_parquet(path) if str(path).endswith(".parquet") else pd.read_csv(path)


def load_and_join(activity_path, crosswalk_path, frame_path):
    """Join activity(run_accession) -> crosswalk -> iAtlas burden/response.
    Returns (merged_df, set_cols, meta). Only pre-treatment rows with a burden
    row are retained (iatlas_burden_available)."""
    act = _read_any(activity_path)
    if "run_accession" not in act.columns:
        act = act.reset_index().rename(columns={act.index.name or "index": "run_accession"})
    # activity set columns = numeric cols that aren't the key or a stray label
    drop = {"run_accession", "resp", "RESPONDER"}
    set_cols = [c for c in act.columns
                if c not in drop and pd.api.types.is_numeric_dtype(act[c])]

    xw = _read_any(crosswalk_path)
    if "iatlas_burden_available" in xw.columns:
        xw = xw[xw["iatlas_burden_available"] == True]  # noqa: E712  (pre-treatment only)
    # required keys only; carry optional descriptors (arm/timepoint/cohort) if present
    keep = ["run_accession", "iatlas_patientId"] + [
        c for c in ("arm", "timepoint", "cohort") if c in xw.columns]
    xw = xw[keep].drop_duplicates("run_accession")

    frame = _read_any(frame_path)
    bcols = (["SPLICE_NEOANTIGEN", "ERV_NEOANTIGEN", "FUSION_NEOANTIGEN",
              "INDEL_NEOANTIGEN", "SNV_NEOANTIGEN"] + BURDEN_COVARS
             + ["cohort", "RESPONDER"])
    bcols = [c for c in bcols if c in frame.columns]
    frame = frame[["patientId"] + bcols].rename(columns={"patientId": "iatlas_patientId"})

    # iAtlas patient IDs (e.g. 'Pt1') are NOT unique across cohorts (hugo & riaz collide),
    # so join on (cohort, iatlas_patientId) whenever the crosswalk carries cohort.
    if "cohort" in xw.columns and "cohort" in frame.columns:
        m = act.merge(xw, on="run_accession", how="inner").merge(
            frame, on=["cohort", "iatlas_patientId"], how="inner")
    else:
        xw = xw.drop(columns=[c for c in ("cohort",) if c in xw.columns])
        m = act.merge(xw, on="run_accession", how="inner").merge(
            frame, on="iatlas_patientId", how="inner")
    meta = dict(n_activity=len(act), n_crosswalk=len(xw), n_joined=len(m),
                set_cols=set_cols,
                cohorts=sorted(m["cohort"].dropna().unique().tolist()) if "cohort" in m else [])
    return m, set_cols, meta


# ---------------------------------------------------------------------------
# T1  mechanistic: regulator activity vs matched-source burden
# ---------------------------------------------------------------------------
def test_mechanistic(m, set_cols):
    rows = []
    pairs = ([p for p in MATCHED if p[0] in set_cols]
             + [p for p in OFFTARGET if p[0] in set_cols]
             + [(s, "TMB_NONSYNONYMOUS") for s in set_cols])
    for reg, src in pairs:
        if reg not in m or src not in m:
            continue
        sub = m[[reg, src]].dropna()
        kind = ("matched" if (reg, src) in MATCHED
                else "control_tmb" if src == "TMB_NONSYNONYMOUS" else "offtarget")
        if sub[src].nunique() < 2 or sub[reg].nunique() < 2 or len(sub) < 4:
            rows.append(dict(regulator=reg, source=src, kind=kind, n=len(sub),
                             rho=np.nan, p=np.nan, note="constant/insufficient"))
            continue
        rho, p = stats.spearmanr(sub[reg], sub[src])
        rows.append(dict(regulator=reg, source=src, kind=kind, n=len(sub),
                         rho=round(float(rho), 4), p=round(float(p), 4), note=""))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# T2  shared-state: PC1 %-variance of burden-adjusted activity vs perm null
# ---------------------------------------------------------------------------
def _residualize(Y, covar_df):
    """OLS-residualize each column of Y on covariates (incl. intercept)."""
    X = np.column_stack([np.ones(len(covar_df))] + [covar_df[c].values for c in covar_df.columns])
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    return Y - X @ beta


def _pc1_frac(Z):
    Zc = Z - Z.mean(0)
    if Zc.shape[1] < 2:
        return np.nan
    s = np.linalg.svd(Zc, full_matrices=False)[1]
    return float((s[0] ** 2) / (s ** 2).sum())


def test_shared_state(m, set_cols, n_perm=5000, seed=0):
    use = [c for c in set_cols if c in m and m[c].notna().sum() > 0]
    if len(use) < 2:
        return dict(status="skipped", reason="need >=2 activity sets", pc1=None)
    # covariates: log-TMB, log-mutation-count (where present), cohort dummies
    cov = pd.DataFrame(index=m.index)
    if "TMB_NONSYNONYMOUS" in m and m["TMB_NONSYNONYMOUS"].notna().any():
        cov["logTMB"] = np.log1p(m["TMB_NONSYNONYMOUS"])
    if "MUTATION_COUNT" in m and m["MUTATION_COUNT"].notna().any():
        cov["logMUT"] = np.log1p(m["MUTATION_COUNT"])
    if "cohort" in m and m["cohort"].nunique() > 1:
        d = pd.get_dummies(m["cohort"], drop_first=True, prefix="ch").astype(float)
        cov = pd.concat([cov, d], axis=1)
    block = m[use].copy()
    keep = block.notna().all(1)
    if cov.shape[1]:
        keep &= cov.notna().all(1)
    block = block[keep]
    covk = cov[keep] if cov.shape[1] else None
    n = len(block)
    if n < len(use) + 2:
        return dict(status="skipped", reason=f"n={n} too small for {len(use)} sets",
                    pc1=None, n=n)
    Y = block.values.astype(float)
    R = _residualize(Y, covk) if (covk is not None and covk.shape[1]) else (Y - Y.mean(0))
    R = R / (R.std(0) + 1e-12)
    obs = _pc1_frac(R)
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for i in range(n_perm):
        Rp = np.column_stack([rng.permutation(R[:, j]) for j in range(R.shape[1])])
        null[i] = _pc1_frac(Rp)
    p = float((null >= obs).mean())
    return dict(status="ok", n=n, n_sets=len(use), sets=use,
                pc1_observed=round(obs, 4), null_mean=round(float(null.mean()), 4),
                null_sd=round(float(null.std()), 4), p_value=round(p, 4),
                covariates=list(covk.columns) if covk is not None else [])


# ---------------------------------------------------------------------------
# T3  prediction: leave-one-cohort-out logistic AUROC
# ---------------------------------------------------------------------------
def _loco_auroc(m, feat_cols):
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.metrics import roc_auc_score
    d = m.dropna(subset=feat_cols + ["RESPONDER", "cohort"]).copy()
    d["y"] = d["RESPONDER"].astype(int)
    cohorts = [c for c in d["cohort"].unique() if d[d.cohort == c]["y"].nunique() == 2]
    oof_y, oof_p = [], []
    for held in d["cohort"].unique():
        tr = d[d.cohort != held]
        te = d[d.cohort == held]
        if tr["y"].nunique() < 2 or len(te) == 0:
            continue
        sc = StandardScaler().fit(tr[feat_cols])
        clf = LogisticRegression(max_iter=1000, C=1.0)
        clf.fit(sc.transform(tr[feat_cols]), tr["y"])
        p = clf.predict_proba(sc.transform(te[feat_cols]))[:, 1]
        oof_y += te["y"].tolist(); oof_p += p.tolist()
    if len(set(oof_y)) < 2:
        return dict(auroc=None, n=len(oof_y), note="pooled OOF single-class (need >=2 cohorts)")
    return dict(auroc=round(float(roc_auc_score(oof_y, oof_p)), 4), n=len(oof_y), note="")


def test_prediction(m, set_cols):
    have_tmb = "TMB_NONSYNONYMOUS" in m and m["TMB_NONSYNONYMOUS"].notna().any()
    act = [c for c in set_cols if c in m]
    out = {}
    if have_tmb:
        out["TMB_alone"] = _loco_auroc(m[m["TMB_NONSYNONYMOUS"].notna()], ["TMB_NONSYNONYMOUS"])
        out["TMB_plus_activity"] = _loco_auroc(
            m[m["TMB_NONSYNONYMOUS"].notna()], ["TMB_NONSYNONYMOUS"] + act)
    out["activity_alone"] = _loco_auroc(m, act)
    n_coh = m["cohort"].nunique() if "cohort" in m else 1
    if n_coh < 2:
        out["_warning"] = f"only {n_coh} cohort(s) -> LOCO undefined; run at scale."
    return out


# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--activity", required=True)
    ap.add_argument("--crosswalk", required=True)
    ap.add_argument("--frame", required=True)
    ap.add_argument("--out", default="activity_test_results")
    ap.add_argument("--n-perm", type=int, default=5000)
    a = ap.parse_args(argv)

    m, set_cols, meta = load_and_join(a.activity, a.crosswalk, a.frame)
    t1 = test_mechanistic(m, set_cols)
    t2 = test_shared_state(m, set_cols, n_perm=a.n_perm)
    t3 = test_prediction(m, set_cols)

    t1.to_csv(f"{a.out}_T1_mechanistic.csv", index=False)
    report = dict(meta=meta, T2_shared_state=t2, T3_prediction=t3)
    with open(f"{a.out}_report.json", "w") as fh:
        json.dump(report, fh, indent=2, default=str)

    print(f"[join] activity={meta['n_activity']} crosswalk_pre={meta['n_crosswalk']} "
          f"joined={meta['n_joined']} cohorts={meta['cohorts']}")
    print("\n[T1] mechanistic regulator ~ matched-source burden")
    print(t1.to_string(index=False))
    print("\n[T2] shared-state PC1 vs permutation null")
    print(json.dumps(t2, indent=2, default=str))
    print("\n[T3] LOCO AUROC (bar = TMB floor ~0.62)")
    print(json.dumps(t3, indent=2, default=str))
    print(f"\n[out] {a.out}_T1_mechanistic.csv, {a.out}_report.json")
    return report


if __name__ == "__main__":
    main()
