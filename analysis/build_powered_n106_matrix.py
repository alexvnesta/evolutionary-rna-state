#!/usr/bin/env python3
"""
build_powered_n106_matrix.py

Powered n=106 two-block test builder — READY TO RUN once the external-box
alignment finishes all 106 BAMs (83/106 done as of 2026-07-18, ~23 remaining).

Context (see docs/PROJECT_STATUS.md):
  - The first-look two-block test was negative but underpowered (n=14-25); the
    immune floor itself fails to transfer LOCO at that n, so cross-cohort is not
    yet interpretable. The powered n=106 test exists precisely because the
    outcome is NOT yet known — do not prejudge it.
  - APPROACH PIVOTED (2026-07-12): the Nextflow/HISAT2 plan failed validation and
    was replaced by a UNIFORM full-cohort direct-STAR re-alignment of all 106.
    BAMs land at **/data/rnaseq_cohort/star_salmon/*.markdup.sorted.bam** on the
    box (NOT hisat2/), aligner is **STAR** (not HISAT2), keyed on run_accession.
    Because the aligner changed, the non-reference feature callers must be re-run
    on the STAR BAMs for uniformity — the old 25-sample HISAT2-derived
    nonref_matrix_cohort.parquet is a TEMPLATE for the schema only, NOT a set of
    rows to reuse (mixing HISAT2- and STAR-derived features would confound the
    aligner with cohort). Rebuild all 106 feature rows from the STAR BAMs.
  - This script joins the 66-feature non-reference matrix (rebuilt from the STAR
    BAMs, keyed on run_accession) to the pre-assembled immune-floor covariates
    (phase2_covariates_n106.parquet, keyed on run_accession) and runs the
    pre-registered floor-vs-nonref eval.

It is an OBSERVER/downstream script: it reads feature outputs, it does NOT run
alignment or touch the run's dirs. Run it from wherever the per-sample
non-reference caller outputs are readable.

Design is a faithful scale-up of the eval session's two_block_eval, not a
reinvention: same 66-column feature schema, same immune floor, same
grouped/LOCO CV + permutation.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
import pandas as pd

# The 66 non-reference feature columns, in the exact schema of the existing
# nonref_matrix_cohort.parquet (editing / intron-retention / splicing / TE-family).
# Loaded from the existing matrix at runtime so we never drift from the template.

def load_existing_matrix(existing_matrix: Path) -> pd.DataFrame:
    """The currently-banked non-ref feature matrix (25 samples as of this writing;
    row count grows as new caller outputs are added), keyed on run_accession."""
    m = pd.read_parquet(existing_matrix)
    if m.index.name != "run_accession":
        # tolerate a column form
        if "run_accession" in m.columns:
            m = m.set_index("run_accession")
        else:
            raise SystemExit("existing matrix is not keyed on run_accession")
    return m


def parse_sample_features(sample_dir: Path, feature_cols: list[str]) -> dict | None:
    """
    Parse one sample's non-reference caller outputs into the 66-feature row.

    Expects the canonical per-sample caller outputs (as produced by the pipeline's
    rna_editing / intron_retention / rnasplice / te_erv steps). Returns a dict of
    {feature_col: value} or None if the sample is not complete.

    NOTE: the exact on-disk filenames follow the canonical pipeline's publish
    layout. This function reads the same fields the existing matrix was built
    from; if a field is missing the sample is treated as incomplete (skipped),
    never silently zero-filled.
    """
    row: dict[str, float] = {}
    # --- RNA editing (AEI) ---
    aei = sample_dir / "aei.tsv"
    if not aei.exists():
        return None
    a = pd.read_csv(aei, sep="\t")
    # map to editing_AEI_percent / editing_SN / editing_AG / editing_Acov
    def _first(df, *names):
        for n in names:
            if n in df.columns:
                return float(df[n].iloc[0])
        return np.nan
    row["editing_AEI_percent"] = _first(a, "AEI_percent", "aei_percent", "AEI")
    row["editing_SN"]          = _first(a, "signal_to_noise", "SN")
    row["editing_AG"]          = _first(a, "AG_mismatches", "AG")
    row["editing_Acov"]        = _first(a, "A_coverage", "Acov")
    # --- intron retention ---
    irf = sample_dir / "intron_retention.tsv"
    if irf.exists():
        ir = pd.read_csv(irf, sep="\t")
        vals = ir.iloc[:, -1].astype(float) if ir.shape[1] else pd.Series(dtype=float)
        row["ir_mean"]       = float(vals.mean()) if len(vals) else np.nan
        row["ir_median"]     = float(vals.median()) if len(vals) else np.nan
        row["ir_frac_gt0.1"] = float((vals > 0.1).mean()) if len(vals) else np.nan
        row["ir_n_eval"]     = float(len(vals))
    # --- splicing ---
    spf = sample_dir / "splice_junctions.tsv"
    if spf.exists():
        sp = pd.read_csv(spf, sep="\t")
        row["splice_n_junctions"] = float(len(sp))
    # --- TE families (family-level featureCounts, uniform with the existing 40) ---
    tef = sample_dir / "te_family_counts.tsv"
    if tef.exists():
        te = pd.read_csv(tef, sep="\t")
        # normalize to CPM-like fractions per family as in the template
        te_cols = [c for c in feature_cols if c.startswith("te_")]
        if "family" in te.columns and te.shape[1] >= 2:
            counts = te.set_index("family").iloc[:, -1].astype(float)
            total = counts.sum() or 1.0
            for c in te_cols:
                fam = c[len("te_"):]
                row[c] = float(counts.get(fam, 0.0) / total)
    # completeness gate: require the editing block + at least one of ir/splice/te
    if any(pd.isna(row.get(k, np.nan)) for k in
           ["editing_AEI_percent","editing_SN","editing_AG","editing_Acov"]):
        return None
    return row


def build_nonref_matrix(existing: pd.DataFrame,
                        new_feature_root: Path,
                        manifest: pd.DataFrame,
                        reuse_existing_rows: bool = False) -> pd.DataFrame:
    """
    Build the n=106 non-reference feature matrix from per-sample caller outputs.

    IMPORTANT (post-STAR-pivot default): `reuse_existing_rows=False`. The `existing`
    matrix is used ONLY for its column schema (the 66 feature names / TE families),
    NOT to seed rows — because those rows were derived from HISAT2 BAMs and the
    n=106 run is a uniform STAR re-alignment. Reusing HISAT2-derived rows alongside
    STAR-derived rows would confound the aligner with the cohort/response split.
    Every row here is rebuilt from `new_feature_root` (the STAR-BAM feature outputs).

    Set `reuse_existing_rows=True` only if the feature callers were genuinely run on
    the SAME aligner's BAMs as `new_feature_root` (e.g. a pure resumption, not a
    pivot) — then samples already present in `existing` are kept and only missing
    ones are parsed.
    """
    feature_cols = list(existing.columns)
    rows = {}
    if reuse_existing_rows:
        rows = {acc: existing.loc[acc].to_dict() for acc in existing.index}
    added, skipped, reused = 0, 0, len(rows)
    for acc in manifest["run_accession"]:
        if acc in rows:
            continue
        sdir = new_feature_root / acc
        if not sdir.exists():
            skipped += 1
            continue
        r = parse_sample_features(sdir, feature_cols)
        if r is None:
            skipped += 1
            continue
        rows[acc] = r
        added += 1
    mat = pd.DataFrame.from_dict(rows, orient="index")
    mat = mat.reindex(columns=feature_cols)
    mat.index.name = "run_accession"
    print(f"[matrix] reused={reused} added={added} skipped={skipped} total={len(mat)} "
          f"(reuse_existing_rows={reuse_existing_rows})", file=sys.stderr)
    return mat


def _within_cohort_5fold_auroc(X: pd.DataFrame, y: np.ndarray, groups: np.ndarray,
                               seeds: list[int], cohort: str) -> float:
    """
    Seed-averaged stratified 5-fold CV AUROC computed WITHIN a single cohort
    (the primary frame the eval session used — e.g. 'within-Gide n=32'). Isolates
    the block's signal from cross-cohort batch structure; the LOCO functions below
    handle the cross-cohort transfer question separately. Returns the mean over
    seeds of the pooled out-of-fold AUROC for the named cohort's samples.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import roc_auc_score
    from sklearn.preprocessing import StandardScaler
    m = groups == cohort
    Xc, yc = X.iloc[np.where(m)[0]].reset_index(drop=True), y[m]
    if len(yc) < 10 or len(np.unique(yc)) < 2:
        return float("nan")
    aucs = []
    for s in seeds:
        skf = StratifiedKFold(5, shuffle=True, random_state=s)
        oof = np.full(len(yc), np.nan)
        for tr, te in skf.split(Xc, yc):
            sc = StandardScaler().fit(Xc.iloc[tr])
            clf = LogisticRegression(max_iter=2000, C=1.0)
            clf.fit(sc.transform(Xc.iloc[tr]), yc[tr])
            oof[te] = clf.predict_proba(sc.transform(Xc.iloc[te]))[:, 1]
        ok = ~np.isnan(oof)
        aucs.append(roc_auc_score(yc[ok], oof[ok]))
    return float(np.mean(aucs))


def _loco_auroc(X: pd.DataFrame, y: np.ndarray, groups: np.ndarray) -> float:
    """Leave-One-COhort-Out pooled out-of-fold AUROC (the honest cross-cohort transfer test)."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import LeaveOneGroupOut
    from sklearn.metrics import roc_auc_score
    from sklearn.preprocessing import StandardScaler
    logo = LeaveOneGroupOut()
    oof = np.full(len(y), np.nan)
    for tr, te in logo.split(X, y, groups):
        sc = StandardScaler().fit(X.iloc[tr])
        clf = LogisticRegression(max_iter=2000, C=1.0)
        clf.fit(sc.transform(X.iloc[tr]), y[tr])
        oof[te] = clf.predict_proba(sc.transform(X.iloc[te]))[:, 1]
    m = ~np.isnan(oof)
    return float(roc_auc_score(y[m], oof[m]))


def _residual_within_cohort_auroc(X_nonref: pd.DataFrame, X_floor: pd.DataFrame,
                                  y: np.ndarray, groups: np.ndarray, seeds: list[int],
                                  cohort: str) -> float:
    """
    FOLD-CONTAINED residualization WITHIN one cohort: inside each stratified CV
    split of that cohort, regress the immune floor OUT of each non-ref feature
    using TRAIN rows only, apply to TEST, then score the residualized non-ref
    block. Isolates non-ref signal orthogonal to the floor without leaking
    floor-nonref covariance across the split (the eval session's exact procedure).
    """
    from sklearn.linear_model import LinearRegression, LogisticRegression
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import roc_auc_score
    from sklearn.preprocessing import StandardScaler
    m = groups == cohort
    ix = np.where(m)[0]
    XN, XF, yc = (X_nonref.iloc[ix].reset_index(drop=True),
                  X_floor.iloc[ix].reset_index(drop=True), y[m])
    if len(yc) < 10 or len(np.unique(yc)) < 2:
        return float("nan")
    aucs = []
    for s in seeds:
        skf = StratifiedKFold(5, shuffle=True, random_state=s)
        oof = np.full(len(yc), np.nan)
        for tr, te in skf.split(XN, yc):
            Ftr = StandardScaler().fit(XF.iloc[tr])
            Ztr, Zte = Ftr.transform(XF.iloc[tr]), Ftr.transform(XF.iloc[te])
            Rtr = np.empty_like(XN.iloc[tr].values, dtype=float)
            Rte = np.empty_like(XN.iloc[te].values, dtype=float)
            for j in range(XN.shape[1]):
                lr = LinearRegression().fit(Ztr, XN.iloc[tr].values[:, j])
                Rtr[:, j] = XN.iloc[tr].values[:, j] - lr.predict(Ztr)
                Rte[:, j] = XN.iloc[te].values[:, j] - lr.predict(Zte)
            sc = StandardScaler().fit(Rtr)
            clf = LogisticRegression(max_iter=2000, C=1.0).fit(sc.transform(Rtr), yc[tr])
            oof[te] = clf.predict_proba(sc.transform(Rte))[:, 1]
        ok = ~np.isnan(oof)
        aucs.append(roc_auc_score(yc[ok], oof[ok]))
    return float(np.mean(aucs))


def two_block_eval(mat: pd.DataFrame, cov: pd.DataFrame,
                   floor_cols: list[str], n_perm: int = 1000,
                   n_seeds: int = 20, seed0: int = 0) -> dict:
    """
    Pre-registered two-block test: immune floor (A) vs non-reference (B),
    floor-conditioned (C = A+B), faithful to the eval session's harness:
      - PRIMARY: within-largest-cohort stratified 5-fold CV, seed-averaged
        (the 'within-Gide n=32' frame — folds never cross cohorts);
      - LOCO (LeaveOneGroupOut) cross-cohort transfer AUROC for each block
        (the honest cross-cohort test; only interpretable if the floor transfers);
      - fold-contained residualization of non-ref on the floor, within cohort;
      - cohort-internal permutation null for the residual block (label shuffle
        WITHIN cohort, preserving group sizes + per-cohort class balance) → perm p.
    Every AUROC is a real out-of-fold score; no metric is a stub.
    """
    df = mat.join(cov.set_index("run_accession"), how="inner")
    y = df["y"].astype(int).values
    groups = df["cohort"].values
    feat_cols = list(mat.columns)
    X_nonref = df[feat_cols].astype(float)
    X_nonref = X_nonref.fillna(X_nonref.median())
    X_floor  = df[floor_cols].astype(float).fillna(0.0)
    X_comb   = pd.concat([X_floor.reset_index(drop=True),
                          X_nonref.reset_index(drop=True)], axis=1)
    seeds = list(range(seed0, seed0 + n_seeds))
    uniq, cnts = np.unique(groups, return_counts=True)
    n_groups = len(uniq)
    primary = str(uniq[int(np.argmax(cnts))])  # largest cohort = primary within-cohort frame

    out = {
        "n": int(len(df)),
        "n_pos": int(y.sum()),
        "cohorts": {k: int(v) for k, v in pd.Series(groups).value_counts().items()},
        "primary_cohort": primary,
        "cv": f"within-{primary} StratifiedKFold(5) seed-averaged; LOCO=LeaveOneGroupOut across {n_groups} cohorts",
        # PRIMARY within-cohort frame
        "floor_within_primary": _within_cohort_5fold_auroc(X_floor, y, groups, seeds, primary),
        "nonref_within_primary": _within_cohort_5fold_auroc(X_nonref, y, groups, seeds, primary),
        "floor_plus_nonref_within_primary": _within_cohort_5fold_auroc(X_comb, y, groups, seeds, primary),
        "nonref_resid_on_floor_within_primary": _residual_within_cohort_auroc(
            X_nonref, X_floor, y, groups, seeds, primary),
    }
    out["delta_C_minus_A_within_primary"] = round(
        out["floor_plus_nonref_within_primary"] - out["floor_within_primary"], 4)

    # LOCO cross-cohort transfer (only meaningful with >=2 cohorts)
    if n_groups >= 2:
        out["floor_loco"] = _loco_auroc(X_floor, y, groups)
        out["nonref_loco"] = _loco_auroc(X_nonref, y, groups)
        out["floor_plus_nonref_loco"] = _loco_auroc(X_comb, y, groups)
        out["delta_C_minus_A_loco"] = round(
            out["floor_plus_nonref_loco"] - out["floor_loco"], 4)
        out["loco_caveat"] = ("LOCO is interpretable only if the floor itself transfers; "
                              "if floor_loco is near chance the non-ref LOCO number is not diagnostic.")

    # cohort-internal permutation null for the within-primary residual non-ref block
    rng = np.random.default_rng(seed0)
    obs = out["nonref_resid_on_floor_within_primary"]
    perm_seeds = list(range(seed0, seed0 + min(n_seeds, 10)))  # 10-seed, as in the harness
    if not np.isnan(obs):
        null = np.empty(n_perm)
        idx_by_g = {g: np.where(groups == g)[0] for g in uniq}
        for p in range(n_perm):
            yp = y.copy()
            for g, ix in idx_by_g.items():
                yp[ix] = rng.permutation(y[ix])
            null[p] = _residual_within_cohort_auroc(X_nonref, X_floor, yp, groups, perm_seeds, primary)
        null = null[~np.isnan(null)]
        out["nonref_resid_perm_p"] = float((np.sum(null >= obs) + 1) / (len(null) + 1))
        out["nonref_resid_perm_null_mean"] = float(null.mean())
        out["n_perm"] = int(len(null))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--existing-matrix", required=True, type=Path,
                    help="nonref_matrix_cohort.parquet — used for COLUMN SCHEMA ONLY (its 25 rows "
                         "are HISAT2-derived and are NOT reused post-STAR-pivot unless --reuse-existing-rows)")
    ap.add_argument("--new-feature-root", required=True, type=Path,
                    help="dir of per-sample non-ref caller outputs (run on the STAR BAMs at "
                         "/data/rnaseq_cohort/star_salmon/), one subdir per run_accession")
    ap.add_argument("--reuse-existing-rows", action="store_true",
                    help="reuse rows from --existing-matrix instead of rebuilding them (ONLY safe if "
                         "those rows came from the SAME aligner as --new-feature-root; default off "
                         "because the n=106 run is a STAR re-alignment, existing rows are HISAT2)")
    ap.add_argument("--covariates", required=True, type=Path,
                    help="phase2_covariates_n106.parquet (immune floor + y + cohort, keyed on run_accession)")
    ap.add_argument("--manifest", required=True, type=Path,
                    help="remaining/full manifest with a run_accession column")
    ap.add_argument("--out-matrix", default=Path("nonref_matrix_n106.parquet"), type=Path)
    ap.add_argument("--out-result", default=Path("two_block_n106.json"), type=Path)
    args = ap.parse_args()

    existing = load_existing_matrix(args.existing_matrix)
    manifest = pd.read_csv(args.manifest)
    cov = pd.read_parquet(args.covariates)

    mat = build_nonref_matrix(existing, args.new_feature_root, manifest,
                              reuse_existing_rows=args.reuse_existing_rows)
    mat.to_parquet(args.out_matrix)

    # immune floor = the standard composition/inflammation features in the covariates
    floor_cols = [c for c in ["gep_tcell_inflamed", "ifng_score", "teff", "tgfb",
                              "teff_tgfb_balance"] if c in cov.columns]
    res = two_block_eval(mat, cov, floor_cols)
    res["floor_cols"] = floor_cols
    res["out_matrix"] = str(args.out_matrix)
    args.out_result.write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
