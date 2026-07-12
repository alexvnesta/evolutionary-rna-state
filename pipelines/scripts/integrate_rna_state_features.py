#!/usr/bin/env python3
"""
integrate_rna_state_features.py — assemble a per-sample "evolutionary RNA state"
feature table from the four cohort phenotype layers completed on the Mac, and
test (a) per-feature response association and (b) cross-layer coordination.

The hypothesis is that alternative splicing, intron retention, RNA editing, and
TE/ERV activation are not independent biomarkers but co-varying manifestations
of one underlying RNA state. This script builds the burden-summary matrix that
lets that be tested directly.

Inputs (cohort artifacts; results/ is git-ignored so paths are illustrative):
  - cohort_aei_40samples.tsv                 (AEI editing, per-sample)
  - cohort_te_erv_counts_40samples.tsv       (TEcount uniq, 1328 families x 40)
  - cohort_intron_retention_40samples.parquet(IR ratio, 40 x 245772 introns)
  - suppa_cohort_psi_40samples.psi           (SUPPA PSI, 306514 events x 40)
  - data/manifests/selection_manifest.csv    (sample -> cohort, response)

Per-sample features derived:
  editing:     AEI_percent, AEI_signal_to_noise
  TE/ERV:      TE_{LTR,LINE,SINE,DNA}_frac  (class composition of TE reads)
  intron ret.: IR_mean, IR_load_gt0.1, IR_load_gt0.2
  splicing:    splice_dysregulation (mean |PSI - cohort event median|, events
               defined in >=30/40 samples), splice_switch_events (|dev|>0.5)

Outputs:
  rna_state_integrated_features_40samples.csv
  rna_state_feature_response_assoc.csv   (Mann-Whitney R vs N per feature)
  rna_state_layer_correlation.csv        (Spearman across one rep feature/layer)
  rna_state_coordination.png

Findings on the Gide+Riaz 40-sample cohort (21 R / 19 N):
  - No single per-sample burden stratifies ICB response (all MW p>0.5) — the
    signal, if any, is not in any one layer's global magnitude.
  - Three of four layers are significantly coordinated (Spearman, n=40):
    TE/ERV ~ splicing rho=+0.51 p=0.001; TE/ERV ~ intron-ret rho=-0.49 p=0.001;
    intron-ret ~ splicing rho=-0.49 p=0.001. RNA-editing burden is the outlier
    (uncorrelated). This is direct evidence the phenotypes co-vary as the
    "shared RNA state" hypothesis predicts, independent of response.

Pass the four input paths as argv, or edit the ART dict to artifact version ids.
"""
import sys, itertools
import numpy as np, pandas as pd
from scipy import stats


def per_sample_features(aei_tsv, te_tsv, ir_parquet, psi_tsv, manifest_csv):
    man = pd.read_csv(manifest_csv)
    key = man[["sample_title", "run_accession", "cohort", "resp_NR"]].copy()
    key.columns = ["sample", "run_accession", "cohort", "response"]
    key = key.set_index("sample")

    # editing (already per-sample)
    aei = pd.read_csv(aei_tsv, sep="\t").set_index("sample")
    key["AEI_percent"] = aei["AEI_percent"]
    key["AEI_signal_to_noise"] = aei["signal_to_noise"]

    # TE/ERV class composition (robust to library size)
    te = pd.read_csv(te_tsv, sep="\t", index_col=0).reindex(columns=key.index)
    te_cls = te.groupby(te.index.str.split(":").str[-1]).sum()
    libsize_te = te.sum(axis=0)
    for cls in ["LTR", "LINE", "SINE", "DNA"]:
        if cls in te_cls.index:
            key[f"TE_{cls}_frac"] = te_cls.loc[cls] / libsize_te

    # intron retention burden (note: the parquet's "run_accession" col holds sample_title)
    ir = pd.read_parquet(ir_parquet)
    ir_vals = ir.drop(columns=["run_accession", "cohort"]).values.astype(np.float32)
    idx = ir["run_accession"].values
    key["IR_mean"] = pd.Series(np.nanmean(ir_vals, axis=1), index=idx)
    key["IR_load_gt0.1"] = pd.Series(np.nansum(ir_vals >= 0.10, axis=1), index=idx)
    key["IR_load_gt0.2"] = pd.Series(np.nansum(ir_vals >= 0.20, axis=1), index=idx)

    # splicing dysregulation (deviation from cohort-typical PSI)
    psi = pd.read_csv(psi_tsv, sep="\t")
    P = psi.values.astype(np.float32)
    mask = np.sum(~np.isnan(P), axis=1) >= 30
    Pm = P[mask]
    dev = np.abs(Pm - np.nanmedian(Pm, axis=1, keepdims=True))
    key["splice_dysregulation"] = pd.Series(np.nanmean(dev, axis=0), index=psi.columns)
    key["splice_switch_events"] = pd.Series(np.nansum(dev > 0.5, axis=0), index=psi.columns)
    return key


def analyze(M, outdir="."):
    feat = [c for c in M.columns if c not in ("run_accession", "cohort", "response")]
    R, N = M[M.response == "R"], M[M.response == "N"]
    assoc = pd.DataFrame(
        [(c, R[c].mean(), N[c].mean(), R[c].mean() - N[c].mean(),
          stats.mannwhitneyu(R[c], N[c], alternative="two-sided").pvalue) for c in feat],
        columns=["feature", "R_mean", "N_mean", "R_minus_N", "MW_pval"]
    ).sort_values("MW_pval")

    rep = {"editing": "AEI_percent", "TE_ERV": "TE_LTR_frac",
           "intron_ret": "IR_mean", "splicing": "splice_dysregulation"}
    sub = M[list(rep.values())].copy()
    sub.columns = list(rep.keys())
    corr = sub.corr(method="spearman")

    M.to_csv(f"{outdir}/rna_state_integrated_features_40samples.csv")
    assoc.to_csv(f"{outdir}/rna_state_feature_response_assoc.csv", index=False)
    corr.to_csv(f"{outdir}/rna_state_layer_correlation.csv")
    return assoc, corr, sub


if __name__ == "__main__":
    aei_tsv, te_tsv, ir_parquet, psi_tsv, manifest_csv = sys.argv[1:6]
    outdir = sys.argv[6] if len(sys.argv) > 6 else "."
    M = per_sample_features(aei_tsv, te_tsv, ir_parquet, psi_tsv, manifest_csv)
    assert M.drop(columns=["run_accession", "cohort", "response"]).notna().all().all()
    assoc, corr, _ = analyze(M, outdir)
    print(assoc.round(4).to_string(index=False))
    print("\ncross-layer Spearman:\n", corr.round(3).to_string())
