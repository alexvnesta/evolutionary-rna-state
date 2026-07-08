"""
reframing_tests.py — the three tests that adjudicate the RNA-state hypothesis reframing.

After four threads (see SYNTHESIS_hypothesis_reconsideration.md), the antigen-QUANTITY
link is a well-supported negative and the surviving lead is a coordinated, proliferation-
independent RNA-regulator ACTIVITY state that appears to mark IMMUNE-COLD tumors rather
than drive antigenicity. These three tests operationalize that reframing so it is checked
automatically as cohorts (esp. liu2019, n=122 w/ TMB) land and n grows from 40 to ~150.

Inputs: a per-sample gene TPM matrix (rows=run_accession, cols=ENSG...), the regulator
ENSG map (from move2_autorun.regulator_ensembl_map), and a joined response table.

PREDICTIONS being tested (falsifiable):
  P1 coordination holds : activity PC1 var-explained > permutation null at n~150 (p<0.05)
  P2 immune-cold coupling: partial corr(activity PC1, IFN sig | generic expr) stays negative & sig
  P3 no incremental value: CV-AUROC(IFN + activity) - CV-AUROC(IFN) ~ 0
If P1 fails -> even the narrow RNA-state claim dies. If P3 shows a positive delta ->
resurrects a direct role worth chasing.
"""
from __future__ import annotations
import json, urllib.request
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_predict, StratifiedKFold
from sklearn.metrics import roc_auc_score

# canonical marker sets (symbols); resolved to ENSG at runtime and cached by caller
PROLIF = ["MKI67","PCNA","TOP2A","CCNB1","CCNB2","CDK1","AURKA","BUB1","BUB1B","CCNA2",
          "CDC20","MCM2","MCM6","FEN1","RRM2","TYMS","UBE2C","BIRC5","CENPA","KIF11"]
IFN_TCELL = ["IFNG","STAT1","CCR5","CXCL9","CXCL10","CXCL11","IDO1","PRF1","GZMA","MHC2TA",
             "HLA-DRA","CD3D","CD3E","CD2","IL2RG","NKG7","CCL5","LAG3","TAGAP","CD8A","CD27","CXCR6"]


def _ens(symbols):
    req = urllib.request.Request("https://rest.ensembl.org/lookup/symbol/homo_sapiens",
        data=json.dumps({"symbols": symbols}).encode(),
        headers={"Content-Type": "application/json", "Accept": "application/json"})
    r = json.load(urllib.request.urlopen(req, timeout=90))
    return {s: i["id"] for s, i in r.items()
            if isinstance(i, dict) and str(i.get("id", "")).startswith("ENSG")}


def _signature(logX, symbols):
    emap = _ens(symbols); cols = [e for e in emap.values() if e in logX.columns]
    z = (logX[cols] - logX[cols].mean()) / logX[cols].std()
    return z.mean(axis=1).values, len(cols)


def _partial_corr(a, b, c):
    a = a - np.polyval(np.polyfit(c, a, 1), c)
    b = b - np.polyval(np.polyfit(c, b, 1), c)
    return stats.pearsonr(a, b)


def run(tpm, reg_ensg, response_df, run_col="run_accession", resp_col="RESPONDER",
        n_perm=5000, seed=0):
    """tpm: rows=run_accession, cols=ENSG (raw TPM). reg_ensg: list of regulator ENSG ids
    present in tpm. response_df: has [run_accession, RESPONDER]. Returns a dict report."""
    logX = np.log2(tpm + 1)
    reg = [e for e in reg_ensg if e in logX.columns]
    z_reg = ((logX[reg] - logX[reg].mean()) / logX[reg].std()).fillna(0)

    # activity PC1 + variance explained
    pca = PCA(n_components=min(5, z_reg.shape[1])).fit(z_reg)
    act_pc1 = pca.transform(z_reg)[:, 0]
    ve = float(pca.explained_variance_ratio_[0])

    # P1: permutation null for PC1 variance-explained (shuffle genes within sample)
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    A = z_reg.values
    for i in range(n_perm):
        Ap = np.array([rng.permutation(row) for row in A])
        null[i] = PCA(n_components=1).fit(Ap).explained_variance_ratio_[0]
    p1 = float((null >= ve).mean())

    generic = logX.mean(axis=1).values
    prolif, n_pro = _signature(logX, PROLIF)
    ifn, n_ifn = _signature(logX, IFN_TCELL)

    r_pro = stats.pearsonr(act_pc1, prolif)
    r_ifn = _partial_corr(act_pc1, ifn, generic)

    # join response
    sc = pd.DataFrame({run_col: tpm.index.values, "act_pc1": act_pc1,
                       "ifn": ifn, "prolif": prolif})
    sc = sc.merge(response_df[[run_col, resp_col]], on=run_col, how="inner").dropna()
    sc["y"] = sc[resp_col].astype(int)

    out = {"n": int(len(tpm)), "n_regulators": len(reg), "n_prolif": n_pro, "n_ifn": n_ifn,
           "P1_activity_pc1_var": round(ve, 4), "P1_null_mean": round(float(null.mean()), 4),
           "P1_p_value": round(p1, 4),
           "prolif_independence_r": round(float(r_pro[0]), 4),
           "prolif_independence_p": round(float(r_pro[1]), 4),
           "P2_immune_cold_partial_r": round(float(r_ifn[0]), 4),
           "P2_immune_cold_p": round(float(r_ifn[1]), 4),
           "n_with_response": int(len(sc))}

    if len(sc) >= 20 and sc.y.nunique() == 2:
        cv = StratifiedKFold(5, shuffle=True, random_state=seed)
        def auroc(cols):
            p = cross_val_predict(LogisticRegression(max_iter=1000), sc[cols].values,
                                  sc.y.values, cv=cv, method="predict_proba")[:, 1]
            return float(roc_auc_score(sc.y.values, p))
        a_ifn = auroc(["ifn"]); a_both = auroc(["ifn", "act_pc1"]); a_act = auroc(["act_pc1"])
        for feat in ["ifn", "act_pc1", "prolif"]:
            rho, p = stats.spearmanr(sc[feat], sc.y)
            out[f"rho_{feat}_response"] = round(float(rho), 4)
            out[f"p_{feat}_response"] = round(float(p), 4)
        out.update({"P3_auroc_ifn": round(a_ifn, 4), "P3_auroc_ifn_plus_activity": round(a_both, 4),
                    "P3_incremental_delta": round(a_both - a_ifn, 4),
                    "P3_auroc_activity_alone": round(a_act, 4)})
    else:
        out["P3"] = "insufficient (need n>=20, both classes)"
    return out, sc
