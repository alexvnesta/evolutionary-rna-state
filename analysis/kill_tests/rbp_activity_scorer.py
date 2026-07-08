"""
rbp_activity_scorer.py — regulator-activity feature scorer for IO-response modeling.

Computes per-sample activity scores for curated RNA-regulator gene sets
(splicing factors, broad RBPs, ADAR editing enzymes) from a gene-expression
matrix. Built as the upstream-representation alternative to summed per-source
antigen-burden features, which showed no coordination, no marginal prediction,
and no interaction/nonlinear signal over TMB (see io_rna_antigen_analysis.md).

Usage:
    from rbp_activity_scorer import score_regulator_activity, REGULATOR_SETS
    S = score_regulator_activity(tpm_df, gene_id_map)   # samples x sets
    # tpm_df: rows=samples, cols=Ensembl gene IDs (linear TPM)
    # gene_id_map: {HGNC_symbol: ENSG...}  (from BioMart batch_translate)
"""
import numpy as np
import pandas as pd

# Curated canonical regulator gene sets (HGNC symbols). All 56 resolve in Ensembl.
REGULATOR_SETS = {
    "SPLICING_FACTOR": ["SF3B1","SF3B2","SF3B4","SF1","U2AF1","U2AF2","SRSF1","SRSF2",
        "SRSF3","SRSF5","SRSF6","SRSF7","SRSF10","HNRNPA1","HNRNPA2B1","HNRNPC","HNRNPD",
        "HNRNPK","HNRNPM","HNRNPU","PRPF8","PRPF19","SNRNP70","SNRPB","RBM10","RBM39",
        "RBM25","RBFOX2","PTBP1","PTBP2","TRA2B","CELF1","MBNL1","QKI","ESRP1","ESRP2"],
    "RBP_BROAD": ["ELAVL1","IGF2BP1","IGF2BP2","IGF2BP3","YBX1","DDX3X","DHX9","PABPC1",
        "EIF4E","LIN28B","MSI2","PCBP1","PCBP2","KHDRBS1","FUS","TARDBP","MATR3"],
    "ADAR_EDITING": ["ADAR","ADARB1","ADARB2"],
}

def score_regulator_activity(tpm_df, gene_id_map, sets=None, min_genes=3):
    """Mean z-scored log2-TPM activity per regulator set.

    tpm_df       : DataFrame, rows=samples, cols=Ensembl gene IDs, values=linear TPM.
    gene_id_map  : dict HGNC_symbol -> Ensembl gene ID.
    sets         : dict set_name -> [HGNC symbols]; defaults to REGULATOR_SETS.
    min_genes    : minimum member genes present+variable to emit a score.
    Returns      : DataFrame rows=samples, cols=set names (activity z).
    """
    sets = sets or REGULATOR_SETS
    logx = np.log2(tpm_df.astype(float) + 1.0)
    z = (logx - logx.mean(0)) / logx.std(0).replace(0, np.nan)
    out = {}
    for name, syms in sets.items():
        ensg = [gene_id_map[s] for s in syms if s in gene_id_map]
        present = [g for g in ensg if g in z.columns and np.isfinite(z[g]).all()]
        if len(present) >= min_genes:
            out[name] = z[present].mean(axis=1)
    return pd.DataFrame(out)
