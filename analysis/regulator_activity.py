#!/usr/bin/env python
"""
regulator_activity.py — expression-derived RNA-regulator activity features.

Wraps the rbp-activity-scorer skill's REGULATOR_SETS + score_regulator_activity
into a pilot_ingest-compatible builder: read the gene-TPM matrix, score
splicing-factor / broad-RBP / ADAR-editing activity per sample, return a
[run_accession, cohort, <3 activity cols>] frame keyed on the contract.

These are NAMED, computable-now proxies for the differentiated RNA state
(upstream regulator expression, not de-novo antigen burden). On the first
3-cohort pilot (n=52) none beat the TIDE floor — recorded as an honest interim
null, not the hypothesis test.

The skill's kernel plugin injects REGULATOR_SETS + score_regulator_activity into
an interactive kernel; for a standalone module run we re-derive the gene sets
from that plugin if present, else fall back to importing it.
"""
import os
import pandas as pd

ACTIVITY_COLUMNS = ["SPLICING_FACTOR", "RBP_BROAD", "ADAR_EDITING"]
KEY = ["run_accession", "cohort"]


# ---------------------------------------------------------------------------
# Fallback vendored copy of the rbp-activity-scorer skill's REGULATOR_SETS +
# score_regulator_activity. The LOADED SKILL is the source of truth (primary
# path below reads it from builtins); this fallback lets a standalone
# orchestrator subprocess — which does not get the skill's kernel plugin —
# compute the same features. Keep in sync with the skill if the gene sets change.
# ---------------------------------------------------------------------------
_VENDORED_REGULATOR_SETS = {
    "SPLICING_FACTOR": ["SF3B1", "SF3B2", "SF3B4", "SF1", "U2AF1", "U2AF2", "SRSF1", "SRSF2",
        "SRSF3", "SRSF5", "SRSF6", "SRSF7", "SRSF10", "HNRNPA1", "HNRNPA2B1", "HNRNPC", "HNRNPD",
        "HNRNPK", "HNRNPM", "HNRNPU", "PRPF8", "PRPF19", "SNRNP70", "SNRPB", "RBM10", "RBM39",
        "RBM25", "RBFOX2", "PTBP1", "PTBP2", "TRA2B", "CELF1", "MBNL1", "QKI", "ESRP1", "ESRP2"],
    "RBP_BROAD": ["ELAVL1", "IGF2BP1", "IGF2BP2", "IGF2BP3", "YBX1", "DDX3X", "DHX9", "PABPC1",
        "EIF4E", "LIN28B", "MSI2", "PCBP1", "PCBP2", "KHDRBS1", "FUS", "TARDBP", "MATR3"],
    "ADAR_EDITING": ["ADAR", "ADARB1", "ADARB2"],
}


def _vendored_score(tpm_df, gene_id_map, sets=None, min_genes=3):
    import numpy as np
    if sets is None:
        sets = _VENDORED_REGULATOR_SETS
    logx = np.log2(tpm_df.astype(float) + 1.0)
    z = (logx - logx.mean(0)) / logx.std(0).replace(0, np.nan)
    out = {}
    for name, syms in sets.items():
        ensg = [gene_id_map[s] for s in syms if s in gene_id_map]
        present = [g for g in ensg if g in z.columns and np.isfinite(z[g]).all()]
        if len(present) >= min_genes:
            out[name] = z[present].mean(axis=1)
    return pd.DataFrame(out)


def _load_skill_api():
    """Return (REGULATOR_SETS, score_regulator_activity). Prefer the LOADED
    skill (builtins injected by its kernel plugin); fall back to the vendored
    copy for standalone / subprocess runs where the plugin is absent."""
    import builtins
    rs = getattr(builtins, "REGULATOR_SETS", None)
    fn = getattr(builtins, "score_regulator_activity", None)
    if rs is not None and fn is not None:
        return rs, fn
    g = globals()
    if "REGULATOR_SETS" in g and "score_regulator_activity" in g:
        return g["REGULATOR_SETS"], g["score_regulator_activity"]
    return _VENDORED_REGULATOR_SETS, _vendored_score


def regulator_sets():
    """Return the active REGULATOR_SETS (loaded skill if present, else vendored)."""
    rs, _ = _load_skill_api()
    return rs


def build_regulator_activity(tpm_samples_by_gene, gene_symbol_index=None,
                             sample_meta=None):
    """Score regulator activity.

    tpm_samples_by_gene : DataFrame rows=samples (run_accession index), cols=Ensembl
        gene ids (unversioned), values=linear TPM.
    gene_symbol_index : optional dict symbol->ensembl for the regulator genes.
        If None, the caller must have already restricted columns to the regulator
        genes' Ensembl ids.
    sample_meta : optional DataFrame with run_accession + cohort to attach.
    """
    REGULATOR_SETS, score = _load_skill_api()
    if gene_symbol_index is None:
        raise ValueError("gene_symbol_index (symbol->ensembl) required")
    S = score(tpm_samples_by_gene, gene_symbol_index)
    S = S.reset_index().rename(columns={"index": "run_accession"})
    if sample_meta is not None:
        S = S.merge(sample_meta[[c for c in KEY if c in sample_meta.columns]].drop_duplicates(),
                    on="run_accession", how="left")
    return S


def gene_map_from_matrix(expr_with_gene_name, symbols):
    """Build symbol->unversioned-Ensembl map from a matrix that carries a
    'gene_name' column and Ensembl (possibly versioned) row index."""
    tmp = expr_with_gene_name.copy()
    tmp["_ens"] = [str(i).split(".")[0] for i in tmp.index]
    name2ens = tmp.set_index("gene_name")["_ens"].to_dict()
    return {s: name2ens[s] for s in symbols if s in name2ens}
