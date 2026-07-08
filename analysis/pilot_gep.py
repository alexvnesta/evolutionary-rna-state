#!/usr/bin/env python
"""
pilot_gep.py — normalize a real pipeline gene-quantification matrix into the
genes(symbol) x samples orientation that gep_scores.score_all expects.

The pipeline emits ``quant_gene_tpm.parquet`` as samples x genes with Ensembl
gene-id columns (ENSG...) plus id columns run_accession/cohort. gep_scores wants
genes x samples indexed by HGNC symbol. This module bridges the two: detects
orientation, maps the signature-gene Ensembl ids to symbols, and returns
(matrix, sample_meta). Only the ~60 signature/housekeeping genes are mapped
(one small metadata fetch), not the whole transcriptome.
"""
import os, json
import pandas as pd

from analysis.baseline import gep_scores as _gs

KEY = ["run_accession", "cohort"]

# Curated symbol->Ensembl for the union of GEP / IFN-gamma / housekeeping /
# Teff / TGF-beta signature genes. Falls back to a live mygene.info lookup for
# any symbol not present here, so the map stays correct if signatures change.
_STATIC_SYM2ENS = {
    "CTGF": "ENSG00000118523",  # CCN2 — old symbol not in symbol scope
}


def _signature_symbols():
    need = set()
    for k in ("GEP_TCELL_INFLAMED_GENES", "IFNG_GENES", "AYERS_HOUSEKEEPING_GENES",
              "TEFF_GENES", "TGFB_GENES"):
        need |= set(getattr(_gs, k, []))
    return need


def _mygene_symbol_to_ensembl(symbols, contact_email=None):
    import urllib.request, urllib.parse
    body = urllib.parse.urlencode({
        "q": ",".join(sorted(symbols)), "scopes": "symbol",
        "fields": "ensembl.gene", "species": "human",
    }).encode()
    req = urllib.request.Request("https://mygene.info/v3/query", data=body,
                                 headers={"Content-Type": "application/x-www-form-urlencoded"})
    with urllib.request.urlopen(req, timeout=60) as r:
        res = json.load(r)
    out = {}
    for row in res:
        q, ens = row.get("query"), row.get("ensembl")
        if not ens:
            continue
        g = ens["gene"] if isinstance(ens, dict) else (ens[0]["gene"] if isinstance(ens, list) and ens else None)
        if g and q not in out:
            out[q] = g
    return out


def signature_sym2ens(cache_path=None):
    """symbol->Ensembl map covering the signature genes (cached if requested)."""
    if cache_path and os.path.exists(cache_path):
        return json.load(open(cache_path))
    need = _signature_symbols()
    m = _mygene_symbol_to_ensembl(need)
    m.update(_STATIC_SYM2ENS)
    if cache_path:
        json.dump(m, open(cache_path, "w"))
    return m


def to_symbol_gene_matrix(raw, cache_path=None):
    """Return (matrix, sample_meta): matrix is genes(symbol) x samples for the
    signature genes present; sample_meta is a run_accession/cohort frame indexed
    by run_accession. Accepts the real pilot layout (samples x ENSG columns) or an
    already-symbol-indexed genes x samples frame.
    """
    # Already genes x samples with symbols? (index holds signature symbols)
    if raw.index.astype(str).isin(_signature_symbols()).any():
        mat = raw
        sm = pd.DataFrame({"run_accession": list(mat.columns)})
        sm["cohort"] = None
        return mat, sm.set_index("run_accession", drop=False)

    # genes x samples with a 'gene_name' symbol column + Ensembl (versioned) index
    # and run-accession sample columns (the 3-cohort n=52 pilot layout).
    if "gene_name" in raw.columns:
        sample_cols = [c for c in raw.columns
                       if c != "gene_name" and not str(c).startswith("ENSG")]
        mat = raw.set_index("gene_name")[sample_cols]
        mat = mat.groupby(level=0).max()  # collapse duplicate symbols
        sm = pd.DataFrame({"run_accession": sample_cols})
        sm["cohort"] = None
        return mat, sm.set_index("run_accession", drop=False)

    # Real pilot: samples x genes, ENSG columns + id columns
    id_cols = [c for c in KEY if c in raw.columns]
    gene_cols = [c for c in raw.columns if str(c).startswith("ENSG")]
    sym2ens = signature_sym2ens(cache_path)
    ens2sym = {v: k for k, v in sym2ens.items()}
    keep = [c for c in gene_cols if c in ens2sym]
    sub = raw.set_index("run_accession")
    mat = sub[keep].T
    mat.index = [ens2sym[e] for e in mat.index]
    sm = raw[id_cols].drop_duplicates().set_index("run_accession", drop=False)
    return mat, sm
