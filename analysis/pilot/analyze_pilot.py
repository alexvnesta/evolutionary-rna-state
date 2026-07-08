"""
analysis/pilot/analyze_pilot.py — turn the Salmon pilot quant into (a) a
per-sample de-novo feature matrix in the HANDOFF_CONTRACT format, and (b) the
headline pilot result: does de-novo transcriptomic signal recapitulate the
iAtlas proxy on the SAME samples, and does it carry response information the
proxy misses?

Inputs:  results/pilot_salmon/<run>/quant.sf  (+ quant.genes.sf)
         results/salmon_pilot_manifest.csv    (run -> response/arm)
         data/cbioportal/gide2019_*.clinical.tsv  (proxy neoantigen loads)

Outputs (results/):
  features/quant_gene_tpm.parquet    — per-sample gene TPM (contract format)
  pilot_denovo_features.csv          — interpretable de-novo summaries per sample
  pilot_concordance.json             — de-novo vs proxy + response stats
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
PILOT = REPO / "results" / "pilot_salmon"
FEAT = REPO / "results" / "features"


def load_quant() -> pd.DataFrame:
    """Load per-sample transcript TPM into a genes/tx x samples matrix.

    Returns transcript-level TPM wide frame: index = transcript id, cols = run.
    """
    idx = pd.read_csv(PILOT / "pilot_index.csv")
    ok = idx[idx["status"] == "ok"]["run_accession"].tolist()
    cols = {}
    for run in ok:
        q = PILOT / run / "quant.sf"
        if q.exists():
            df = pd.read_csv(q, sep="\t")
            cols[run] = df.set_index("Name")["TPM"]
    tpm = pd.DataFrame(cols)
    return tpm


def gene_tpm(tpm_tx: pd.DataFrame, tx2gene: Path) -> pd.DataFrame:
    """Collapse transcript TPM to gene TPM (sum over transcripts)."""
    m = pd.read_csv(tx2gene, sep="\t", header=None, names=["tx", "gene"])
    m = m.set_index("tx")["gene"]
    g = tpm_tx.join(m).groupby("gene").sum(numeric_only=True)
    return g


# ---------------------------------------------------------------------------
# Interpretable de-novo features derivable from transcriptome quant alone.
# (The heavier IR / editing / TE matrices come from the pipeline session; this
# pilot proves the read-level machinery + gives features the WES proxy lacks.)
# ---------------------------------------------------------------------------
def denovo_features(tpm_tx: pd.DataFrame) -> pd.DataFrame:
    """Per-sample de-novo transcriptomic summaries computable from Salmon
    transcript TPM alone (no BAM, no annotation beyond GENCODE):

    - tx_effective_number: exp(Shannon entropy of the transcript-TPM
      distribution) = effective number of expressed transcripts, a
      transcriptome-'diffuseness' quantity the WES neoantigen proxy cannot see.
    - n_expressed_tx: count of transcripts with TPM > 1.

    The heavier phenotype matrices (intron retention, RNA editing / AEI,
    TE/ERV locus counts, splicing PSI) come from the pipeline session per
    HANDOFF_CONTRACT.md; this pilot proves the read-level machinery end to end.
    """
    tpm = tpm_tx.copy()
    tpm.index = tpm.index.str.split("|").str[0]  # ENST ids, strip any suffix
    out = pd.DataFrame(index=tpm.columns)
    p = tpm.div(tpm.sum(axis=0), axis=1).replace(0, np.nan)
    shannon = -(p * np.log(p)).sum(axis=0)
    out["tx_effective_number"] = np.exp(shannon)
    out["n_expressed_tx"] = (tpm > 1).sum(axis=0)
    return out


if __name__ == "__main__":
    FEAT.mkdir(parents=True, exist_ok=True)
    tpm_tx = load_quant()
    print("loaded quant: tx x samples =", tpm_tx.shape)
    if tpm_tx.shape[1] == 0:
        raise SystemExit("no quant.sf found yet")
    tx2gene = REPO / "refs" / "tx2gene.tsv"
    g = gene_tpm(tpm_tx, tx2gene)
    # write contract-format gene TPM matrix (rows=run_accession)
    man = pd.read_csv(REPO / "results" / "salmon_pilot_manifest.csv")
    gene_wide = g.T.reset_index().rename(columns={"index": "run_accession"})
    gene_wide.insert(1, "cohort", "gide2019")
    gene_wide.to_parquet(FEAT / "quant_gene_tpm.parquet", index=False)
    print("wrote features/quant_gene_tpm.parquet", gene_wide.shape)
