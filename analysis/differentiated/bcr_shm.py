"""
analysis/differentiated/bcr_shm.py
==================================

DIFFERENTIATED / RNA-native axis — B-cell receptor repertoire, somatic
hypermutation (SHM), isotype class-switch, and clonality from bulk RNA-seq.

This is the repertoire-level, RNA-observable reading of B-cell affinity
maturation: the same biology as the population-level TLS/B-cell expression
score (analysis/baseline/tls_bcell_scores.py), measured at the level of
reconstructed immunoglobulin sequences. It is built on TRUST4
(Song et al., Nat Methods 2021, doi:10.1038/s41592-021-01142-2), which
assembles BCR/TCR contigs directly from bulk RNA-seq reads.

Inputs (per sample), produced by TRUST4:
    <prefix>_cdr3.out    -- per-consensus CDR3 + V/D/J/C genes +
                            CDR3_germline_similarity + read_fragment_count
    <prefix>_report.tsv  -- clonotype report (read_count, frequency, CDR3, VDJC)
    <prefix>_airr.tsv    -- AIRR rearrangements (v_identity when present) —
                            preferred SHM source when available

Per-sample features (keyed on run_accession, cohort):
    bcr_shm_rate          -- mean somatic-hypermutation rate over BCR (IGH)
                             clonotypes = 1 - V-gene germline identity. The
                             direct molecular read-out of affinity maturation.
    bcr_igg_fraction      -- class-switched IgG proportion (IGHG / all-isotype)
    bcr_switched_fraction -- all class-switched (IGHG+IGHA+IGHE) / total
    bcr_clonality         -- 1 - normalised Shannon entropy of IGH clonotype
                             frequencies (0 = polyclonal, 1 = monoclonal)
    bcr_n_clonotypes      -- number of distinct IGH clonotypes recovered
    bcr_n_reads           -- total IGH-assigned reads (assembly-depth / QC)

DESIGN / HONESTY NOTES
----------------------
* SHM source priority: AIRR v_identity if present and non-null, else
  (1 - CDR3_germline_similarity) from _cdr3.out. Which source was used is
  recorded per sample in the returned `_shm_source` column so it is auditable.
* SHM is averaged over IGH clonotypes. For the cdr3_germline_similarity source
  the mean is read_fragment_count-weighted by default (a read-abundant,
  well-supported contig is a more reliable SHM estimate than a singleton);
  weight=False switches that source to an unweighted mean. The AIRR
  v_identity source is ALWAYS an unweighted mean over rearrangements (the AIRR
  tsv is not parsed for per-rearrangement read support), so the `weight` flag
  does not affect it.
* DEPTH SENSITIVITY (the honest caveat): BCR contig yield and therefore the
  stability of every feature here scales with sequencing depth and with how
  B-cell-rich the sample is. A sample with too few IGH reads gives an
  unreliable SHM rate / clonality. `min_clonotypes` (default 3) sets the
  contigs below which features are returned as NaN rather than a noisy point
  estimate — NEVER imputed. bcr_n_clonotypes / bcr_n_reads are reported so the
  reliability of each sample is visible downstream.
* This module PARSES TRUST4 output; it does not run TRUST4. Running is done by
  tools/run_trust4_poc.sh (pilot) / the wired pipeline stage on the raw reads.
"""
from __future__ import annotations

import glob
import os
from typing import Iterable

import numpy as np
import pandas as pd

# Isotype prefixes (heavy-chain constant regions).
_IGH_C = ("IGHM", "IGHD", "IGHG", "IGHA", "IGHE")
_SWITCHED = ("IGHG", "IGHA", "IGHE")  # class-switched (post-germinal-center)

FEATURE_COLUMNS = (
    "bcr_shm_rate", "bcr_igg_fraction", "bcr_switched_fraction",
    "bcr_clonality", "bcr_n_clonotypes", "bcr_n_reads",
)

_CDR3_COLS = [
    "consensus_id", "index_within_consensus", "V_gene", "D_gene", "J_gene",
    "C_gene", "CDR1", "CDR2", "CDR3", "CDR3_score", "read_fragment_count",
    "CDR3_germline_similarity", "complete_vdj_assembly",
]


def _isotype(c_gene: str) -> str | None:
    """Map a TRUST4 C_gene string to a heavy-chain isotype prefix, else None."""
    if not isinstance(c_gene, str) or c_gene in ("", "*", ".", "None"):
        return None
    for iso in _IGH_C:
        if c_gene.startswith(iso):
            return iso
    return None


def _read_cdr3(path: str) -> pd.DataFrame:
    """Parse a TRUST4 *_cdr3.out file (no header) into a typed frame."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame(columns=_CDR3_COLS)
    df = pd.read_csv(path, sep="\t", header=None, names=_CDR3_COLS, dtype=str)
    df["read_fragment_count"] = pd.to_numeric(df["read_fragment_count"], errors="coerce").fillna(0.0)
    df["CDR3_germline_similarity"] = pd.to_numeric(df["CDR3_germline_similarity"], errors="coerce")
    return df


def _airr_v_identity(airr_path: str) -> pd.Series | None:
    """Return per-row v_identity (fraction, 0-1) from an AIRR tsv if usable."""
    if not os.path.exists(airr_path) or os.path.getsize(airr_path) == 0:
        return None
    try:
        a = pd.read_csv(airr_path, sep="\t", dtype=str)
    except Exception:
        return None
    if "v_identity" not in a.columns:
        return None
    v = pd.to_numeric(a["v_identity"], errors="coerce").dropna()
    if v.empty:
        return None
    # AIRR v_identity may be a fraction (0-1) or a percentage (0-100); normalise.
    if v.max() > 1.5:
        v = v / 100.0
    return v


def sample_bcr_features(
    prefix_dir: str,
    prefix: str,
    *,
    weight: bool = True,
    min_clonotypes: int = 3,
) -> dict:
    """Compute BCR/SHM features for ONE sample from its TRUST4 output dir.

    prefix_dir/prefix_cdr3.out (+ optional _airr.tsv) must exist. Returns a dict
    of FEATURE_COLUMNS (+ _shm_source, _shm_n). Features below min_clonotypes
    IGH contigs are NaN (not imputed).
    """
    cdr3 = _read_cdr3(os.path.join(prefix_dir, f"{prefix}_cdr3.out"))
    # restrict to heavy chain (IGH) — the SHM/isotype-bearing chain
    igh = cdr3[cdr3["C_gene"].map(_isotype).notna() | cdr3["V_gene"].str.startswith("IGHV", na=False)].copy()
    igh["isotype"] = igh["C_gene"].map(_isotype)
    n_clono = int(len(igh))
    n_reads = float(igh["read_fragment_count"].sum())

    out = {c: np.nan for c in FEATURE_COLUMNS}
    out["bcr_n_clonotypes"] = float(n_clono)
    out["bcr_n_reads"] = n_reads
    out["_shm_source"] = "none"
    out["_shm_n"] = 0.0
    if n_clono < min_clonotypes:
        return out  # too shallow — leave SHM/isotype/clonality NaN

    # ---- SHM rate ----
    airr = _airr_v_identity(os.path.join(prefix_dir, f"{prefix}_airr.tsv"))
    if airr is not None and len(airr) >= min_clonotypes:
        shm_vals = (1.0 - airr).clip(lower=0.0)
        out["bcr_shm_rate"] = float(shm_vals.mean())
        out["_shm_source"] = "airr_v_identity"
        out["_shm_n"] = float(len(shm_vals))
    else:
        sim = igh["CDR3_germline_similarity"].dropna().clip(upper=1.0, lower=0.0)
        if len(sim) >= min_clonotypes:
            shm = (1.0 - sim)
            if weight:
                w = igh.loc[sim.index, "read_fragment_count"].clip(lower=1.0)
                out["bcr_shm_rate"] = float(np.average(shm, weights=w))
            else:
                out["bcr_shm_rate"] = float(shm.mean())
            out["_shm_source"] = "cdr3_germline_similarity"
            out["_shm_n"] = float(len(shm))

    # ---- isotype fractions (read-weighted) ----
    iso_reads = igh.dropna(subset=["isotype"]).groupby("isotype")["read_fragment_count"].sum()
    total_iso = float(iso_reads.sum())
    if total_iso > 0:
        igg = float(iso_reads.get("IGHG", 0.0))
        switched = float(sum(iso_reads.get(i, 0.0) for i in _SWITCHED))
        out["bcr_igg_fraction"] = igg / total_iso
        out["bcr_switched_fraction"] = switched / total_iso

    # ---- clonality: 1 - normalised Shannon entropy over IGH clonotype reads ----
    freqs = igh["read_fragment_count"].values
    freqs = freqs[freqs > 0]
    if len(freqs) >= 2:
        p = freqs / freqs.sum()
        H = -np.sum(p * np.log(p))
        Hmax = np.log(len(p))
        out["bcr_clonality"] = float(1.0 - H / Hmax) if Hmax > 0 else np.nan
    elif len(freqs) == 1:
        out["bcr_clonality"] = 1.0
    return out


def build_bcr_features(
    trust4_root: str,
    sample_index: pd.DataFrame,
    *,
    weight: bool = True,
    min_clonotypes: int = 3,
) -> pd.DataFrame:
    """Assemble the per-sample BCR feature table.

    trust4_root : directory containing one sub-dir per run_accession, each with
        <run>_cdr3.out (TRUST4 --od <root>/<run> -o <run>).
    sample_index : frame with 'run_accession' and 'cohort' columns; defines
        output rows and order. Samples with no TRUST4 output get all-NaN
        features (never imputed).
    """
    rows = []
    for _, r in sample_index.iterrows():
        run = str(r["run_accession"]); cohort = r.get("cohort")
        pdir = os.path.join(trust4_root, run)
        feats = sample_bcr_features(pdir, run, weight=weight, min_clonotypes=min_clonotypes)
        rows.append({"run_accession": run, "cohort": cohort, **feats})
    cols = ["run_accession", "cohort", *FEATURE_COLUMNS, "_shm_source", "_shm_n"]
    return pd.DataFrame(rows)[cols]
