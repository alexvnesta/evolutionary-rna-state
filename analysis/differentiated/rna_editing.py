"""
analysis/differentiated/rna_editing.py

DIFFERENTIATED BUCKET — A-to-I RNA-editing index + recoding neo-epitopes.

Two named, interpretable per-sample features on top of the pipeline session's
RNA-editing quantification (JACUSA2 per-site A-to-I calls + the per-sample Alu
Editing Index computed by pipelines/rna_editing):

    alu_editing_index          (float, [0,1]) — the ADAR A-to-I editing "dial"
                                     for a sample: A>G mismatches divided by
                                     total adenosine coverage, pooled over ALL
                                     Alu elements genome-wide (Roth, Levanon
                                     et al., Nat Methods 2019). This is the
                                     PRIMARY, most batch/platform-reproducible
                                     editing feature (see BATCH_ROBUSTNESS_NOTE).
    editing_neoantigen_burden  (int)         — MHC-I binder count over the
                                     ALTERED (I=G recoded) 8-11mer peptides
                                     produced by NONSYNONYMOUS editing sites in
                                     CDS (REDIportal recoding catalog), scored
                                     by the SHARED MHCflurry engine so it is
                                     directly comparable to the splice / TE /
                                     IR / fusion / SNV burdens.

Design (per FEATURE_CONTRACT_v2.md and feature_registry.json):
  * Upstream inputs consumed (pipeline session, pipelines/rna_editing):
      - cohort_aei.tsv / per-sample *.aei.tsv   (compute_aei.py output; cols:
        sample, AEI_percent, AG_mismatches, A_coverage, signal_to_noise,
        noise_floor_percent, cov_A.. , n_AtoG..)  -> alu_editing_index
      - per-sample *.editing_sites.tsv           (filter_editing_sites.py output;
        cols: chrom, pos, strand, ref, alt, edit_freq, edited_reads, coverage)
        -> which recoding sites are actually edited in each sample
      - a RECODING-SITE annotation catalog (REDIportal coding/recoding sites,
        built once into reference/): per site the reference CDS codon window +
        edit offset + strand, so the altered (I=G) peptide is generated
        deterministically. See RecodingSite / load_recoding_catalog.
      - the sample's HLA-I alleles (hla_typing.parquet / build_hla_table).
  * Shared engine: analysis.antigen_core.mhc_binding.count_binders — the SAME
    peptide scorer every antigen module uses. We NEVER touch MHCflurry directly.

Why edited (I=G) CDS sites generate neoantigens
-----------------------------------------------
ADAR deaminates adenosine to inosine; the ribosome and the sequencer both read
inosine as guanosine (I=G). A recoding edit therefore changes a codon and, when
nonsynonymous, swaps an amino acid in the translated protein (the textbook
examples: GRIA2 Q607R, CYFIP2 K320E, BLCAP Q5R, FLNA Q2341R, IGFBP7 R95G).
The recoded residue is absent from the germline-encoded proteome, so 8-11mers
straddling it are tumour-associated antigenic sequences invisible to a DNA/WES
neoantigen pipeline. We tile the altered protein window into 8-11mers that
SPAN the recoded residue and hand them to the shared MHC-I engine.

Batch robustness (headline: AEI is the reproducible feature; site calls are not)
--------------------------------------------------------------------------------
See BATCH_ROBUSTNESS_NOTE. In short: the AEI is a single genome-wide RATIO
pooled over millions of Alu adenosines, so it is a coverage-weighted mean
editing level that is stable at modest depth and does not require calling (and
thresholding) individual sites — the property that makes it the batch-tolerant
summary. Per-SITE editing calls, by contrast, are depth-sensitive (detection
power scales with coverage) and therefore more batch-variable; we blunt that for
the recoding burden by (1) scoring only a FIXED catalog of known recoding sites
(no de-novo site discovery, so detection is not depth-driven), (2) requiring an
edit-frequency threshold rather than mere presence, and (3) using MHCflurry's
allele-calibrated percentile rank (batch-invariant) for the binder call. Report
per clinical context / z-scored within cohort; never pool raw counts.

NO real-cohort feature values are fabricated here. This module is a runnable,
unit-validated feature builder; the full-cohort run is deferred to the pilot.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd

# ---------------------------------------------------------------------------
# Shared antigen core — import, do NOT reimplement.
# The module lives at analysis/antigen_core/mhc_binding.py. Support both an
# installed-package import (analysis.antigen_core...) and a direct-path import
# (when the antigen_core dir is on sys.path, as in its own unit tests).
# ---------------------------------------------------------------------------
try:  # package-style import (repo root on sys.path)
    from analysis.antigen_core.mhc_binding import (  # type: ignore
        count_binders, STRONG_BINDER_RANK, WEAK_BINDER_RANK,
    )
    from analysis.antigen_core.hla_typing import ALLELE_COLS  # type: ignore
except Exception:  # pragma: no cover - fallback for path-style import
    import sys
    _CORE = Path(__file__).resolve().parents[1] / "antigen_core"
    if str(_CORE) not in sys.path:
        sys.path.insert(0, str(_CORE))
    from mhc_binding import (  # type: ignore
        count_binders, STRONG_BINDER_RANK, WEAK_BINDER_RANK,
    )
    from hla_typing import ALLELE_COLS  # type: ignore

# ---------------------------------------------------------------------------
# Documented constants (all overridable at call time).
# ---------------------------------------------------------------------------
EDIT_FREQ_THRESHOLD = 0.10   # min per-site editing frequency for a recoding site
#   to count as "edited" in a sample. Matches the pipeline's editing_min_freq
#   default (filter_editing_sites.py --min-edit-freq 0.1). Presence-only calls
#   are depth-driven; a frequency floor is the batch-robust form.
PEP_MIN, PEP_MAX = 8, 11     # MHC-I peptide length window (shared-engine window).
CODON_FLANK = 10             # codons of reference CDS kept EACH SIDE of the
#   recoded codon when building the peptide window (>= PEP_MAX-1 so every
#   8-11mer that spans the recoded residue can be tiled).

STOP = "*"
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# Standard genetic code (DNA codon -> amino acid; '*' = stop). Inosine reads as
# guanosine, so an I=G edit is applied at the DNA level as A -> G before this
# table is consulted.
_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": STOP, "TAG": STOP,
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": STOP, "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

BATCH_ROBUSTNESS_NOTE = (
    "PRIMARY feature alu_editing_index is the Alu Editing Index (Roth/Levanon "
    "2019): AG_mismatches / A_coverage pooled over ALL Alu adenosines "
    "genome-wide. It is batch/platform ROBUST by construction: (1) it is a "
    "single genome-wide RATIO (edited signal / adenosine coverage), i.e. a "
    "coverage-weighted mean editing level, so it is self-normalising for "
    "library size and does not scale with depth; (2) pooling millions of Alu "
    "positions averages out per-site sampling noise, so it is stable at modest "
    "depth where per-site calls are not; (3) it requires NO per-site calling or "
    "thresholding, removing the depth-dependent detection step that makes "
    "site-level editing batch-sensitive; (4) the companion aei_signal_to_noise "
    "(A>G rate / mean non-A>G mismatch rate) is a built-in QC flag for "
    "library/alignment artefacts, and alu_editing_index_cohortz is the "
    "within-cohort robust z for pooling. "
    "SECONDARY feature editing_neoantigen_burden derives from per-SITE calls, "
    "which ARE depth-sensitive (detection power scales with coverage). "
    "Mitigations built in: (a) score only a FIXED catalog of known REDIportal "
    "recoding sites (no de-novo discovery -> detection is not depth-driven), "
    f"(b) require edit frequency >= {EDIT_FREQ_THRESHOLD} rather than mere "
    "presence, (c) the binder call uses MHCflurry percentile rank, calibrated "
    "against a fixed per-allele background and therefore batch/platform "
    "invariant. Report per clinical context / z-scored within cohort; never "
    "pool raw site counts across platforms."
)


# ===========================================================================
# Sequence utilities (stdlib translation — no BioPython dependency needed, and
# fully deterministic for unit tests).
# ===========================================================================
def revcomp(seq: str) -> str:
    """Reverse complement of a DNA string."""
    return seq.translate(_COMPLEMENT)[::-1]


def translate(seq: str, frame: int = 0) -> str:
    """Translate ``seq`` in reading ``frame`` (0/1/2). Unknown codons -> 'X'."""
    prot = []
    for i in range(frame, len(seq) - 2, 3):
        prot.append(_CODON_TABLE.get(seq[i:i + 3].upper(), "X"))
    return "".join(prot)


def apply_ItoG(cds: str, offset: int) -> str:
    """Apply one A-to-I(=G) edit at ``offset`` on a coding-strand CDS window.

    Inosine is read as guanosine, so on the coding strand the edited base
    A becomes G. ``offset`` is 0-based into ``cds`` and MUST point at an 'A'
    (the coding-strand reference base of an A-to-I site); raises ValueError
    otherwise so a mis-oriented / mis-annotated site fails loudly rather than
    silently producing a wrong peptide.
    """
    cds = cds.upper()
    if not (0 <= offset < len(cds)):
        raise ValueError(f"edit offset {offset} out of range for window len {len(cds)}")
    if cds[offset] != "A":
        raise ValueError(
            f"coding-strand base at offset {offset} is {cds[offset]!r}, expected 'A' "
            "(A-to-I edit). Check site strand / window orientation."
        )
    return cds[:offset] + "G" + cds[offset + 1:]


def _kmers_spanning(prot: str, residue_idx: int,
                    kmin: int = PEP_MIN, kmax: int = PEP_MAX) -> list[str]:
    """kmin..kmax substrings of ``prot`` that CONTAIN ``residue_idx``.

    A peptide contains the recoded residue if start <= residue_idx <= end,
    i.e. window [i, i+k) with i <= residue_idx < i+k. Stop-containing windows
    are dropped (a peptide cannot span a stop codon).
    """
    out: list[str] = []
    n = len(prot)
    for k in range(kmin, kmax + 1):
        lo = max(0, residue_idx - k + 1)
        hi = min(residue_idx, n - k)
        for i in range(lo, hi + 1):
            pep = prot[i:i + k]
            if STOP in pep:
                continue
            out.append(pep)
    return out


# ===========================================================================
# Feature 1 — alu_editing_index (the primary, batch-robust editing summary)
# ===========================================================================
def _robust_z(x: pd.Series) -> pd.Series:
    """Median/MAD robust z-score. Returns 0 where MAD==0 (or <2 samples)."""
    x = x.astype(float)
    med = x.median()
    mad = (x - med).abs().median()
    if not mad or pd.isna(mad):
        return pd.Series(0.0, index=x.index)
    return (x - med) / (1.4826 * mad)


def compute_alu_editing_index(
    aei_data: "str | Path | pd.DataFrame",
    cohort: "str | dict[str, str] | None" = None,
    sample_col: str = "sample",
) -> pd.DataFrame:
    """Per-sample Alu Editing Index from the pipeline's AEI table.

    Consumes the ``compute_aei.py`` / ``cohort_aei.tsv`` output. The named
    feature ``alu_editing_index`` is the editing RATIO in [0,1] — recomputed
    directly from the raw numerator/denominator (AG_mismatches / A_coverage)
    when those columns are present (robust to any percent formatting), else
    from AEI_percent / 100.

    Parameters
    ----------
    aei_data : path to cohort_aei.tsv (or a single *.aei.tsv), or a DataFrame
        with at least ``sample`` and either (AG_mismatches, A_coverage) or
        AEI_percent.
    cohort : cohort label. Scalar (single-cohort run) applied to every row; or
        a {sample: cohort} mapping; or None if aei_data already has a ``cohort``
        column.
    sample_col : the sample-id column in the AEI table (default 'sample'); it is
        renamed to the contract key ``run_accession``.

    Returns
    -------
    DataFrame keyed on (run_accession, cohort) with columns:
        alu_editing_index          float [0,1]  — NAMED FEATURE (headline)
        alu_editing_index_percent  float        — same, as percent (0-100)
        aei_ag_mismatches          int          — numerator (QC)
        aei_a_coverage             int          — denominator (QC)
        aei_signal_to_noise        float        — A>G rate / mean control-mismatch rate (QC)
        alu_editing_index_cohortz  float        — within-cohort robust z (pooling form)
    """
    if isinstance(aei_data, (str, Path)):
        df = pd.read_csv(aei_data, sep="\t")
    else:
        df = aei_data.copy()

    if sample_col not in df.columns:
        raise KeyError(f"AEI table has no {sample_col!r} column; got {list(df.columns)}")
    df = df.rename(columns={sample_col: "run_accession"})
    df["run_accession"] = df["run_accession"].astype(str)

    # cohort resolution
    if "cohort" not in df.columns:
        if isinstance(cohort, dict):
            df["cohort"] = df["run_accession"].map(cohort)
        elif cohort is not None:
            df["cohort"] = cohort
        else:
            df["cohort"] = pd.NA

    # editing ratio in [0,1] — prefer raw counts (most robust), else percent/100
    if {"AG_mismatches", "A_coverage"}.issubset(df.columns):
        ag = pd.to_numeric(df["AG_mismatches"], errors="coerce")
        acov = pd.to_numeric(df["A_coverage"], errors="coerce")
        idx = (ag / acov).where(acov > 0)
        df["aei_ag_mismatches"] = ag.astype("Int64")
        df["aei_a_coverage"] = acov.astype("Int64")
    elif "AEI_percent" in df.columns:
        idx = pd.to_numeric(df["AEI_percent"], errors="coerce") / 100.0
        df["aei_ag_mismatches"] = pd.NA
        df["aei_a_coverage"] = pd.NA
    else:
        raise KeyError(
            "AEI table needs either (AG_mismatches, A_coverage) or AEI_percent; "
            f"got {list(df.columns)}"
        )

    df["alu_editing_index"] = idx.astype(float)
    df["alu_editing_index_percent"] = (idx * 100.0).astype(float)
    df["aei_signal_to_noise"] = (
        pd.to_numeric(df["signal_to_noise"], errors="coerce")
        if "signal_to_noise" in df.columns else pd.NA
    )
    df["alu_editing_index_cohortz"] = (
        df.groupby("cohort")["alu_editing_index"].transform(_robust_z)
        if df["cohort"].notna().any() else _robust_z(df["alu_editing_index"])
    )

    cols = ["run_accession", "cohort", "alu_editing_index",
            "alu_editing_index_percent", "aei_ag_mismatches", "aei_a_coverage",
            "aei_signal_to_noise", "alu_editing_index_cohortz"]
    return df[cols].reset_index(drop=True)


# ===========================================================================
# Feature 2 — editing_neoantigen_burden (recoding I=G peptides -> shared engine)
# ===========================================================================
@dataclass
class RecodingSite:
    """One known CDS recoding editing site with enough context to build the
    altered (I=G) peptide deterministically.

    Fields
    ------
    site_id     : stable id (e.g. 'GRIA2_Q607R' or 'chr4_157336723').
    chrom, pos  : genomic coordinate of the edited adenosine (1-based, matching
                  the pipeline's editing_sites.tsv ``pos``).
    strand      : gene strand ('+'/'-'); the pipeline reports A-to-I on the
                  sense strand, so an edited '+' site has ref 'A' and a '-' site
                  ref 'T' on the genomic forward strand.
    gene        : gene symbol (annotation only).
    cds_window  : reference CODING-STRAND CDS nucleotide sequence, a window of
                  (2*CODON_FLANK + 1) codons centred on the recoded codon, IN
                  FRAME (len % 3 == 0, frame 0). At pilot time this is filled
                  from the genome FASTA + transcript CDS model
                  (build_recoding_catalog_from_genome); in tests it is supplied
                  directly.
    edit_offset : 0-based offset of the edited adenosine within cds_window
                  (coding strand); cds_window[edit_offset] must be 'A'.
    """
    site_id: str
    chrom: str
    pos: int
    strand: str
    gene: str
    cds_window: str
    edit_offset: int


def recoding_peptides(site: RecodingSite,
                      kmin: int = PEP_MIN, kmax: int = PEP_MAX) -> dict:
    """Generate the altered (I=G) peptides for one recoding site.

    Applies the A->G (I=G) edit at ``edit_offset`` on the coding-strand CDS
    window, translates the reference and edited windows in frame 0, confirms the
    edit is NONSYNONYMOUS (recoding), and tiles the kmin..kmax-mers of the
    EDITED protein that SPAN the recoded residue.

    Returns
    -------
    dict: {ref_aa, alt_aa, aa_change, is_recoding (bool), peptides (list[str])}.
    peptides is empty when the edit is synonymous or introduces a stop at the
    recoded residue (no spanning antigenic peptide).
    """
    ref_cds = site.cds_window.upper()
    edt_cds = apply_ItoG(ref_cds, site.edit_offset)
    ref_prot = translate(ref_cds, 0)
    edt_prot = translate(edt_cds, 0)
    residue_idx = site.edit_offset // 3
    ref_aa = ref_prot[residue_idx] if residue_idx < len(ref_prot) else "X"
    alt_aa = edt_prot[residue_idx] if residue_idx < len(edt_prot) else "X"
    is_recoding = (ref_aa != alt_aa) and (alt_aa != STOP) and (ref_aa != "X")
    peptides = (_kmers_spanning(edt_prot, residue_idx, kmin, kmax)
                if is_recoding else [])
    return {
        "site_id": site.site_id,
        "ref_aa": ref_aa,
        "alt_aa": alt_aa,
        "aa_change": f"{ref_aa}->{alt_aa}",
        "is_recoding": bool(is_recoding),
        "peptides": sorted(set(peptides)),
    }


def _edited_site_keys(editing_sites: pd.DataFrame,
                      freq_threshold: float) -> "dict[tuple, set[tuple[str, int]]]":
    """{(run_accession, cohort): set of (chrom, pos) edited above threshold}."""
    df = editing_sites.copy()
    df["edit_freq"] = pd.to_numeric(df["edit_freq"], errors="coerce").fillna(0.0)
    df = df[df["edit_freq"] >= freq_threshold]
    out: dict[tuple, set[tuple[str, int]]] = {}
    for (samp, cohort), g in df.groupby(["run_accession", "cohort"], sort=False):
        out[(str(samp), cohort)] = {
            (str(c), int(p)) for c, p in zip(g["chrom"], g["pos"])
        }
    return out


def compute_editing_neoantigen_burden(
    editing_sites: pd.DataFrame,
    recoding_catalog: Sequence[RecodingSite],
    hla_by_sample: "dict[str, Sequence[str]]",
    freq_threshold: float = EDIT_FREQ_THRESHOLD,
    rank_threshold: float = WEAK_BINDER_RANK,
) -> pd.DataFrame:
    """Per-sample recoding neoantigen burden via the SHARED engine.

    For each sample: intersect its per-site A-to-I calls (edit freq >=
    ``freq_threshold``) with the fixed recoding catalog; for every edited
    recoding site, build the altered (I=G) spanning peptides; pool the unique
    peptides across sites and hand them + the sample's HLA-I alleles to
    ``count_binders`` (the one shared MHCflurry engine). The result is the count
    of UNIQUE binder peptides — directly comparable to the other
    *_neoantigen_burden features.

    Parameters
    ----------
    editing_sites : long table of per-sample A-to-I calls, columns
        (run_accession, cohort, chrom, pos, strand, edit_freq, ...) — the
        pooled pipeline *.editing_sites.tsv.
    recoding_catalog : fixed list of known CDS recoding sites (RecodingSite).
    hla_by_sample : {run_accession: [allele, ...]} (up to 6 HLA-I alleles).
        Samples absent from this map get <NA> burden (HLA not typed).

    Returns
    -------
    DataFrame (run_accession, cohort, editing_neoantigen_burden,
    n_recoding_sites_edited, n_candidate_peptides). editing_neoantigen_burden is
    <NA> where HLA is missing.
    """
    # precompute each catalog site's peptides + genomic key once (fixed catalog)
    site_by_key: dict[tuple[str, int], RecodingSite] = {}
    peps_by_key: dict[tuple[str, int], list[str]] = {}
    for s in recoding_catalog:
        key = (str(s.chrom), int(s.pos))
        site_by_key[key] = s
        peps_by_key[key] = recoding_peptides(s)["peptides"]

    edited = _edited_site_keys(editing_sites, freq_threshold)

    # Enumerate ALL assayed samples (present in editing_sites), so a sample that
    # was assayed but has no above-threshold recoding edit gets burden 0 rather
    # than dropping out of the table. edited{} only holds samples with >=1
    # above-threshold site.
    assayed = (editing_sites[["run_accession", "cohort"]]
               .astype({"run_accession": str})
               .drop_duplicates())
    all_keys = [(str(s), c) for s, c in
                zip(assayed["run_accession"], assayed["cohort"])]

    rows = []
    for (samp, cohort) in all_keys:
        keys = edited.get((samp, cohort), set())
        hit_keys = [k for k in keys if k in site_by_key]
        peps: set[str] = set()
        for k in hit_keys:
            peps.update(peps_by_key[k])
        pep_list = sorted(peps)
        alleles = hla_by_sample.get(samp)
        burden = (count_binders(pep_list, list(alleles), rank_threshold=rank_threshold)
                  if alleles else pd.NA)
        rows.append({
            "run_accession": samp,
            "cohort": cohort,
            "editing_neoantigen_burden": burden,
            "n_recoding_sites_edited": len(hit_keys),
            "n_candidate_peptides": len(pep_list),
        })
    out = pd.DataFrame(rows, columns=["run_accession", "cohort",
                                      "editing_neoantigen_burden",
                                      "n_recoding_sites_edited",
                                      "n_candidate_peptides"])
    if "editing_neoantigen_burden" in out:
        out["editing_neoantigen_burden"] = out["editing_neoantigen_burden"].astype("Int64")
    return out


# ===========================================================================
# Recoding catalog IO — build once from REDIportal coords + genome (pilot), or
# load a prebuilt TSV. Kept light; the synthetic test supplies RecodingSites
# directly so the peptide logic is validated without a genome FASTA.
# ===========================================================================
def load_recoding_catalog(tsv_path: "str | Path") -> list[RecodingSite]:
    """Load a prebuilt recoding-site catalog TSV into RecodingSite objects.

    Expected columns: site_id, chrom, pos, strand, gene, cds_window, edit_offset.
    Such a catalog is built once from the REDIportal coding/recoding table
    (genomic coords + gene) plus a transcript CDS model + genome FASTA
    (build_recoding_catalog_from_genome), and cached under reference/ so runs
    are reproducible and offline.
    """
    df = pd.read_csv(tsv_path, sep="\t")
    req = {"site_id", "chrom", "pos", "strand", "gene", "cds_window", "edit_offset"}
    missing = req - set(df.columns)
    if missing:
        raise KeyError(f"recoding catalog missing columns: {sorted(missing)}")
    return [
        RecodingSite(
            site_id=str(r.site_id), chrom=str(r.chrom), pos=int(r.pos),
            strand=str(r.strand), gene=str(r.gene),
            cds_window=str(r.cds_window), edit_offset=int(r.edit_offset),
        )
        for r in df.itertuples(index=False)
    ]


def save_recoding_catalog(catalog: Sequence[RecodingSite], tsv_path: "str | Path") -> None:
    """Serialise a recoding catalog back to the load_recoding_catalog TSV format."""
    pd.DataFrame([{
        "site_id": s.site_id, "chrom": s.chrom, "pos": s.pos, "strand": s.strand,
        "gene": s.gene, "cds_window": s.cds_window, "edit_offset": s.edit_offset,
    } for s in catalog]).to_csv(tsv_path, sep="\t", index=False)


# ===========================================================================
# Convenience: build both feature columns and the contract matrix.
# ===========================================================================
def hla_map_from_table(hla_table: pd.DataFrame) -> "dict[str, list[str]]":
    """{run_accession: [6 HLA-I alleles]} from an hla_typing.parquet-style table."""
    cols = [c for c in ALLELE_COLS if c in hla_table.columns]
    out: dict[str, list[str]] = {}
    for _, r in hla_table.iterrows():
        alleles = [str(r[c]) for c in cols if pd.notna(r[c]) and str(r[c]).strip()]
        out[str(r["run_accession"])] = alleles
    return out


def build_editing_features(
    aei_data: "str | Path | pd.DataFrame",
    editing_sites: pd.DataFrame,
    recoding_catalog: Sequence[RecodingSite],
    hla_table: pd.DataFrame,
    cohort: "str | dict[str, str] | None" = None,
    freq_threshold: float = EDIT_FREQ_THRESHOLD,
) -> pd.DataFrame:
    """End-to-end: alu_editing_index + editing_neoantigen_burden per sample.

    Returns a tidy feature frame keyed on (run_accession, cohort) carrying the
    two NAMED feature columns plus their companions/QC. This is what lands in
    results/features/ (merged into antigen_features.parquet per the v2 contract).
    Outer-join so a sample with an AEI but no recoding sites (or vice versa) is
    retained with NA in the missing feature (graceful degradation, per v1).
    """
    aei = compute_alu_editing_index(aei_data, cohort=cohort)
    hla_map = hla_map_from_table(hla_table)
    burden = compute_editing_neoantigen_burden(
        editing_sites, recoding_catalog, hla_map, freq_threshold=freq_threshold)
    feat = aei.merge(burden, on=["run_accession", "cohort"], how="outer")
    return feat


__all__ = [
    "EDIT_FREQ_THRESHOLD", "PEP_MIN", "PEP_MAX", "CODON_FLANK",
    "BATCH_ROBUSTNESS_NOTE",
    "revcomp", "translate", "apply_ItoG",
    "RecodingSite", "recoding_peptides",
    "compute_alu_editing_index",
    "compute_editing_neoantigen_burden",
    "load_recoding_catalog", "save_recoding_catalog",
    "hla_map_from_table", "build_editing_features",
]
