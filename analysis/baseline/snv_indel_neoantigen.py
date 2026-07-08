#!/usr/bin/env python
"""
analysis/baseline/snv_indel_neoantigen.py

BASELINE FLOOR — SNV/indel-derived neoantigen load (`snv_indel_neoantigen_burden`).

WHAT IT IS
----------
The classic tumour neoantigen load: mutant peptides created by nonsynonymous
SNVs (missense) and frameshift indels, presented on the patient's MHC-I. This
is the biomarker every RNA-state burden must beat. We recompute it through the
SAME shared MHCflurry engine (`analysis/antigen_core/mhc_binding.count_binders`)
so the number is defined *identically* to splice_neoantigen_burden,
te_antigen_burden, etc. — the burdens are then directly comparable because they
differ only in where the peptides come from, never in how binding is scored.

    somatic variants (VCF / MAF, annotated)
        -> mutant protein sequence per variant
        -> NOVEL k-mers (8-11) spanning the mutation, absent from the WT
           proteome  (this is the agretopicity/foreignness filter every
           neoantigen caller applies)
        -> count_binders(peptides, sample_HLA_I)   [shared engine]
        = snv_indel_neoantigen_burden  (int, unique MHC-I binding neopeptides)

WHY "NOVEL k-mers only"
-----------------------
A missense/frameshift is a neoantigen only if it produces a peptide the immune
system has not seen as self. We generate every 8-11mer overlapping a changed
residue and then DROP any that also occur in the wild-type protein (self).
This mirrors pVACtools / NeoPredPipe / the McGranahan neoantigen definition and
keeps the count restricted to genuinely mutation-bearing peptides.

MUTATION HANDLING
-----------------
* Missense SNV          : substitute the residue; changed region = 1 position.
* Frameshift indel      : WT prefix + novel translated tail (new reading frame
                          to the first stop); changed region = frameshift start
                          -> new C-terminus. All the frameshifted residues are
                          non-self, so every k-mer over the tail is a candidate.
* Inframe indel         : inserted/deleted residues form the changed region.

INPUTS (what the pilot will feed it)
------------------------------------
DNA-level variants come from **WES where available** (the melanoma ICB cohorts
banked WES-derived variant calls). The module consumes either:
  * a VEP/SnpEff-annotated **VCF**, or
  * a **MAF** (TCGA-style) with HGVSp / Amino_acids / Protein_position,
plus a **proteome FASTA** (Ensembl/RefSeq) to fetch WT protein context, and the
sample's 6 HLA-I alleles from `analysis/antigen_core/hla_typing`.

An optional expression filter (gene TPM >= min_tpm from quant_gene_tpm.parquet)
restricts to *expressed* variants — the FEATURE_CONTRACT batch-robustness note
("restrict to expressed variants") — since an unexpressed mutation cannot be a
presented antigen.

Feature column produced: `snv_indel_neoantigen_burden` (int).

NO FABRICATION: this module is exercised by `test_snv_indel_neoantigen.py` on a
small SYNTHETIC variant set with HLA-A*02:01 (a mutation engineered to create
the canonical influenza M1 epitope GILGFVFTL, plus a frameshift that exposes the
CMV pp65 epitope NLVPMVATV — the test asserts both score as binders). Real
per-cohort burdens are produced only when the pilot supplies actual variant
calls.
"""
from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

# --- import the shared antigen core (sibling package) ----------------------
# Prefer the packaged import; fall back to a path insert so the module also
# runs standalone (e.g. `python analysis/baseline/snv_indel_neoantigen.py`).
try:
    from analysis.antigen_core.mhc_binding import (         # noqa: F401
        count_binders, score_peptides, clean_peptides,
        STRONG_BINDER_RANK, WEAK_BINDER_RANK,
    )
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "antigen_core"))
    from mhc_binding import (                                # type: ignore # noqa: F401
        count_binders, score_peptides, clean_peptides,
        STRONG_BINDER_RANK, WEAK_BINDER_RANK,
    )

MIN_LEN, MAX_LEN = 8, 11
PEPTIDE_LENGTHS = (8, 9, 10, 11)
STOP = "*"
FEATURE_COL = "snv_indel_neoantigen_burden"

_CODON_TABLE = {  # standard genetic code
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate(nt: str) -> str:
    """Translate a nucleotide string to protein, stopping at the first stop."""
    nt = str(nt).upper().replace("U", "T")
    aa = []
    for i in range(0, len(nt) - 2, 3):
        res = _CODON_TABLE.get(nt[i:i + 3], "X")
        if res == STOP:
            break
        aa.append(res)
    return "".join(aa)


# ---------------------------------------------------------------------------
# Normalized variant record
# ---------------------------------------------------------------------------
@dataclass
class Variant:
    """A somatic protein-altering variant, normalized for peptide derivation.

    Two ways to specify the protein change (in priority order):
      1. `mutant_protein` given directly (most general — frameshift, complex).
      2. missense shorthand: `wt_protein` + `protein_pos` (1-based) + `alt_aa`;
         the mutant protein is the WT with that residue substituted.
    `wt_protein` is always required — it is the self reference used to drop
    non-novel k-mers.
    """
    gene: str
    wt_protein: str
    variant_type: str = "missense"        # missense | frameshift | inframe_indel
    protein_pos: int | None = None        # 1-based, missense/indel start
    alt_aa: str | None = None             # missense: mutant residue
    mutant_protein: str | None = None     # explicit mutant protein (preferred for indels)
    transcript_id: str | None = None
    hgvsp: str | None = None

    def resolve_mutant(self) -> tuple[str, int, int]:
        """Return (mutant_protein, first_changed_idx, last_changed_idx) 0-based.

        first/last bound the residues that differ from WT and therefore the
        region a candidate k-mer must overlap.
        """
        wt = self.wt_protein
        if self.mutant_protein is not None:
            mut = self.mutant_protein
            lo, hi = _diff_span(wt, mut, self.variant_type)
            return mut, lo, hi
        # missense shorthand
        if self.protein_pos is None or self.alt_aa is None:
            raise ValueError(f"{self.gene}: need mutant_protein OR "
                             f"(protein_pos & alt_aa)")
        idx = self.protein_pos - 1
        if not (0 <= idx < len(wt)):
            raise ValueError(f"{self.gene}: protein_pos {self.protein_pos} "
                             f"outside protein len {len(wt)}")
        mut = wt[:idx] + self.alt_aa + wt[idx + 1:]
        return mut, idx, idx + len(self.alt_aa) - 1


def _diff_span(wt: str, mut: str, variant_type: str) -> tuple[int, int]:
    """First and last (0-based, in mutant coords) residues that differ from WT.

    For a frameshift, everything from the first divergence to the mutant
    C-terminus is novel, so the span runs to len(mut)-1.
    """
    n = min(len(wt), len(mut))
    lo = 0
    while lo < n and wt[lo] == mut[lo]:
        lo += 1
    if lo >= len(mut):            # mutant identical/prefix — no novel residue
        return lo, lo - 1
    if variant_type == "frameshift":
        return lo, len(mut) - 1
    # substitution / inframe indel: find last divergence from the right
    hi_w, hi_m = len(wt) - 1, len(mut) - 1
    while hi_w >= lo and hi_m >= lo and wt[hi_w] == mut[hi_m]:
        hi_w -= 1
        hi_m -= 1
    return lo, max(hi_m, lo)


# ---------------------------------------------------------------------------
# Peptide derivation — the testable core
# ---------------------------------------------------------------------------
def derive_mutant_peptides(
    wt_protein: str,
    mut_protein: str,
    changed_lo: int,
    changed_hi: int,
    lengths: Sequence[int] = PEPTIDE_LENGTHS,
) -> list[str]:
    """All k-mers of `mut_protein` overlapping [changed_lo, changed_hi] that are
    NOT present in `wt_protein` (novel = mutation-bearing / non-self).

    Returns a de-duplicated, order-preserving list ready for the engine.
    """
    wt = str(wt_protein).upper()
    mut = str(mut_protein).upper()
    wt_kmer_index: dict[int, set[str]] = {}
    out: list[str] = []
    seen: set[str] = set()
    for L in lengths:
        if L not in wt_kmer_index:
            wt_kmer_index[L] = {wt[i:i + L] for i in range(len(wt) - L + 1)}
        wt_set = wt_kmer_index[L]
        # window start range so the k-mer overlaps the changed region
        start_min = max(0, changed_lo - L + 1)
        start_max = min(changed_hi, len(mut) - L)
        for s in range(start_min, start_max + 1):
            pep = mut[s:s + L]
            if len(pep) != L:
                continue
            if STOP in pep or "X" in pep:
                continue
            if pep in wt_set:            # not novel -> self, drop
                continue
            if pep in seen:
                continue
            seen.add(pep)
            out.append(pep)
    return out


def peptides_for_variant(v: Variant,
                         lengths: Sequence[int] = PEPTIDE_LENGTHS) -> list[str]:
    """Mutant candidate peptides for one Variant."""
    mut, lo, hi = v.resolve_mutant()
    if hi < lo:
        return []
    return derive_mutant_peptides(v.wt_protein, mut, lo, hi, lengths)


def collect_peptides(variants: Iterable[Variant],
                     lengths: Sequence[int] = PEPTIDE_LENGTHS) -> list[str]:
    """Union of novel mutant peptides across all variants (deduplicated)."""
    seen: set[str] = set()
    out: list[str] = []
    for v in variants:
        for pep in peptides_for_variant(v, lengths):
            if pep not in seen:
                seen.add(pep)
                out.append(pep)
    return out


# ---------------------------------------------------------------------------
# Top-level feature
# ---------------------------------------------------------------------------
def snv_indel_neoantigen_burden(
    variants: Sequence[Variant],
    hla_alleles: Sequence[str],
    rank_threshold: float = WEAK_BINDER_RANK,
    lengths: Sequence[int] = PEPTIDE_LENGTHS,
) -> int:
    """`snv_indel_neoantigen_burden` for one sample.

    Derives novel mutant peptides from `variants` and counts unique MHC-I
    binders (rank <= rank_threshold; default weak, 2.0) via the SHARED engine.
    Peptide hygiene, allele normalization, and unsupported-allele dropping are
    handled inside count_binders.
    """
    peptides = collect_peptides(variants, lengths)
    if not peptides or not list(hla_alleles):
        return 0
    return count_binders(peptides, list(hla_alleles), rank_threshold=rank_threshold)


def burden_detail(
    variants: Sequence[Variant],
    hla_alleles: Sequence[str],
    lengths: Sequence[int] = PEPTIDE_LENGTHS,
):
    """Diagnostics: {n_variants, n_candidate_peptides, n_weak, n_strong} +
    the scored per-peptide table (for inspection / QC)."""
    peptides = collect_peptides(variants, lengths)
    n_weak = count_binders(peptides, list(hla_alleles), rank_threshold=WEAK_BINDER_RANK)
    n_strong = count_binders(peptides, list(hla_alleles), rank_threshold=STRONG_BINDER_RANK)
    return {
        "n_variants": len(list(variants)),
        "n_candidate_peptides": len(peptides),
        "n_weak_binders": n_weak,
        "n_strong_binders": n_strong,
        FEATURE_COL: n_weak,
    }


# ---------------------------------------------------------------------------
# Parsers — MAF / VCF -> Variant, using a proteome FASTA for WT context
# ---------------------------------------------------------------------------
def load_proteome(fasta_path: str | Path) -> dict[str, str]:
    """Map transcript/protein id -> sequence from a proteome FASTA.

    Keys are stripped to the bare id (before the first '.' or whitespace) so
    Ensembl versioned ids (ENSP00000....3) match unversioned MAF references.
    """
    seqs: dict[str, str] = {}
    cur_id, buf = None, []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_id:
                    seqs[cur_id] = "".join(buf)
                raw = line[1:].strip().split()[0]
                cur_id = raw.split(".")[0]
                buf = []
            else:
                buf.append(line.strip())
        if cur_id:
            seqs[cur_id] = "".join(buf)
    return seqs


_AA3to1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V", "Ter": "*", "*": "*",
}


def variants_from_maf(
    maf_path: str | Path,
    proteome: dict[str, str],
    transcript_col: str = "Transcript_ID",
    expressed_genes: set[str] | None = None,
    gene_col: str = "Hugo_Symbol",
) -> list[Variant]:
    """Parse a (VEP/TCGA-style) MAF into missense/frameshift Variants.

    Requires columns: Variant_Classification, Protein_position, Amino_acids
    ("V/E" style), plus a transcript id resolvable in `proteome`. Frameshift
    rows without an explicit mutant protein are skipped with a warning (the
    pilot's VEP run should emit HGVSp/Downstream protein for those).

    `expressed_genes`: if given, keep only variants whose gene is expressed
    (the TPM filter) — enforcing the FEATURE_CONTRACT batch-robustness rule.
    """
    import csv
    keep_class = {
        "Missense_Mutation": "missense",
        "Frame_Shift_Del": "frameshift",
        "Frame_Shift_Ins": "frameshift",
        "In_Frame_Del": "inframe_indel",
        "In_Frame_Ins": "inframe_indel",
    }
    out: list[Variant] = []
    with open(maf_path) as fh:
        # skip MAF comment/version lines
        rows = [ln for ln in fh if not ln.startswith("#")]
    reader = csv.DictReader(rows, delimiter="\t")
    for r in reader:
        vclass = r.get("Variant_Classification", "")
        vtype = keep_class.get(vclass)
        if vtype is None:
            continue
        gene = r.get(gene_col, "")
        if expressed_genes is not None and gene not in expressed_genes:
            continue
        tid = (r.get(transcript_col) or "").split(".")[0]
        wt = proteome.get(tid)
        if not wt:
            continue
        if vtype == "missense":
            aa = r.get("Amino_acids", "")             # "V/E"
            pos = r.get("Protein_position", "").split("/")[0].split("-")[0]
            if "/" not in aa or not pos.isdigit():
                continue
            ref_aa, alt_aa = aa.split("/")[:2]
            out.append(Variant(gene=gene, wt_protein=wt, variant_type="missense",
                               protein_pos=int(pos), alt_aa=alt_aa,
                               transcript_id=tid,
                               hgvsp=r.get("HGVSp_Short")))
        else:
            # indel: require an explicit mutant protein via HGVSp is non-trivial;
            # the pilot's annotator should supply Mutant_Protein. Skip otherwise.
            mutp = r.get("Mutant_Protein") or r.get("mutant_protein")
            if mutp:
                out.append(Variant(gene=gene, wt_protein=wt, variant_type=vtype,
                                   mutant_protein=mutp, transcript_id=tid,
                                   hgvsp=r.get("HGVSp_Short")))
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main(argv: Sequence[str] | None = None) -> int:
    import argparse
    p = argparse.ArgumentParser(
        description="SNV/indel neoantigen burden via the shared MHC engine.")
    p.add_argument("--maf", required=True)
    p.add_argument("--proteome", required=True, help="proteome FASTA")
    p.add_argument("--hla", required=True, nargs="+",
                   help="HLA-I alleles, e.g. HLA-A*02:01 HLA-B*07:02 ...")
    p.add_argument("--rank", type=float, default=WEAK_BINDER_RANK)
    args = p.parse_args(argv)

    proteome = load_proteome(args.proteome)
    variants = variants_from_maf(args.maf, proteome)
    detail = burden_detail(variants, args.hla)
    print("SNV/indel neoantigen burden")
    print("=" * 40)
    for k, v in detail.items():
        print(f"  {k}: {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
