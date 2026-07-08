#!/usr/bin/env python
"""
test_snv_indel_neoantigen.py — unit validation for the SNV/indel neoantigen module.

Strategy (no fabrication of real-cohort burdens): build a SYNTHETIC variant set
whose mutant proteins are engineered to expose textbook HLA-A*02:01 immunodominant
epitopes, then assert the shared engine flags them as binders.

  * missense    : substitute one residue so a novel 9-mer == GILGFVFTL
                  (influenza M1 58-66, canonical A*02:01 strong binder).
  * frameshift  : a 1-bp deletion (translated from a real CDS) shifts the frame
                  to read NLVPMVATV (CMV pp65 495-503, canonical A*02:01 strong).

Also checks the mechanics that make the count meaningful: novelty filtering
(WT k-mers dropped), overlap windowing, translate(), and empty-input guards.

Run:  python analysis/baseline/test_snv_indel_neoantigen.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "antigen_core"))

import snv_indel_neoantigen as sn  # noqa: E402

ALLELE = ["HLA-A*02:01"]

# --- synthetic missense: create GILGFVFTL by an A->G substitution -----------
_MISS_WT = "MKTAYIAKQRAILGFVFTLQRSTVWYKKDNH"   # residue 11 (A) -> G yields ...GILGFVFTL...
MISSENSE = sn.Variant(gene="SYN_MISS", wt_protein=_MISS_WT,
                      variant_type="missense", protein_pos=11, alt_aa="G")

# --- synthetic frameshift: 1-bp deletion, translated from a real CDS --------
_PREFIX = "ATGAAAACTGCTTATATTGCTAAACAACGT"      # MKTAYIAKQR (in frame)
_TAIL = "AATCTTGTTCCTATGGTTGCTACTGTTTAA"        # -> NLVPMVATV* when in this frame
_WT_CDS = _PREFIX + "C" + _TAIL                 # reference frame (self tail)
_MUT_CDS = _PREFIX + _TAIL                      # 1-bp del -> frameshift
FRAMESHIFT = sn.Variant(gene="SYN_FS", wt_protein=sn.translate(_WT_CDS),
                        variant_type="frameshift",
                        mutant_protein=sn.translate(_MUT_CDS))


def test_translate():
    assert sn.translate("ATGAAA") == "MK"
    assert sn.translate("AATCTTGTTCCTATGGTTGCTACTGTTTAA") == "NLVPMVATV"  # stops at TAA
    assert sn.translate("ATGTAAAAA") == "M"                             # early stop
    print("  translate: OK")


def test_missense_creates_known_epitope():
    peps = sn.peptides_for_variant(MISSENSE)
    assert "GILGFVFTL" in peps, "missense did not generate the GILGFVFTL epitope"
    # every candidate must be novel (absent from WT protein) and in-length
    for p in peps:
        assert 8 <= len(p) <= 11
        assert p not in _MISS_WT, f"non-novel peptide leaked: {p}"
    print(f"  missense: GILGFVFTL present, {len(peps)} novel candidates")


def test_frameshift_creates_known_epitope():
    assert FRAMESHIFT.mutant_protein.endswith("NLVPMVATV")
    peps = sn.peptides_for_variant(FRAMESHIFT)
    assert "NLVPMVATV" in peps, "frameshift did not expose NLVPMVATV"
    print(f"  frameshift: NLVPMVATV present, {len(peps)} novel candidates")


def test_novelty_filter():
    # a silent/synonymous-like change producing an identical protein => no peptides
    v = sn.Variant(gene="SYN_SILENT", wt_protein=_MISS_WT,
                   variant_type="missense", protein_pos=11, alt_aa="A")  # A->A
    assert sn.peptides_for_variant(v) == [], "identical protein yielded peptides"
    print("  novelty filter: identical protein -> 0 candidates")


def test_known_peptides_are_binders():
    """The engine must flag both engineered epitopes as A*02:01 binders."""
    from mhc_binding import score_peptides  # shared engine
    scored = score_peptides(["GILGFVFTL", "NLVPMVATV"], ALLELE)
    assert not scored.empty
    for pep in ("GILGFVFTL", "NLVPMVATV"):
        row = scored[scored["peptide"] == pep].iloc[0]
        assert row["affinity_percentile"] <= 2.0, f"{pep} not even a weak binder"
        assert bool(row["is_strong"]), f"{pep} expected to be a strong binder"
    print("  engine: GILGFVFTL & NLVPMVATV are strong A*02:01 binders")


def test_burden_counts_binders():
    detail = sn.burden_detail([MISSENSE, FRAMESHIFT], ALLELE)
    print("  burden detail:", detail)
    # both engineered epitopes should surface as strong binders
    assert detail["n_strong_binders"] >= 2, "expected >=2 strong neoepitope binders"
    assert detail["n_weak_binders"] >= detail["n_strong_binders"]
    assert detail[sn.FEATURE_COL] == detail["n_weak_binders"]
    burden = sn.snv_indel_neoantigen_burden([MISSENSE, FRAMESHIFT], ALLELE)
    assert burden == detail["n_weak_binders"]
    assert isinstance(burden, int)


def test_empty_guards():
    assert sn.snv_indel_neoantigen_burden([], ALLELE) == 0
    assert sn.snv_indel_neoantigen_burden([MISSENSE], []) == 0
    assert sn.collect_peptides([]) == []
    print("  empty-input guards: OK")


if __name__ == "__main__":
    print("SNV/indel neoantigen module validation:")
    test_translate()
    test_missense_creates_known_epitope()
    test_frameshift_creates_known_epitope()
    test_novelty_filter()
    test_empty_guards()
    test_known_peptides_are_binders()
    test_burden_counts_binders()
    print("\nALL SNV/INDEL NEOANTIGEN TESTS PASSED")
