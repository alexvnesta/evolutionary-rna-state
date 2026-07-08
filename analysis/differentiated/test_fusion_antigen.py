#!/usr/bin/env python
"""
test_fusion_antigen.py — unit validation for the fusion-transcript antigen module.

Validates the caller-agnostic logic on SYNTHETIC in-frame fusion examples (no
real cohort data, no MHCflurry download required for the peptide-logic tests).
The one test that touches the shared MHCflurry engine asserts the burden is a
non-negative integer and that a fusion junction carrying a textbook A*02:01
epitope scores as a binder.

Run:  PYTHONPATH=<repo> python analysis/differentiated/test_fusion_antigen.py
"""
from __future__ import annotations

import sys
from pathlib import Path

# make both the package path and the module dir importable
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from analysis.differentiated import fusion_antigen as fa


# ---------------------------------------------------------------------------
# 1. Junction-crossing peptide enumeration
# ---------------------------------------------------------------------------
def test_window_peptides_cross_junction():
    # 5 residues from gene1, 5 from gene2; junction at index 5
    prot = "AAAAABBBBB"
    peps = fa._window_peptides(prot, junction_idx=5, k_min=8, k_max=11)
    # every returned 8-11mer must contain >=1 'A' and >=1 'B' (crosses junction)
    assert peps, "no junction peptides generated"
    for p in peps:
        assert "A" in p and "B" in p, f"{p} does not cross the junction"
        assert 8 <= len(p) <= 11
    # a pure-gene1 8mer must NOT appear
    assert "AAAAABBB"[:8] in peps  # sanity: this one does cross (3 B's)
    print(f"  window_peptides: {len(peps)} junction-crossing k-mers OK")


def test_marked_peptide_parsing():
    # Arriba-style peptide_sequence: '|' junction, '...' truncation, lower-case
    # frameshift region, trailing '*' stop
    marked = "MKLVA|rqtpgsseqf*"
    peps = fa.peptides_from_marked_peptide(marked, 8, 11)
    assert peps, "no peptides from marked peptide"
    # all peptides cross the junction: contain >=1 upper-left residue and >=1
    # right-side residue (right side upper-cased internally)
    left = "MKLVA"
    for p in peps:
        assert any(c in p for c in left), f"{p} missing gene1 side"
    # stop codon truncates: no peptide should extend past the '*'
    assert all("*" not in p for p in peps)
    print(f"  marked_peptide: {len(peps)} peptides, stop-truncated OK")


def test_marked_peptide_stop_before_junction():
    # stop codon at/left of junction -> no viable crossing peptide
    assert fa.peptides_from_marked_peptide("MK*LV|ABCDEFGH", 8, 11) == []
    assert fa.peptides_from_marked_peptide("NOJUNCTIONHERE", 8, 11) == []
    print("  stop-before-junction / no-junction -> [] OK")


# ---------------------------------------------------------------------------
# 2. Nucleotide fusion-transcript translation
# ---------------------------------------------------------------------------
def test_translate_fusion_transcript():
    # gene1 CDS 'ATG GCC GCC' (M A A) | gene2 'GAT GAA TTT' (D E F), frame 0
    tx = "ATGGCCGCC|GATGAATTT"
    prot, junc = fa.translate_fusion_transcript(tx, frame=0)
    assert prot == "MAADEF", prot
    assert junc == 3, junc  # 9 nt left / 3 = 3 aa
    print(f"  translate_fusion_transcript: {prot!r} junction@{junc} OK")


# ---------------------------------------------------------------------------
# 3. FusionCall in-frame logic + Arriba/STAR-Fusion parsing
# ---------------------------------------------------------------------------
def test_inframe_flag():
    assert fa.FusionCall("A", "B", reading_frame="in-frame").is_inframe()
    assert fa.FusionCall("A", "B", reading_frame="INFRAME").is_inframe()
    assert not fa.FusionCall("A", "B", reading_frame="out-of-frame").is_inframe()
    assert not fa.FusionCall("A", "B", reading_frame=".").is_inframe()
    print("  is_inframe flag OK")


def _write_synthetic_arriba(path: Path, junction_peptide: str) -> None:
    """Write a minimal Arriba fusions.tsv with one in-frame + one out-of-frame
    + one low-confidence row."""
    header = ["gene1", "gene2", "breakpoint1", "breakpoint2", "confidence",
              "reading_frame", "peptide_sequence", "fusion_transcript",
              "split_reads1", "split_reads2", "discordant_mates"]
    rows = [
        # in-frame, high confidence -> should produce peptides
        ["EML4", "ALK", "2:42522656", "2:29446394", "high", "in-frame",
         junction_peptide, "atggcc|gataaa", "12", "8", "5"],
        # out-of-frame, high -> counts to n_fusions but not n_inframe
        ["BCR", "ABL1", "22:x", "9:y", "high", "out-of-frame",
         "MKL|rqtp*", "", "6", "4", "2"],
        # low confidence -> filtered out entirely by default
        ["FOO", "BAR", "1:a", "1:b", "low", "in-frame",
         "AAAA|BBBBCCCC", "", "1", "0", "0"],
    ]
    with open(path, "w") as fh:
        fh.write("#" + "\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")


def test_read_arriba_and_confidence_filter(tmp: Path):
    p = tmp / "fusions.tsv"
    # junction peptide carrying the influenza M1 A*02:01 epitope GILGFVFTL on
    # the gene2 side of the breakpoint
    _write_synthetic_arriba(p, "MKLVATP|GILGFVFTLQEE")
    calls = fa.read_arriba(p)
    assert len(calls) == 2, f"confidence filter failed: {len(calls)} kept"  # low dropped
    assert calls[0].name == "EML4--ALK"
    assert calls[0].is_inframe() and not calls[1].is_inframe()
    assert calls[0].support == 25  # 12+8+5
    print(f"  read_arriba: {len(calls)} calls (low-conf dropped), support OK")
    return calls


def test_read_starfusion(tmp: Path):
    p = tmp / "starfusion.tsv"
    header = ["#FusionName", "LeftBreakpoint", "RightBreakpoint",
              "PROT_FUSION_TYPE", "FUSION_TRANSLATION", "CDS_LEFT_RANGE",
              "JunctionReadCount", "SpanningFragCount"]
    # 15-residue left CDS (45 nt -> aa_junction = 15) then gene2 side
    prot = "MKLVATPQRSTUVWX" + "GILGFVFTLQEE"  # note: 'U'/'X' dropped by engine
    rows = [["EML4--ALK", "2:a", "2:b", "INFRAME", prot, "1-45", "9", "4"],
            ["BCR--ABL1", "22:a", "9:b", "FRAMESHIFT", "MKLrqtp", "1-9", "3", "1"]]
    with open(p, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    calls = fa.read_starfusion(p)
    assert len(calls) == 2
    assert calls[0].is_inframe() and not calls[1].is_inframe()
    assert "|" in calls[0].peptide, "STAR-Fusion junction not marked"
    peps = fa.fusion_peptides(calls[0])
    assert peps, "no peptides from STAR-Fusion in-frame call"
    print(f"  read_starfusion: junction reconstructed, {len(peps)} peptides OK")


# ---------------------------------------------------------------------------
# 4. End-to-end burden through the SHARED engine (uses MHCflurry)
# ---------------------------------------------------------------------------
def test_end_to_end_burden(tmp: Path):
    p = tmp / "fusions_e2e.tsv"
    # gene2 side carries GILGFVFTL (influenza M1, textbook HLA-A*02:01 strong
    # binder) so at least one junction peptide MUST be a binder for an A*02:01
    # sample -> burden >= 1.
    _write_synthetic_arriba(p, "MKLVATP|GILGFVFTLQEE")
    hla = ["HLA-A*02:01", "HLA-A*01:01", "HLA-B*07:02",
           "HLA-B*08:01", "HLA-C*07:01", "HLA-C*07:02"]
    row = fa.fusion_features_for_sample(
        run_accession="ERRTEST01", cohort="gide2019",
        hla_alleles=hla, arriba_tsv=p, caller_version="2.4.0",
    )
    print("  end-to-end row:", {k: row[k] for k in
          ("n_fusions", "n_inframe_fusions", "fusion_neoantigen_burden",
           "fusion_neoantigen_burden_strong")})
    # ASSERTIONS required by the task:
    # (a) breakpoint peptides generated + (b) burden is an integer
    assert isinstance(row["fusion_neoantigen_burden"], int)
    assert isinstance(row["n_inframe_fusions"], int)
    assert row["n_fusions"] == 2           # 2 high-conf (low dropped)
    assert row["n_inframe_fusions"] == 1   # only EML4--ALK is in-frame
    assert row["fusion_neoantigen_burden"] >= 1, \
        "GILGFVFTL junction peptide should bind A*02:01"
    assert row["fusion_neoantigen_burden_strong"] >= 1
    assert row["fusion_neoantigen_burden"] >= row["fusion_neoantigen_burden_strong"]
    assert row["caller"] == "arriba" and row["caller_version"] == "2.4.0"
    print("  end-to-end burden OK")
    return row


def test_fixed_caller_requirement(tmp: Path):
    p = tmp / "x.tsv"
    _write_synthetic_arriba(p, "MKLVATP|GILGFVFTLQEE")
    # both callers -> error; neither -> error
    for kwargs in ({"arriba_tsv": p, "starfusion_tsv": p}, {}):
        try:
            fa.fusion_features_for_sample("R", "c", ["HLA-A*02:01"], **kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError("fixed-caller requirement not enforced")
    print("  fixed-caller requirement enforced OK")


def test_table_assembly():
    df = fa.build_fusion_feature_table([
        {"run_accession": "R1", "cohort": "gide2019", "n_fusions": 3,
         "n_inframe_fusions": 1, "fusion_neoantigen_burden": 2,
         "fusion_neoantigen_burden_strong": 1, "caller": "arriba",
         "caller_version": "2.4.0"},
    ])
    assert list(df.columns) == fa.FEATURE_COLS
    assert df.loc[0, "fusion_neoantigen_burden"] == 2
    print("  table assembly OK")


if __name__ == "__main__":
    import tempfile
    tmp = Path(tempfile.mkdtemp())
    print("== fusion_antigen unit tests ==")
    test_window_peptides_cross_junction()
    test_marked_peptide_parsing()
    test_marked_peptide_stop_before_junction()
    test_translate_fusion_transcript()
    test_inframe_flag()
    test_read_arriba_and_confidence_filter(tmp)
    test_read_starfusion(tmp)
    test_table_assembly()
    test_fixed_caller_requirement(tmp)
    print("-- MHCflurry engine test (loads models) --")
    test_end_to_end_burden(tmp)
    print("ALL FUSION ANTIGEN TESTS PASSED")
