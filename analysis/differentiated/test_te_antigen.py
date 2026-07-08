#!/usr/bin/env python
"""
test_te_antigen.py — logic + unit validation for the TE/ERV antigen module.

Runs on a SMALL SYNTHETIC locus-count + sequence example (no real-cohort
values are ever fabricated). It checks the whole path:

    Telescope-style per-locus counts + genomic sequences
        -> within-sample activity filter
        -> 6-frame ORF translation -> 8-11mer peptides
        -> SHARED antigen_core MHCflurry engine
        -> te_antigen_burden family/locus-resolved features.

Key design test: one synthetic ERV (LTR/ERVK) locus is constructed so that a
known HLA-A*02:01 immunodominant epitope (CMV pp65 NLVPMVATV) sits in-frame in
its sequence, so the module MUST recover >=1 ERV binder against an A*02:01
genotype and attribute it to the LTR + ERV buckets. A low-count locus must be
filtered out (activity gate). Family classification + ERV detection are checked
directly. The headline burden must equal count_binders() on the pooled peptide
set (proves we route through the shared engine identically).

Run:  cd analysis/differentiated
      PYTHONPATH="..:.." python test_te_antigen.py
"""
from __future__ import annotations

import sys
from pathlib import Path

# make analysis/ importable (for antigen_core) and this dir
_HERE = Path(__file__).resolve().parent
_ANALYSIS = _HERE.parent
for p in (str(_ANALYSIS), str(_HERE)):
    if p not in sys.path:
        sys.path.insert(0, p)

from Bio.Seq import Seq  # noqa: E402

import te_antigen as te  # noqa: E402
from antigen_core.mhc_binding import count_binders, WEAK_BINDER_RANK  # noqa: E402

ALLELES = ["HLA-A*02:01", "HLA-A*01:01", "HLA-B*07:02",
           "HLA-B*08:01", "HLA-C*07:01", "HLA-C*07:02"]

# A known strong A*02:01 epitope we will bury, in-frame, inside an ERV locus.
EPITOPE = "NLVPMVATV"     # CMV pp65 495-503


def _back_translate(pep: str) -> str:
    """Naive codon back-translation (one codon per residue) for a fixture."""
    codon = {
        'A': 'GCT', 'R': 'CGT', 'N': 'AAT', 'D': 'GAT', 'C': 'TGT',
        'Q': 'CAA', 'E': 'GAA', 'G': 'GGT', 'H': 'CAT', 'I': 'ATT',
        'L': 'CTT', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT', 'P': 'CCT',
        'S': 'TCT', 'T': 'ACT', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTT',
    }
    return "".join(codon[a] for a in pep)


def make_fixture():
    """Build synthetic Telescope counts + locus sequences + annotation."""
    # ERV locus: start codon + flank + epitope + flank + stop, all in frame 0
    erv_orf = ("ATG"
               + _back_translate("GADETRICH")     # flank residues
               + _back_translate(EPITOPE)         # the buried epitope
               + _back_translate("KLPWQASTV")     # more flank
               + "TAA")
    erv_seq = "GG" * 5 + erv_orf + "CC" * 5        # some UTR-ish padding

    # LINE-1 locus: a benign ORF (random-ish residues, no known epitope)
    line_seq = ("ATG" + _back_translate("GGSSGGSSGGSSGGSSGGSSGGSS") + "TAA")

    # SINE/Alu locus: short + expressed (should still yield some peptides)
    sine_seq = ("ATG" + _back_translate("PLPLPLPLPLPLPLPL") + "TGA")

    # A low-count locus that must be filtered out by the activity gate
    quiet_seq = ("ATG" + _back_translate(EPITOPE + EPITOPE) + "TAA")

    counts = {
        "ERVK_locus_1": 480.0,     # active, ERV
        "L1_locus_1":   250.0,     # active, LINE
        "Alu_locus_1":  120.0,     # active, SINE
        "ERVK_locus_quiet": 3.0,   # BELOW min_reads=10 -> filtered
    }
    locus_seqs = {
        "ERVK_locus_1": erv_seq,
        "L1_locus_1":   line_seq,
        "Alu_locus_1":  sine_seq,
        "ERVK_locus_quiet": quiet_seq,
    }
    import pandas as pd
    annotation = pd.DataFrame([
        {"locus_id": "ERVK_locus_1", "repeat_class": "LTR/ERVK",
         "chrom": "chr1", "start": 1000, "end": 1000 + len(erv_seq), "strand": "+"},
        {"locus_id": "L1_locus_1", "repeat_class": "LINE/L1",
         "chrom": "chr1", "start": 5000, "end": 5000 + len(line_seq), "strand": "+"},
        {"locus_id": "Alu_locus_1", "repeat_class": "SINE/Alu",
         "chrom": "chr2", "start": 200, "end": 200 + len(sine_seq), "strand": "+"},
        {"locus_id": "ERVK_locus_quiet", "repeat_class": "LTR/ERVK",
         "chrom": "chr3", "start": 10, "end": 10 + len(quiet_seq), "strand": "+"},
    ])
    return counts, locus_seqs, annotation


# ---------------------------------------------------------------------------
def test_family_classification():
    assert te.classify_family("LINE/L1") == "LINE"
    assert te.classify_family("SINE/Alu") == "SINE"
    assert te.classify_family("LTR/ERVK") == "LTR"
    assert te.classify_family("DNA/TcMar-Tigger") == "DNA"
    assert te.classify_family("Simple_repeat") == "OTHER"
    assert te.classify_family(None) == "OTHER"
    assert te.is_erv("LTR/ERVK") is True
    assert te.is_erv("LTR/ERV1") is True
    assert te.is_erv("LTR/Gypsy") is True
    assert te.is_erv("LINE/L1") is False    # LINE is not an ERV
    assert te.is_erv("SINE/Alu") is False
    print("  family classification + ERV detection: OK")


def test_activity_gate():
    counts = {"a": 480.0, "b": 3.0, "c": 0.0, "d": 15.0}
    active = te.select_expressed_loci(counts, min_reads=10.0, min_cpm=1.0)
    assert "a" in active and "d" in active
    assert "b" not in active            # below min_reads
    assert "c" not in active            # zero
    # ordered by descending count
    assert active[0] == "a"
    print("  within-sample activity gate:", active)


def test_orf_peptides_recover_epitope():
    _, locus_seqs, _ = make_fixture()
    peps = te.peptides_from_sequence(locus_seqs["ERVK_locus_1"], strand="+")
    assert EPITOPE in peps, "buried epitope not recovered from ORF translation"
    # every peptide is a valid 8-11mer with standard AAs
    assert all(8 <= len(p) <= 11 for p in peps)
    print(f"  ORF translation recovered epitope; {len(peps)} candidate peptides")


def test_burden_recovers_erv_binder():
    counts, locus_seqs, annotation = make_fixture()
    r = te.compute_te_antigen_burden(
        counts, locus_seqs, ALLELES, annotation=annotation)

    # quiet locus filtered -> 3 expressed loci
    assert r["te_antigen_n_expressed_loci"] == 3, r["te_antigen_n_expressed_loci"]
    # the buried A*02:01 epitope -> at least one binder overall and in ERV/LTR
    assert r["te_antigen_burden"] >= 1
    assert r["te_antigen_burden_ERV"] >= 1, "ERV binder not attributed"
    assert r["te_antigen_burden_LTR"] >= r["te_antigen_burden_ERV"]
    # top locus should be the ERV one that carries the epitope
    assert r["te_antigen_top_locus"] == "ERVK_locus_1", r["te_antigen_top_locus"]
    # binder loci QC
    assert r["te_antigen_n_binder_loci"] >= 1
    print("  burden features:",
          {k: v for k, v in r.items()
           if k not in ("locus_contributions", "scored")})


def test_headline_equals_shared_engine():
    """Headline burden must equal count_binders() on the pooled peptide set —
    proves we route through the shared engine with identical semantics."""
    counts, locus_seqs, annotation = make_fixture()
    active = te.select_expressed_loci(counts)
    # reconstruct the pooled candidate set EXACTLY as compute does — honouring
    # the annotation strand (3-frame for '+'/'-', 6-frame if unknown).
    strand_map = dict(zip(annotation["locus_id"], annotation["strand"]))
    pep_by_locus = te.peptides_by_locus(
        {l: locus_seqs[l] for l in active if l in locus_seqs},
        strand_map=strand_map)
    pooled = set().union(*pep_by_locus.values()) if pep_by_locus else set()
    expected = count_binders(sorted(pooled), ALLELES,
                             rank_threshold=WEAK_BINDER_RANK)
    r = te.compute_te_antigen_burden(counts, locus_seqs, ALLELES,
                                     annotation=annotation)
    assert r["te_antigen_burden"] == expected, (r["te_antigen_burden"], expected)
    print(f"  headline burden == shared count_binders() == {expected}: OK")


def test_row_and_table_shape():
    counts, locus_seqs, annotation = make_fixture()
    row = te.te_antigen_row("ERRTEST01", "gide2019", counts, locus_seqs,
                            ALLELES, annotation=annotation)
    assert row["run_accession"] == "ERRTEST01" and row["cohort"] == "gide2019"
    assert set(te._ROW_COLS) <= set(row.keys())
    df = te.build_te_antigen_table([row])
    assert list(df.columns) == te._ROW_COLS
    assert len(df) == 1
    print("  tidy row + table shape:", list(df.columns))


def test_empty_and_degenerate_inputs():
    # no expressed loci -> all-zero burden, no error
    r = te.compute_te_antigen_burden({"x": 1.0}, {"x": "ATGATG"}, ALLELES)
    assert r["te_antigen_burden"] == 0
    assert r["te_antigen_top_locus"] == ""
    # no HLA alleles -> zero (engine guards)
    counts, locus_seqs, ann = make_fixture()
    r2 = te.compute_te_antigen_burden(counts, locus_seqs, [], annotation=ann)
    assert r2["te_antigen_burden"] == 0
    print("  empty / degenerate-input guards: OK")


if __name__ == "__main__":
    print("TE/ERV antigen module validation:")
    test_family_classification()
    test_activity_gate()
    test_orf_peptides_recover_epitope()
    test_empty_and_degenerate_inputs()
    test_burden_recovers_erv_binder()
    test_headline_equals_shared_engine()
    test_row_and_table_shape()
    print("\nALL TE-ANTIGEN TESTS PASSED")
