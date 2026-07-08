#!/usr/bin/env python
"""
test_intron_retention.py — logic + engine validation for the retained-intron
antigen module (analysis/differentiated/intron_retention.py).

What is validated (on SYNTHETIC inputs — no real-cohort feature values fabricated):
  1. retained_intron_load counting + depth-normalised fraction + within-cohort z.
  2. Reverse-complement / 3-frame translation helpers (deterministic).
  3. Cryptic-ORF peptide generation from a synthetic intron whose sequence is
     constructed so a KNOWN HLA-A*02:01 epitope (influenza M1 GILGFVFTL) is
     encoded straddling the exon/intron junction — the peptide must appear in
     the junction-spanning set.
  4. ir_neoantigen_burden through the SHARED MHCflurry engine: with the planted
     A*02:01 epitope the burden is > 0; a sample with no HLA typed is <NA>.
  5. build_ir_features end-to-end wiring on a tiny synthetic genome + SAF.

Run:  python analysis/differentiated/test_intron_retention.py
      (or via pytest)
"""
from __future__ import annotations

import os
import tempfile
from pathlib import Path

import pandas as pd

import intron_retention as ir  # path-style import (dir on sys.path when run direct)
try:  # pytest from repo root
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import intron_retention as ir  # noqa: F811
except Exception:
    pass


# --- reverse-translate a protein to a DNA coding sequence (test helper) ------
_AA2CODON = {
    "F": "TTT", "L": "CTG", "I": "ATT", "M": "ATG", "V": "GTG", "S": "TCG",
    "P": "CCG", "T": "ACG", "A": "GCG", "Y": "TAT", "H": "CAT", "Q": "CAG",
    "N": "AAT", "K": "AAA", "D": "GAT", "E": "GAG", "C": "TGT", "W": "TGG",
    "R": "CGT", "G": "GGT",
}


def back_translate(prot: str) -> str:
    return "".join(_AA2CODON[a] for a in prot)


# ---------------------------------------------------------------------------
# 1. retained_intron_load
# ---------------------------------------------------------------------------
def test_retained_intron_load():
    # tidy-wide matrix: 3 samples x 4 introns, threshold 0.10
    mat = pd.DataFrame({
        "run_accession": ["S1", "S2", "S3"],
        "cohort": ["c", "c", "c"],
        "g1__intron_1": [0.30, 0.05, 0.02],   # S1 retained
        "g1__intron_2": [0.12, 0.50, None],   # S1,S2 retained; S3 NA (not evaluated)
        "g2__intron_1": [0.00, 0.90, 0.15],   # S2,S3 retained
        "g2__intron_2": [0.08, 0.08, 0.11],   # S3 retained
    })
    out = ir.compute_retained_intron_load(mat, threshold=0.10).set_index("run_accession")

    assert out.loc["S1", "retained_intron_load"] == 2, out.loc["S1"]
    assert out.loc["S2", "retained_intron_load"] == 2
    assert out.loc["S3", "retained_intron_load"] == 2
    # S3 has one NA intron -> only 3 evaluated
    assert out.loc["S3", "n_introns_evaluated"] == 3
    assert out.loc["S1", "n_introns_evaluated"] == 4
    # depth-robust fraction
    assert abs(out.loc["S3", "retained_intron_fraction"] - 2 / 3) < 1e-9
    assert abs(out.loc["S1", "retained_intron_fraction"] - 2 / 4) < 1e-9
    # weighted sum: S1 retained introns are 0.30 + 0.12 = 0.42
    assert abs(out.loc["S1", "retained_intron_load_weighted"] - 0.42) < 1e-6
    # cohort z present (finite)
    assert out["retained_intron_load_cohortz"].notna().all()
    print("  [1] retained_intron_load OK:",
          out[["retained_intron_load", "n_introns_evaluated",
               "retained_intron_fraction"]].to_dict("index"))


# ---------------------------------------------------------------------------
# 2. sequence helpers
# ---------------------------------------------------------------------------
def test_seq_helpers():
    assert ir.revcomp("ACGT") == "ACGT"
    assert ir.revcomp("AAAC") == "GTTT"
    # translate a back-translated protein round-trips
    prot = "MGILGFVFTLK"
    assert ir.translate_frame(back_translate(prot), 0) == prot
    print("  [2] revcomp + translate round-trip OK")


# ---------------------------------------------------------------------------
# 3 + 4. cryptic-ORF peptides + burden with a planted A*02:01 epitope
# ---------------------------------------------------------------------------
EPITOPE = "GILGFVFTL"  # influenza M1 58-66, textbook HLA-A*02:01 strong binder


def _synthetic_genome_and_saf(tmpdir: Path):
    """Build a tiny + strand gene where the intron read-through encodes the
    epitope straddling the exon/intron junction.

    Layout on chrTEST (1-based):
        exon flank : positions 1..30   (10 codons of exon, in frame 0)
        intron     : positions 31..    (read-through continues the frame)
    We place the junction inside the epitope: the last exon codon + intron codons
    together spell ...GILGFVFTL... so GILGFVFTL is junction-spanning.
    """
    # exon: 9 filler codons (M + 8xA) then start of epitope's first residue 'G'
    exon_prot = "M" + "A" * 8 + "G"          # 10 codons = 30 nt, ends in G (epitope[0])
    exon_dna = back_translate(exon_prot)      # 30 nt, frame 0
    # intron: rest of epitope (ILGFVFTL) then a few codons then STOP
    intron_prot = EPITOPE[1:] + "KRKR"        # ILGFVFTL + padding
    intron_dna = back_translate(intron_prot) + "TAA"  # explicit stop
    genome_seq = exon_dna + intron_dna
    # so full frame-0 translation = M A..A G ILGFVFTL KRKR *  -> contains GILGFVFTL
    fa = tmpdir / "mini_genome.fa"
    fa.write_text(f">chrTEST\n{genome_seq}\n")

    intron_start = len(exon_dna) + 1          # 1-based first intronic base = 31
    intron_end = len(genome_seq)
    saf = tmpdir / "introns.saf"
    saf.write_text(
        "GeneID\tChr\tStart\tEnd\tStrand\n"
        f"gTEST__intron_1\tchrTEST\t{intron_start}\t{intron_end}\t+\n"
    )
    return fa, saf


def test_cryptic_orf_peptides_contains_epitope():
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        fa, saf = _synthetic_genome_and_saf(tmp)
        coords = ir.load_intron_saf(saf)
        genome = ir.GenomeFasta(fa)
        peps = ir.cryptic_orf_peptides(coords["gTEST__intron_1"], genome)
        genome.close()
    assert EPITOPE in peps, f"planted junction epitope {EPITOPE} not in candidate peptides"
    # all peptides in the 8-11mer window
    assert all(ir.PEP_MIN <= len(p) <= ir.PEP_MAX for p in peps)
    print(f"  [3] cryptic_orf_peptides OK: {len(peps)} candidates, epitope present")


def test_ir_neoantigen_burden_engine():
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        fa, saf = _synthetic_genome_and_saf(tmp)
        coords = ir.load_intron_saf(saf)
        genome = ir.GenomeFasta(fa)
        ir_long = pd.DataFrame({
            "run_accession": ["S_typed", "S_notyped"],
            "cohort": ["c", "c"],
            "intron_id": ["gTEST__intron_1", "gTEST__intron_1"],
            "IR_ratio": [0.40, 0.40],
        })
        hla_map = {"S_typed": ["HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:02"]}
        out = ir.compute_ir_neoantigen_burden(
            ir_long, coords, genome, hla_map, threshold=0.10).set_index("run_accession")
        genome.close()
    # planted A*02:01 epitope -> burden > 0 for the typed sample
    assert int(out.loc["S_typed", "ir_neoantigen_burden"]) >= 1, out.loc["S_typed"]
    # untyped sample -> NA (graceful degradation)
    assert pd.isna(out.loc["S_notyped", "ir_neoantigen_burden"])
    assert out.loc["S_typed", "n_candidate_peptides"] > 0
    print("  [4] ir_neoantigen_burden via shared engine OK:",
          f"typed burden={int(out.loc['S_typed','ir_neoantigen_burden'])}, "
          f"untyped=NA")


# ---------------------------------------------------------------------------
# 5. end-to-end build_ir_features
# ---------------------------------------------------------------------------
def test_build_ir_features_end_to_end():
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        fa, saf = _synthetic_genome_and_saf(tmp)
        ir_matrix = pd.DataFrame({
            "run_accession": ["S_typed"],
            "cohort": ["c"],
            "gTEST__intron_1": [0.40],
        })
        hla_table = pd.DataFrame([{
            "run_accession": "S_typed", "cohort": "c",
            "HLA_A_1": "A*02:01", "HLA_A_2": "A*01:01",
            "HLA_B_1": "B*07:02", "HLA_B_2": "B*08:01",
            "HLA_C_1": "C*07:02", "HLA_C_2": "C*07:01",
        }])
        feat = ir.build_ir_features(ir_matrix, saf, fa, hla_table, threshold=0.10)
    row = feat.set_index("run_accession").loc["S_typed"]
    assert row["retained_intron_load"] == 1
    assert int(row["ir_neoantigen_burden"]) >= 1
    # both NAMED feature columns present
    for col in ("retained_intron_load", "ir_neoantigen_burden"):
        assert col in feat.columns
    print("  [5] build_ir_features end-to-end OK:",
          {"retained_intron_load": int(row["retained_intron_load"]),
           "ir_neoantigen_burden": int(row["ir_neoantigen_burden"])})


if __name__ == "__main__":
    test_retained_intron_load()
    test_seq_helpers()
    test_cryptic_orf_peptides_contains_epitope()
    test_ir_neoantigen_burden_engine()
    test_build_ir_features_end_to_end()
    print("\nAll intron_retention module tests passed.")
