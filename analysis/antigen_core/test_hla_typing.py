#!/usr/bin/env python
"""
test_hla_typing.py — logic validation for the HLA typing module.

Full arcasHLA typing runs on the Linux pipeline host during the pilot (needs
kallisto + the git-lfs IMGT/HLA reference). Here we validate the parsing and
heterozygosity logic — the part that must be exactly correct — on synthetic
genotype JSONs, so the module is trustworthy the moment real genotypes land.
NO real-sample allele calls are fabricated.
"""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

from hla_typing import (
    normalize_allele, is_heterozygous_locus, summarize_genotype,
    parse_arcashla_genotype, parse_optitype_result, build_hla_table,
)


def test_normalize():
    assert normalize_allele("HLA-A*02:01") == "A*02:01"
    assert normalize_allele("A*02:01:01:02") == "A*02:01"   # collapse to 2-field
    assert normalize_allele("A_02_01") == "A*02:01"
    assert normalize_allele("") == ""
    assert normalize_allele(None) == ""
    print("  normalize_allele: OK")


def test_het_locus():
    assert is_heterozygous_locus("A*02:01", "A*01:01") is True
    assert is_heterozygous_locus("A*02:01", "A*02:01") is False   # homozygous
    assert is_heterozygous_locus("A*02:01:01", "A*02:01:99") is False  # same 2-field
    assert is_heterozygous_locus("A*02:01", "") is False          # missing
    print("  is_heterozygous_locus: OK")


def test_summarize_fully_het():
    # heterozygous at all three loci -> Chowell favorable
    g = {"A": ["A*02:01", "A*01:01"],
         "B": ["B*07:02", "B*08:01"],
         "C": ["C*07:01", "C*07:02"]}
    row = summarize_genotype(g, "SAMPLE_HET", "gide2019", tool="synthetic")
    assert row["n_het_loci"] == 3
    assert row["HLA_I_heterozygous"] is True
    assert row["HLA_A_1"] == "A*02:01" and row["HLA_A_2"] == "A*01:01"
    print("  summarize (fully het): OK  ->", row["HLA_I_heterozygous"])


def test_summarize_one_homozygous():
    # homozygous at A (arcasHLA reports a single-element list) -> NOT fully het
    g = {"A": ["A*02:01"],                # homozygous, single element
         "B": ["B*07:02", "B*08:01"],
         "C": ["C*07:01", "C*07:02"]}
    row = summarize_genotype(g, "SAMPLE_HOM_A", "gide2019", tool="synthetic")
    assert row["HLA_A_1"] == "A*02:01" and row["HLA_A_2"] == "A*02:01"
    assert row["n_het_loci"] == 2
    assert row["HLA_I_heterozygous"] is False
    print("  summarize (homozygous A): OK  ->", row["HLA_I_heterozygous"])


def test_parse_arcashla(tmpdir: Path):
    # arcasHLA writes 3-field calls and may include non-class-I loci
    data = {"A": ["A*03:01:01", "A*11:01:01"],
            "B": ["B*35:01:01", "B*44:02:01"],
            "C": ["C*04:01:01", "C*05:01:01"],
            "DRB1": ["DRB1*15:01"]}   # should be ignored
    p = tmpdir / "SRRX.genotype.json"
    p.write_text(json.dumps(data))
    g = parse_arcashla_genotype(p)
    assert set(g) == {"A", "B", "C"}
    row = summarize_genotype(g, "SRRX", "hugo2016", tool="arcasHLA")
    assert row["HLA_I_heterozygous"] is True
    assert row["HLA_A_1"] == "A*03:01"
    print("  parse_arcashla_genotype + summarize: OK")


def test_parse_optitype(tmpdir: Path):
    tsv = tmpdir / "opt_result.tsv"
    tsv.write_text("\tA1\tA2\tB1\tB2\tC1\tC2\tReads\tObjective\n"
                   "0\tA*02:01\tA*02:01\tB*07:02\tB*08:01\tC*07:01\tC*07:02\t100\t95.0\n")
    g = parse_optitype_result(tsv)
    row = summarize_genotype(g, "SRRY", "liu2019", tool="OptiType")
    assert row["n_het_loci"] == 2   # homozygous A
    assert row["HLA_I_heterozygous"] is False
    print("  parse_optitype_result + summarize: OK")


def test_build_table():
    rows = [
        summarize_genotype({"A": ["A*02:01", "A*01:01"], "B": ["B*07:02", "B*08:01"],
                            "C": ["C*07:01", "C*07:02"]}, "S1", "gide2019", tool="synthetic"),
        summarize_genotype({"A": ["A*02:01"], "B": ["B*07:02", "B*08:01"],
                            "C": ["C*07:01", "C*07:02"]}, "S2", "gide2019", tool="synthetic"),
    ]
    df = build_hla_table(rows)
    assert list(df.columns)[:2] == ["run_accession", "cohort"]
    assert df["HLA_I_heterozygous"].tolist() == [True, False]
    print("  build_hla_table: OK\n", df[["run_accession", "n_het_loci",
                                          "HLA_I_heterozygous"]].to_string(index=False))


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        print("HLA typing logic validation:")
        test_normalize()
        test_het_locus()
        test_summarize_fully_het()
        test_summarize_one_homozygous()
        test_parse_arcashla(tmp)
        test_parse_optitype(tmp)
        test_build_table()
        print("\nALL HLA LOGIC TESTS PASSED")
