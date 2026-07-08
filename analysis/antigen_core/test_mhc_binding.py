#!/usr/bin/env python
"""
test_mhc_binding.py — unit validation for the shared MHC-I binding engine.

Sanity check: well-characterized HLA-A*02:01 immunodominant epitopes (influenza
M1 GILGFVFTL, CMV pp65 NLVPMVATV) must score as strong binders and must rank
far better than random decapeptides against A*02:01. This confirms the engine
is wired up and returning biologically sensible ranks before any antigen module
depends on it.
"""
from __future__ import annotations

import random

import mhc_binding as mb

# Textbook HLA-A*02:01 immunodominant epitopes (9-mers)
KNOWN_A0201 = ["GILGFVFTL",   # influenza M1 58-66
               "NLVPMVATV",   # CMV pp65 495-503
               "GLCTLVAML"]   # EBV BMLF1 259-267
ALLELE = ["HLA-A*02:01"]


def random_decamers(n: int, seed: int = 0) -> list[str]:
    rng = random.Random(seed)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    return ["".join(rng.choice(aa) for _ in range(10)) for _ in range(n)]


def test_cleaning():
    peps = mb.clean_peptides(
        ["GILGFVFTL", "gilgfvftl",         # dupe (case) -> 1
         "TOOSHORT1", "X" * 9,             # non-standard AA -> dropped
         "ACDEFGHIKLMNOP", "SHORT", ""]    # too long / too short -> dropped
    )
    assert "GILGFVFTL" in peps
    assert peps.count("GILGFVFTL") == 1
    assert all(set(p) <= set("ACDEFGHIKLMNPQRSTVWY") for p in peps)
    print("  clean_peptides:", peps)


def test_allele_norm():
    assert mb.normalize_allele("A*02:01") == "HLA-A*02:01"
    assert mb.normalize_allele("A0201") == "HLA-A*02:01"
    assert mb.normalize_allele("HLA-A*02:01") == "HLA-A*02:01"
    print("  normalize_allele: OK")


def test_known_vs_random():
    decoys = random_decamers(40)
    scored = mb.score_peptides(KNOWN_A0201 + decoys, ALLELE)
    assert not scored.empty, "engine returned no rows"

    known = scored[scored["peptide"].isin(KNOWN_A0201)]
    decoy = scored[~scored["peptide"].isin(KNOWN_A0201)]

    med_known = known["affinity_percentile"].median()
    med_decoy = decoy["affinity_percentile"].median()
    print(f"  median affinity %rank  known={med_known:.3f}  decoy={med_decoy:.2f}")
    print(known[["peptide", "affinity_nM", "affinity_percentile",
                 "presentation_score", "is_strong"]].to_string(index=False))

    # known epitopes must rank far better (lower) than random decoys
    assert med_known < med_decoy, "known binders did not outrank random decoys"
    # each known epitope should be at least a weak binder (rank <= 2.0)
    assert (known["affinity_percentile"] <= 2.0).all(), \
        "a known A*02:01 epitope failed the weak-binder threshold"
    # GILGFVFTL / NLVPMVATV are canonical strong binders
    assert known["is_strong"].sum() >= 2, "expected >=2 strong binders among knowns"


def test_count_binders():
    decoys = random_decamers(40, seed=1)
    peps = KNOWN_A0201 + decoys
    n_weak = mb.count_binders(peps, ALLELE, rank_threshold=mb.WEAK_BINDER_RANK)
    n_strong = mb.count_binders(peps, ALLELE, rank_threshold=mb.STRONG_BINDER_RANK)
    counts = mb.binder_counts(peps, ALLELE)
    print(f"  count_binders  weak={n_weak}  strong={n_strong}  detail={counts}")
    assert n_strong <= n_weak
    assert n_weak >= 2       # at least the canonical strong binders
    assert counts["n_scored"] == len(mb.clean_peptides(peps))


def test_empty_inputs():
    assert mb.count_binders([], ALLELE) == 0
    assert mb.count_binders(KNOWN_A0201, []) == 0
    assert mb.score_peptides([], ALLELE).empty
    print("  empty-input guards: OK")


if __name__ == "__main__":
    print("MHC-I binding engine validation:")
    test_cleaning()
    test_allele_norm()
    test_empty_inputs()
    test_known_vs_random()
    test_count_binders()
    print("\nALL MHC ENGINE TESTS PASSED")
