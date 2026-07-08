#!/usr/bin/env python
"""
test_rna_editing.py — logic + unit validation for the RNA-editing module.

Runs on SMALL SYNTHETIC AEI + editing-site + recoding-catalog fixtures (no
real-cohort values are ever fabricated). It checks the whole path:

    pipeline AEI table (compute_aei.py cols)
        -> alu_editing_index in [0,1]  (+ percent, QC, cohort z)
    pipeline editing_sites.tsv (per-site A-to-I calls)
        + fixed REDIportal recoding catalog
        -> apply I=G edit -> nonsynonymous recoding check
        -> altered 8-11mers spanning the recoded residue
        -> SHARED antigen_core MHCflurry engine
        -> editing_neoantigen_burden.

Key design test: one synthetic recoding site is constructed so that applying
the A-to-I(=G) edit turns the local CDS window into one whose translation
CONTAINS a known HLA-A*02:01 immunodominant epitope (CMV pp65 NLVPMVATV) that
is NOT present in the reference (germline) translation. The module must:
  * recover >=1 editing binder against an A*02:01 genotype;
  * make the headline editing_neoantigen_burden EQUAL count_binders() on the
    pooled altered-peptide set (proves we route through the SHARED engine);
  * count that site only when its per-sample edit frequency clears the
    threshold (depth/batch-robust gating), and never for a synonymous edit.

Run:  cd analysis/differentiated
      PYTHONPATH=".." python test_rna_editing.py
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

import pandas as pd  # noqa: E402

import rna_editing as re_mod  # noqa: E402
from antigen_core.mhc_binding import count_binders, WEAK_BINDER_RANK  # noqa: E402

ALLELES = ["HLA-A*02:01", "HLA-A*01:01", "HLA-B*07:02",
           "HLA-B*08:01", "HLA-C*07:01", "HLA-C*07:02"]

# A known strong A*02:01 epitope. We will engineer a recoding site whose EDITED
# window translates to a protein CONTAINING this peptide, spanning the edit.
EPITOPE = "NLVPMVATV"     # CMV pp65 495-503


_CODON = {
    'A': 'GCT', 'R': 'CGT', 'N': 'AAT', 'D': 'GAT', 'C': 'TGT',
    'Q': 'CAA', 'E': 'GAA', 'G': 'GGT', 'H': 'CAT', 'I': 'ATT',
    'L': 'CTT', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT', 'P': 'CCT',
    'S': 'TCT', 'T': 'ACT', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTT',
}


def _back_translate(pep: str) -> str:
    """Naive codon back-translation (one codon per residue) for a fixture."""
    return "".join(_CODON[a] for a in pep)


# ---------------------------------------------------------------------------
# Fixture construction
# ---------------------------------------------------------------------------
def make_recoding_site() -> re_mod.RecodingSite:
    """Build a recoding site whose I=G edit CREATES the A*02:01 epitope.

    Strategy: the recoded residue is the central residue of the epitope. We take
    NLVPMVATV, whose central residue (index 4) is 'M' (codon ATG). We choose a
    reference residue whose codon differs from ATG by exactly one A->G change:
    'I' = ATA -> (A->G at the 3rd base) -> ATG = 'M'. So reference residue = I,
    edited residue = M, and the edited window translates to ...NLVP M VATV...
    containing the epitope, while the reference translates to ...NLVP I VATV...
    which does NOT contain it. This exercises apply_ItoG + nonsynonymous check +
    spanning-kmer tiling end to end.
    """
    center = len(EPITOPE) // 2          # index 4 -> 'M'
    assert EPITOPE[center] == "M"
    ref_res = "I"                       # ATA differs from ATG (M) by one A->G
    # build residue windows with CODON_FLANK codons of filler each side so all
    # 8-11mers spanning the recoded residue can be tiled.
    flank = re_mod.CODON_FLANK
    left_fill = "G" * flank             # glycine filler (codon GGT, no stops)
    right_fill = "G" * flank
    # edited protein around the edit: left_fill + EPITOPE(with M center) + right_fill
    # reference protein: same but center residue = ref_res (I)
    epi_left = EPITOPE[:center]
    epi_right = EPITOPE[center + 1:]
    ref_res_window = left_fill + epi_left + ref_res + epi_right + right_fill
    # back-translate the REFERENCE window, then force the recoded codon to 'ATA'
    # (Ile) so the A->G(=I) edit at its 3rd base yields 'ATG' (Met). Our generic
    # back-translation maps I->ATT, which would not recode on an A->G edit.
    cds = list(_back_translate(ref_res_window))
    recoded_codon_idx = flank + center
    base = recoded_codon_idx * 3
    cds[base:base + 3] = list("ATA")
    cds = "".join(cds)
    edit_offset = base + 2   # 3rd base of the recoded codon (ATA -> ATG)
    assert cds[edit_offset] == "A", cds[edit_offset]
    return re_mod.RecodingSite(
        site_id="SYN_I2M_epitope",
        chrom="chrTEST", pos=1000, strand="+",
        gene="SYNGENE", cds_window=cds, edit_offset=edit_offset,
    )


def make_synonymous_site() -> re_mod.RecodingSite:
    """A recoding-catalog entry whose edit is SYNONYMOUS (must yield 0 peptides).

    Leucine CTA -> (A->G at 3rd base) -> CTG, still Leucine. So the edit changes
    the codon but not the amino acid; recoding_peptides must return is_recoding
    False and an empty peptide list.
    """
    flank = re_mod.CODON_FLANK
    window_res = "G" * flank + "L" + "G" * flank
    cds = list(_back_translate(window_res))
    # CTT (our back-translation of L) -> make it CTA so A->G gives CTG (still L)
    codon_idx = flank
    base = codon_idx * 3
    cds[base:base + 3] = list("CTA")
    cds = "".join(cds)
    edit_offset = base + 2   # 3rd base 'A' -> 'G'
    assert cds[edit_offset] == "A"
    return re_mod.RecodingSite(
        site_id="SYN_synonymous", chrom="chrTEST", pos=2000, strand="+",
        gene="SYNGENE2", cds_window=cds, edit_offset=edit_offset,
    )


def make_aei_table() -> pd.DataFrame:
    """A 3-sample AEI table in the pipeline's compute_aei.py column format."""
    return pd.DataFrame({
        "sample":        ["S_hi", "S_lo", "S_mid"],
        "AEI_percent":   [1.500000, 0.200000, 0.800000],
        "AG_mismatches": [15000, 2000, 8000],
        "A_coverage":    [1_000_000, 1_000_000, 1_000_000],
        "signal_to_noise": [30.0, 4.0, 12.0],
        "noise_floor_percent": [0.05, 0.05, 0.05],
    })


def make_editing_sites() -> pd.DataFrame:
    """Per-sample A-to-I calls (editing_sites.tsv format), long over samples.

    S_hi   edits the epitope site at freq 0.40 (above threshold)  -> counts
    S_lo   edits the epitope site at freq 0.05 (below threshold)  -> excluded
    S_mid  edits the epitope site at 0.30 AND the synonymous site -> 1 recoding
           hit (synonymous contributes no peptide).
    """
    rows = [
        # run_accession, cohort, chrom, pos, strand, ref, alt, edit_freq, edited_reads, coverage
        ("S_hi",  "synth", "chrTEST", 1000, "+", "A", "G", 0.40, 40, 100),
        ("S_lo",  "synth", "chrTEST", 1000, "+", "A", "G", 0.05,  5, 100),
        ("S_mid", "synth", "chrTEST", 1000, "+", "A", "G", 0.30, 30, 100),
        ("S_mid", "synth", "chrTEST", 2000, "+", "A", "G", 0.50, 50, 100),
        # a non-recoding Alu site nobody in the catalog knows about (ignored)
        ("S_hi",  "synth", "chrTEST", 9999, "+", "A", "G", 0.90, 90, 100),
    ]
    return pd.DataFrame(rows, columns=[
        "run_accession", "cohort", "chrom", "pos", "strand", "ref", "alt",
        "edit_freq", "edited_reads", "coverage"])


def make_hla_table() -> pd.DataFrame:
    return pd.DataFrame([{
        "run_accession": s, "cohort": "synth",
        "HLA_A_1": "A*02:01", "HLA_A_2": "A*01:01",
        "HLA_B_1": "B*07:02", "HLA_B_2": "B*08:01",
        "HLA_C_1": "C*07:01", "HLA_C_2": "C*07:02",
    } for s in ("S_hi", "S_lo", "S_mid")])


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
def test_sequence_utils():
    assert re_mod.revcomp("ATGC") == "GCAT"
    assert re_mod.translate("ATGGGT") == "MG"
    # apply_ItoG: A->G at offset, and it must reject a non-A base
    assert re_mod.apply_ItoG("ATA", 2) == "ATG"
    try:
        re_mod.apply_ItoG("ATG", 2)   # base at offset 2 is 'G', not 'A'
    except ValueError:
        pass
    else:
        raise AssertionError("apply_ItoG should reject a non-'A' edit base")
    print("  sequence utils (revcomp / translate / apply_ItoG guard): OK")


def test_alu_editing_index():
    aei = re_mod.compute_alu_editing_index(make_aei_table(), cohort="synth")
    print(aei[["run_accession", "alu_editing_index", "alu_editing_index_percent",
               "aei_signal_to_noise", "alu_editing_index_cohortz"]].to_string(index=False))
    # index must be a ratio in [0,1]
    assert aei["alu_editing_index"].between(0, 1).all(), "AEI outside [0,1]"
    # recomputed from raw counts: 15000/1e6 = 0.015
    hi = aei.set_index("run_accession").loc["S_hi", "alu_editing_index"]
    assert abs(hi - 0.015) < 1e-9, hi
    # percent column consistent
    assert abs(aei.set_index("run_accession").loc["S_hi", "alu_editing_index_percent"]
               - 1.5) < 1e-6
    # highest-editing sample has the highest cohort z
    top = aei.sort_values("alu_editing_index_cohortz").iloc[-1]["run_accession"]
    assert top == "S_hi", top
    print("  alu_editing_index in [0,1], raw-count recompute, cohort-z ordering: OK")


def test_recoding_peptides_create_epitope():
    site = make_recoding_site()
    res = re_mod.recoding_peptides(site)
    print(f"  recoding {site.site_id}: {res['aa_change']}  "
          f"is_recoding={res['is_recoding']}  n_pep={len(res['peptides'])}")
    assert res["is_recoding"], "engineered I->M edit should be nonsynonymous"
    assert res["aa_change"] == "I->M", res["aa_change"]
    # the epitope must appear among the altered spanning peptides, and NOT be
    # derivable from the reference translation
    assert EPITOPE in res["peptides"], "edited peptides must contain the epitope"
    ref_prot = re_mod.translate(site.cds_window, 0)
    assert EPITOPE not in ref_prot, "epitope must be absent from reference protein"
    # synonymous site yields no peptides
    syn = re_mod.recoding_peptides(make_synonymous_site())
    assert not syn["is_recoding"] and syn["peptides"] == [], syn
    print("  I=G edit creates the epitope; synonymous edit yields none: OK")


def test_editing_neoantigen_burden():
    catalog = [make_recoding_site(), make_synonymous_site()]
    sites = make_editing_sites()
    hla = re_mod.hla_map_from_table(make_hla_table())
    burden = re_mod.compute_editing_neoantigen_burden(
        sites, catalog, hla, freq_threshold=re_mod.EDIT_FREQ_THRESHOLD)
    b = burden.set_index("run_accession")
    print(burden.to_string(index=False))

    # S_hi edits the epitope site above threshold -> >=1 binder
    assert b.loc["S_hi", "editing_neoantigen_burden"] >= 1, "S_hi should have >=1 editing binder"
    # S_hi hit exactly one catalog recoding site (the 0.90 site is not in catalog)
    assert b.loc["S_hi", "n_recoding_sites_edited"] == 1, b.loc["S_hi"]
    # S_lo edits it only at 0.05 (below threshold) -> no recoding site, 0 burden
    assert b.loc["S_lo", "n_recoding_sites_edited"] == 0
    assert b.loc["S_lo", "editing_neoantigen_burden"] == 0
    # S_mid hits the epitope site + the synonymous site; synonymous adds 0 peptides
    assert b.loc["S_mid", "n_recoding_sites_edited"] == 2
    assert b.loc["S_mid", "editing_neoantigen_burden"] >= 1

    # headline burden EQUALS count_binders() on the pooled altered peptides
    # (proves we route through the shared engine identically)
    peps = re_mod.recoding_peptides(make_recoding_site())["peptides"]
    expected = count_binders(peps, ALLELES, rank_threshold=WEAK_BINDER_RANK)
    assert int(b.loc["S_hi", "editing_neoantigen_burden"]) == expected, \
        (b.loc["S_hi", "editing_neoantigen_burden"], expected)
    print(f"  editing_neoantigen_burden == shared count_binders() = {expected}: OK")


def test_missing_hla_is_na():
    catalog = [make_recoding_site()]
    sites = make_editing_sites()
    burden = re_mod.compute_editing_neoantigen_burden(sites, catalog, hla_by_sample={})
    assert burden["editing_neoantigen_burden"].isna().all(), \
        "no HLA -> NA burden (never fabricate)"
    print("  missing HLA -> <NA> burden (graceful degradation): OK")


def test_build_editing_features_end_to_end():
    feat = re_mod.build_editing_features(
        aei_data=make_aei_table(),
        editing_sites=make_editing_sites(),
        recoding_catalog=[make_recoding_site(), make_synonymous_site()],
        hla_table=make_hla_table(),
        cohort="synth",
    )
    print(feat[["run_accession", "cohort", "alu_editing_index",
                "editing_neoantigen_burden", "n_recoding_sites_edited"]].to_string(index=False))
    # both named feature columns present and keyed on the contract key
    assert {"run_accession", "cohort", "alu_editing_index",
            "editing_neoantigen_burden"}.issubset(feat.columns)
    assert feat["alu_editing_index"].between(0, 1).all()
    assert len(feat) == 3
    print("  build_editing_features end-to-end (both named columns): OK")


def test_catalog_roundtrip(tmp="/tmp/_recoding_catalog.tsv"):
    cat = [make_recoding_site(), make_synonymous_site()]
    re_mod.save_recoding_catalog(cat, tmp)
    back = re_mod.load_recoding_catalog(tmp)
    assert len(back) == 2
    assert back[0].cds_window == cat[0].cds_window
    assert back[0].edit_offset == cat[0].edit_offset
    print("  recoding catalog save/load round-trip: OK")


if __name__ == "__main__":
    print("RNA-editing module validation:")
    test_sequence_utils()
    test_alu_editing_index()
    test_recoding_peptides_create_epitope()
    test_editing_neoantigen_burden()
    test_missing_hla_is_na()
    test_build_editing_features_end_to_end()
    test_catalog_roundtrip()
    print("\nALL RNA-EDITING MODULE TESTS PASSED")
