#!/usr/bin/env python
"""Tests for pilot_crosswalk + pilot_gep — the pipeline->frame bridge.

Run: MHCFLURRY_DATA_DIR=reference/mhcflurry_models python analysis/test_pilot_crosswalk.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pandas as pd
from analysis.pilot_crosswalk import (
    _gide_alias_to_sampleid, _riaz_title_to_patientid, build_crosswalk, attach_labels)


def test_gide_alias_rule():
    assert _gide_alias_to_sampleid("PD1_13_PRE") == "PD13_Pre"
    assert _gide_alias_to_sampleid("PD1_5_PRE") == "PD05_Pre"      # zero-pad
    assert _gide_alias_to_sampleid("ipiPD1_24_PRE") == "iP24_Pre"  # combo arm
    assert _gide_alias_to_sampleid("ipiPD1_2_PRE") == "iP02_Pre"
    assert _gide_alias_to_sampleid("garbage") is None
    print("PASS test_gide_alias_rule")


def test_riaz_title_rule():
    assert _riaz_title_to_patientid("Pt72_Pre_AD793922-5") == "Pt72"
    assert _riaz_title_to_patientid("Pt18_Pre_E9024732-6") == "Pt18"
    assert _riaz_title_to_patientid("nope") is None
    print("PASS test_riaz_title_rule")


def test_build_and_attach_synthetic():
    # minimal synthetic ENA + frame; verify mapping + label attach + no-guess drop
    ena = pd.DataFrame([
        dict(run_accession="ERR1", project="PRJEB23709", sample_alias="PD1_13_PRE", sample_title="PD1_13_PRE"),
        dict(run_accession="ERR2", project="PRJEB23709", sample_alias="ipiPD1_24_PRE", sample_title="ipiPD1_24_PRE"),
        dict(run_accession="ERR9", project="PRJEB23709", sample_alias="PD1_99_PRE", sample_title="PD1_99_PRE"),  # not in frame
        dict(run_accession="SRR1", project="PRJNA356761", sample_alias="GSM1", sample_title="Pt72_Pre_AD793922-5"),
    ])
    af = pd.DataFrame([
        dict(sampleId="PD13_Pre", cohort="gide2019", patientId="PD13", RESPONDER=True,  TIDE_RESPONDER=1, TMB_NONSYNONYMOUS=5.0),
        dict(sampleId="iP24_Pre", cohort="gide2019", patientId="iP24", RESPONDER=False, TIDE_RESPONDER=0, TMB_NONSYNONYMOUS=3.0),
        dict(sampleId="Pt72_pre", cohort="riaz2017", patientId="Pt72", RESPONDER=True,  TIDE_RESPONDER=1, TMB_NONSYNONYMOUS=8.0),
    ])
    xw = build_crosswalk(ena, af)
    # ERR9 has no frame row -> dropped, never guessed
    assert set(xw["run_accession"]) == {"ERR1", "ERR2", "SRR1"}, xw
    assert xw.set_index("run_accession").loc["ERR1", "sampleId"] == "PD13_Pre"
    assert xw.set_index("run_accession").loc["SRR1", "sampleId"] == "Pt72_pre"

    fm = pd.DataFrame({"run_accession": ["ERR1", "ERR2", "SRR1"],
                       "cohort": ["gide2019", "gide2019", "riaz2017"],
                       "some_feature": [1.0, 2.0, 3.0]})
    ev = attach_labels(fm, xw, af)
    assert ev["RESPONDER"].notna().all()
    assert bool(ev.set_index("run_accession").loc["ERR1", "RESPONDER"]) is True
    print("PASS test_build_and_attach_synthetic")


def test_gep_orientation_normalization():
    from analysis.pilot_gep import to_symbol_gene_matrix, signature_sym2ens
    # tiny samples x ENSG frame with two known signature genes
    s2e = signature_sym2ens()
    cd8a = s2e.get("CD8A"); stat1 = s2e.get("STAT1")
    assert cd8a and stat1
    raw = pd.DataFrame({
        "run_accession": ["r1", "r2"], "cohort": ["c", "c"],
        cd8a: [10.0, 1.0], stat1: [5.0, 0.5],
    })
    mat, sm = to_symbol_gene_matrix(raw)
    assert set(["CD8A", "STAT1"]).issubset(set(mat.index)), mat.index.tolist()
    assert list(mat.columns) == ["r1", "r2"]
    assert list(sm["run_accession"]) == ["r1", "r2"]
    print("PASS test_gep_orientation_normalization")


def test_crosswalk_from_catalog_3cohort():
    """run_catalog path covers all 3 cohorts incl. hugo; gide via alias rule,
    hugo/riaz via patientId; PRE-only; cohort-less feature matrix still labels."""
    from analysis.pilot_crosswalk import build_crosswalk_from_catalog, attach_labels
    cat = pd.DataFrame([
        dict(run_accession="ERRg", cohort="gide2019", patient_id="PD1_13", timepoint="PRE"),
        dict(run_accession="ERRgi", cohort="gide2019", patient_id="ipiPD1_24", timepoint="PRE"),
        dict(run_accession="SRRh", cohort="hugo2016", patient_id="Pt2", timepoint="PRE"),
        dict(run_accession="SRRr", cohort="riaz2017", patient_id="Pt72", timepoint="PRE"),
        dict(run_accession="SRRon", cohort="riaz2017", patient_id="Pt72", timepoint="ON"),  # dropped
    ])
    af = pd.DataFrame([
        dict(sampleId="PD13_Pre", cohort="gide2019", patientId="PD13", RESPONDER=True,  TIDE_RESPONDER=1, TMB_NONSYNONYMOUS=5.0),
        dict(sampleId="iP24_Pre", cohort="gide2019", patientId="iP24", RESPONDER=False, TIDE_RESPONDER=0, TMB_NONSYNONYMOUS=3.0),
        dict(sampleId="Pt2",      cohort="hugo2016", patientId="Pt2",  RESPONDER=True,  TIDE_RESPONDER=1, TMB_NONSYNONYMOUS=9.0),
        dict(sampleId="Pt72_pre", cohort="riaz2017", patientId="Pt72", RESPONDER=True,  TIDE_RESPONDER=1, TMB_NONSYNONYMOUS=8.0),
    ])
    xw = build_crosswalk_from_catalog(cat, af)
    assert set(xw["run_accession"]) == {"ERRg", "ERRgi", "SRRh", "SRRr"}, xw  # ON dropped
    assert set(xw["cohort"]) == {"gide2019", "hugo2016", "riaz2017"}
    # feature matrix with a NULL cohort column must still get labels (join on run_accession)
    fm = pd.DataFrame({"run_accession": ["ERRg", "SRRh", "SRRr"],
                       "cohort": [None, None, None], "feat": [1.0, 2.0, 3.0]})
    ev = attach_labels(fm, xw, af)
    assert ev["RESPONDER"].notna().all(), ev
    assert ev["cohort"].notna().all()  # cohort filled from crosswalk, not the null input
    print("PASS test_crosswalk_from_catalog_3cohort")


def test_regulator_activity_vendored():
    """regulator_activity module scores 3 sets via the vendored fallback."""
    from analysis.regulator_activity import build_regulator_activity, regulator_sets
    sets = regulator_sets()
    assert set(sets) == {"SPLICING_FACTOR", "RBP_BROAD", "ADAR_EDITING"}
    # 6 samples, plant ADAR set as high-variance; map symbols to fake ensembl
    import numpy as np
    syms = sorted({s for v in sets.values() for s in v})
    gmap = {s: f"ENSGFAKE{i}" for i, s in enumerate(syms)}
    rng = np.random.default_rng(0)
    tpm = pd.DataFrame(rng.lognormal(2, 1, size=(6, len(syms))),
                       index=[f"r{i}" for i in range(6)], columns=[gmap[s] for s in syms])
    S = build_regulator_activity(tpm, gene_symbol_index=gmap)
    assert "run_accession" in S.columns
    assert {"SPLICING_FACTOR", "RBP_BROAD", "ADAR_EDITING"}.issubset(S.columns)
    assert len(S) == 6
    print("PASS test_regulator_activity_vendored")


if __name__ == "__main__":
    test_gide_alias_rule()
    test_riaz_title_rule()
    test_build_and_attach_synthetic()
    test_gep_orientation_normalization()
    test_crosswalk_from_catalog_3cohort()
    test_regulator_activity_vendored()
    print("\nALL PILOT-CROSSWALK TESTS PASSED")
