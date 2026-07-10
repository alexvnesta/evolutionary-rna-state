#!/usr/bin/env python
"""
test_bcr_shm.py — unit validation for the TRUST4 BCR/SHM feature module.

Writes SYNTHETIC TRUST4 _cdr3.out files with known germline similarities,
isotypes and read counts into a temp dir, then checks that:
  * SHM rate = 1 - mean germline identity (read-weighted), from cdr3 source;
  * AIRR v_identity is preferred over cdr3 similarity when present;
  * IgG / switched fractions match the isotype read composition;
  * clonality is 0 for perfectly even, ->1 for dominated repertoires;
  * a too-shallow sample (< min_clonotypes IGH contigs) returns NaN features,
    never imputed;
  * build_bcr_features emits the (run_accession, cohort)-keyed contract frame.

Run:  python analysis/differentiated/test_bcr_shm.py
"""
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from analysis.differentiated import bcr_shm as bs  # noqa: E402


def _write_cdr3(d, run, rows):
    """rows: list of (V,C,similarity,readcount). Writes <run>/<run>_cdr3.out."""
    os.makedirs(os.path.join(d, run), exist_ok=True)
    lines = []
    for i, (V, C, sim, rc) in enumerate(rows):
        # cols: cid idx V D J C CDR1 CDR2 CDR3 CDR3_score readcount germ_sim complete
        lines.append("\t".join([
            f"c{i}", "0", V, "*", "IGHJ4", C, "AAA", "BBB", "CARWYFDVW",
            "1.00", str(rc), ("" if sim is None else f"{sim:.4f}"), "1",
        ]))
    with open(os.path.join(d, run, f"{run}_cdr3.out"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


def test_shm_rate_cdr3_source():
    with tempfile.TemporaryDirectory() as d:
        # three IGH clonotypes, germline sim 0.90/0.80/0.85 (rc 100 each)
        # read-weighted SHM = 1 - (0.90+0.80+0.85)/3 = 1 - 0.85 = 0.15
        _write_cdr3(d, "S1", [("IGHV1-2", "IGHG1", 0.90, 100),
                              ("IGHV3-7", "IGHM", 0.80, 100),
                              ("IGHV4-1", "IGHA1", 0.85, 100)])
        f = bs.sample_bcr_features(os.path.join(d, "S1"), "S1")
        assert abs(f["bcr_shm_rate"] - 0.15) < 1e-9, f["bcr_shm_rate"]
        assert f["_shm_source"] == "cdr3_germline_similarity"
        print(f"  SHM (cdr3, read-weighted) = {f['bcr_shm_rate']:.3f} (expect 0.150)")


def test_airr_preferred():
    with tempfile.TemporaryDirectory() as d:
        _write_cdr3(d, "S2", [("IGHV1-2", "IGHG1", 0.90, 100),
                              ("IGHV3-7", "IGHM", 0.90, 100),
                              ("IGHV4-1", "IGHA1", 0.90, 100)])
        # AIRR says v_identity 0.95 for all -> SHM 0.05, must override cdr3's 0.10
        airr = pd.DataFrame({"v_identity": ["0.95", "0.95", "0.95"]})
        airr.to_csv(os.path.join(d, "S2", "S2_airr.tsv"), sep="\t", index=False)
        f = bs.sample_bcr_features(os.path.join(d, "S2"), "S2")
        assert f["_shm_source"] == "airr_v_identity", f["_shm_source"]
        assert abs(f["bcr_shm_rate"] - 0.05) < 1e-9, f["bcr_shm_rate"]
        print(f"  SHM (AIRR preferred) = {f['bcr_shm_rate']:.3f} (expect 0.050)")


def test_isotype_fractions():
    with tempfile.TemporaryDirectory() as d:
        # reads: IGHG 60, IGHA 20, IGHM 20 -> IgG 0.6, switched (G+A) 0.8
        _write_cdr3(d, "S3", [("IGHV1-2", "IGHG1", 0.9, 60),
                              ("IGHV3-7", "IGHA1", 0.9, 20),
                              ("IGHV4-1", "IGHM", 0.9, 20)])
        f = bs.sample_bcr_features(os.path.join(d, "S3"), "S3")
        assert abs(f["bcr_igg_fraction"] - 0.6) < 1e-9, f["bcr_igg_fraction"]
        assert abs(f["bcr_switched_fraction"] - 0.8) < 1e-9, f["bcr_switched_fraction"]
        print(f"  IgG frac = {f['bcr_igg_fraction']:.2f} (0.60), switched = {f['bcr_switched_fraction']:.2f} (0.80)")


def test_clonality():
    with tempfile.TemporaryDirectory() as d:
        # even 4-way -> clonality ~0
        _write_cdr3(d, "EV", [("IGHV1-2", "IGHM", 0.9, 25)] * 1 +
                              [("IGHV3-7", "IGHM", 0.9, 25),
                               ("IGHV4-1", "IGHM", 0.9, 25),
                               ("IGHV2-5", "IGHM", 0.9, 25)])
        fe = bs.sample_bcr_features(os.path.join(d, "EV"), "EV")
        # dominated: one clone 970, three 10 -> clonality high
        _write_cdr3(d, "DM", [("IGHV1-2", "IGHG1", 0.9, 970),
                              ("IGHV3-7", "IGHM", 0.9, 10),
                              ("IGHV4-1", "IGHM", 0.9, 10),
                              ("IGHV2-5", "IGHM", 0.9, 10)])
        fd = bs.sample_bcr_features(os.path.join(d, "DM"), "DM")
        assert fe["bcr_clonality"] < 0.02, fe["bcr_clonality"]
        assert fd["bcr_clonality"] > 0.6, fd["bcr_clonality"]
        print(f"  clonality even={fe['bcr_clonality']:.3f} (~0), dominated={fd['bcr_clonality']:.3f} (>0.6)")


def test_shallow_returns_nan():
    with tempfile.TemporaryDirectory() as d:
        # only 2 IGH contigs, min_clonotypes=3 -> features NaN, counts present
        _write_cdr3(d, "SH", [("IGHV1-2", "IGHG1", 0.9, 5),
                              ("IGHV3-7", "IGHM", 0.8, 5)])
        f = bs.sample_bcr_features(os.path.join(d, "SH"), "SH", min_clonotypes=3)
        assert np.isnan(f["bcr_shm_rate"]) and np.isnan(f["bcr_clonality"])
        assert f["bcr_n_clonotypes"] == 2.0 and f["bcr_n_reads"] == 10.0
        print(f"  shallow sample: SHM={f['bcr_shm_rate']} (NaN), n_clono={f['bcr_n_clonotypes']:.0f} reported")


def test_build_contract():
    with tempfile.TemporaryDirectory() as d:
        _write_cdr3(d, "R1", [("IGHV1-2", "IGHG1", 0.9, 50),
                              ("IGHV3-7", "IGHM", 0.85, 50),
                              ("IGHV4-1", "IGHA1", 0.88, 50)])
        # R2 has no TRUST4 output -> all-NaN row, not dropped
        idx = pd.DataFrame({"run_accession": ["R1", "R2"], "cohort": ["A", "B"]})
        out = bs.build_bcr_features(d, idx)
        assert list(out.columns)[:2] == ["run_accession", "cohort"]
        assert set(bs.FEATURE_COLUMNS).issubset(out.columns)
        assert len(out) == 2
        assert not np.isnan(out.loc[0, "bcr_shm_rate"])
        assert np.isnan(out.loc[1, "bcr_shm_rate"])  # missing sample -> NaN
        print("  build contract OK; missing sample -> NaN row preserved")


if __name__ == "__main__":
    print("TRUST4 BCR/SHM feature module validation:")
    test_shm_rate_cdr3_source()
    test_airr_preferred()
    test_isotype_fractions()
    test_clonality()
    test_shallow_returns_nan()
    test_build_contract()
    print("\nALL BCR/SHM TESTS PASSED")
