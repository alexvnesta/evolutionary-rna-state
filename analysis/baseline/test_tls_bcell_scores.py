#!/usr/bin/env python
"""
test_tls_bcell_scores.py — unit validation for the TLS / B-cell score module.

Checks on small controlled gene x sample matrices:
  * all B-cell-lineage and TLS-chemokine genes resolve (no silent drop);
  * a B-cell-high sample scores above a B-cell-low sample (raw);
  * within-cohort z-score harmonisation collapses a per-cohort additive offset;
  * tls_imprint is the mean of the two standardised arms;
  * score_all returns the (run_accession, cohort)-keyed contract frame.

Run:  python analysis/baseline/test_tls_bcell_scores.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from analysis.baseline import tls_bcell_scores as tls  # noqa: E402
from analysis.baseline.gep_scores import AYERS_HOUSEKEEPING_GENES  # noqa: E402


ALL_GENES = tuple(
    dict.fromkeys(
        tls.BCELL_LINEAGE_GENES + tls.TLS_CHEMOKINE_GENES + AYERS_HOUSEKEEPING_GENES
    )
)


def _matrix(sample_levels, seed=0, hk_level=5.0):
    rng = np.random.default_rng(seed)
    genes = list(ALL_GENES)
    cols = {}
    for sid, lv in sample_levels.items():
        base = np.full(len(genes), lv.get("other", hk_level), dtype=float)
        for i, g in enumerate(genes):
            if g in tls.BCELL_LINEAGE_GENES:
                base[i] = lv.get("bcell", hk_level)
            if g in tls.TLS_CHEMOKINE_GENES:
                base[i] = lv.get("chemokine", hk_level)
            if g in AYERS_HOUSEKEEPING_GENES:
                base[i] = lv.get("hk", hk_level)
        noise = rng.normal(0, 0.05, size=len(genes))
        cols[sid] = np.power(2.0, base + noise) - 1.0
    return pd.DataFrame(cols, index=genes).clip(lower=0)


def test_all_genes_resolve():
    from analysis.baseline.gep_scores import _resolve_genes, _to_log_tpm
    log = _to_log_tpm(_matrix({"s1": {}, "s2": {}}))
    for name, genes in [("bcell", tls.BCELL_LINEAGE_GENES),
                        ("chemokine", tls.TLS_CHEMOKINE_GENES)]:
        present = _resolve_genes(log.index, genes, name)
        assert len(present) == len(genes), f"{name}: dropped {set(genes)-set(present)}"
    print("  all B-cell + TLS-chemokine genes resolve:", len(ALL_GENES), "union genes")


def test_high_vs_low_bcell():
    mat = _matrix({"hi": {"bcell": 9.0}, "lo": {"bcell": 2.0}})
    s = tls.score_bcell_lineage(mat, harmonize="none")
    assert s["hi"] > s["lo"], f"B-cell-high sample not higher: {s.to_dict()}"
    print(f"  high vs low B-cell (raw): hi={s['hi']:.3f} > lo={s['lo']:.3f}")


def test_harmonization_collapses_offset():
    rng = np.random.default_rng(1)
    levels = {}
    for i in range(20):
        base = rng.normal(6.0, 1.0)
        levels[f"A{i}"] = {"bcell": base, "chemokine": base, "hk": 5.0}
    for i in range(20):
        base = rng.normal(6.0, 1.0)
        levels[f"B{i}"] = {"bcell": base + 3.0, "chemokine": base + 3.0, "hk": 8.0}
    mat = _matrix(levels, seed=2)
    batches = pd.Series({s: ("A" if s.startswith("A") else "B") for s in mat.columns})
    s = tls.score_bcell_lineage(mat, batches, harmonize="zscore")
    mean_A = s[[c for c in s.index if c.startswith("A")]].mean()
    mean_B = s[[c for c in s.index if c.startswith("B")]].mean()
    assert abs(mean_A) < 1e-9 and abs(mean_B) < 1e-9, f"z-mean not ~0: A={mean_A}, B={mean_B}"
    print(f"  harmonised per-cohort means ~0: A={mean_A:.2e}, B={mean_B:.2e}")


def test_imprint_is_mean_of_arms():
    mat = _matrix({f"s{i}": {"bcell": 4.0 + i, "chemokine": 3.0 + 0.5 * i} for i in range(6)},
                  seed=5)
    batches = pd.Series({c: "A" for c in mat.columns})
    b = tls.score_bcell_lineage(mat, batches, harmonize="zscore")
    c = tls.score_tls_chemokine(mat, batches, harmonize="zscore")
    imp = tls.score_tls_imprint(mat, batches, harmonize="zscore")
    assert np.allclose(imp.values, ((b + c) / 2.0).values), "imprint != mean of arms"
    print(f"  tls_imprint == mean(bcell_z, chemokine_z): max|Δ|="
          f"{np.max(np.abs(imp.values - ((b+c)/2).values)):.2e}")


def test_score_all_contract():
    mat = _matrix({"r1": {}, "r2": {}, "r3": {}}, seed=4)
    sm = pd.DataFrame(
        {"run_accession": ["r1", "r2", "r3"], "cohort": ["A", "A", "B"]},
        index=["r1", "r2", "r3"],
    )
    out = tls.score_all(mat, sm)
    assert list(out.columns) == ["run_accession", "cohort", *tls.FEATURE_COLUMNS]
    assert len(out) == 3 and out[list(tls.FEATURE_COLUMNS)].notna().all().all()
    print("  score_all contract OK:", list(out.columns))


if __name__ == "__main__":
    print("TLS / B-cell score module validation:")
    test_all_genes_resolve()
    test_high_vs_low_bcell()
    test_harmonization_collapses_offset()
    test_imprint_is_mean_of_arms()
    test_score_all_contract()
    print("\nALL TLS/B-CELL SCORE TESTS PASSED")
