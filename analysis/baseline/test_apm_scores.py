#!/usr/bin/env python
"""
test_apm_scores.py — unit validation for the APM / HLA-expression score module.

Checks numerical properties on small controlled gene x sample matrices:
  * all curated APM genes resolve in a symbol-indexed matrix (no silent drop);
  * a sample with high APM expression scores above a low-APM sample;
  * within-cohort z-score harmonisation collapses a synthetic per-cohort
    additive offset (the platform-robustness guarantee);
  * b2m_hla_abc is driven specifically by HLA-A/B/C + B2M (knocking those down
    drops the floor score while leaving class-II unchanged);
  * exclude_gep_overlap actually removes the GEP-shared gene from scoring;
  * score_all returns the (run_accession, cohort)-keyed contract frame.

Run:  python analysis/baseline/test_apm_scores.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from analysis.baseline import apm_scores as apm  # noqa: E402
from analysis.baseline.gep_scores import AYERS_HOUSEKEEPING_GENES  # noqa: E402


ALL_GENES = tuple(
    dict.fromkeys(
        apm.APM_CLASS1_GENES + apm.APM_CLASS2_GENES + apm.B2M_HLA_ABC_GENES
        + AYERS_HOUSEKEEPING_GENES
    )
)


def _matrix(sample_levels, seed=0, hk_level=5.0):
    """Build a genes x samples log-normal-ish TPM matrix.

    sample_levels: dict sample_id -> dict(class1=, class2=, hk=) mean log2 level
    applied to those gene blocks; everything else gets hk_level. Returns a
    linear-TPM frame (genes x samples).
    """
    rng = np.random.default_rng(seed)
    genes = list(ALL_GENES)
    cols = {}
    for sid, lv in sample_levels.items():
        base = np.full(len(genes), lv.get("other", hk_level), dtype=float)
        for i, g in enumerate(genes):
            if g in apm.APM_CLASS1_GENES:
                base[i] = lv.get("class1", hk_level)
            if g in apm.APM_CLASS2_GENES:
                base[i] = lv.get("class2", hk_level)
            if g in AYERS_HOUSEKEEPING_GENES:
                base[i] = lv.get("hk", hk_level)
        noise = rng.normal(0, 0.05, size=len(genes))
        cols[sid] = np.power(2.0, base + noise) - 1.0
    return pd.DataFrame(cols, index=genes).clip(lower=0)


def test_all_genes_resolve():
    mat = _matrix({"s1": {}, "s2": {}})
    # every curated gene present -> _resolve_genes must not raise and must keep all
    from analysis.baseline.gep_scores import _resolve_genes, _to_log_tpm
    log = _to_log_tpm(mat)
    for name, genes in [("class1", apm.APM_CLASS1_GENES),
                        ("class2", apm.APM_CLASS2_GENES),
                        ("floor", apm.B2M_HLA_ABC_GENES)]:
        present = _resolve_genes(log.index, genes, name)
        assert len(present) == len(genes), f"{name}: dropped {set(genes)-set(present)}"
    print("  all curated APM genes resolve:", len(ALL_GENES), "union genes")


def test_high_vs_low_apm():
    # s_hi high class-I, s_lo low class-I, identical HK -> hi scores higher (raw)
    mat = _matrix({"hi": {"class1": 9.0}, "lo": {"class1": 2.0}})
    s = apm.score_apm_class1(mat, harmonize="none")
    assert s["hi"] > s["lo"], f"high-APM sample not higher: {s.to_dict()}"
    print(f"  high vs low class-I (raw): hi={s['hi']:.3f} > lo={s['lo']:.3f}")


def test_harmonization_collapses_offset():
    # two cohorts, identical biology but +3 log2 additive offset in cohort B;
    # within-batch z-score must remove the cohort mean difference.
    rng = np.random.default_rng(1)
    levels = {}
    for i in range(20):
        base = rng.normal(6.0, 1.0)
        levels[f"A{i}"] = {"class1": base, "class2": base, "hk": 5.0}
    for i in range(20):
        base = rng.normal(6.0, 1.0)
        levels[f"B{i}"] = {"class1": base + 3.0, "class2": base + 3.0, "hk": 5.0 + 3.0}
    mat = _matrix(levels, seed=2)
    batches = pd.Series({s: ("A" if s.startswith("A") else "B") for s in mat.columns})
    s = apm.score_apm_class1(mat, batches, harmonize="zscore")
    mean_A = s[[c for c in s.index if c.startswith("A")]].mean()
    mean_B = s[[c for c in s.index if c.startswith("B")]].mean()
    assert abs(mean_A) < 1e-9 and abs(mean_B) < 1e-9, f"z-mean not ~0: A={mean_A}, B={mean_B}"
    print(f"  harmonised per-cohort means ~0: A={mean_A:.2e}, B={mean_B:.2e}")


def test_floor_is_hla_abc_specific():
    # knock down HLA-A/B/C + B2M only -> b2m_hla_abc drops, class-II unchanged.
    hi = _matrix({"ref": {}}, seed=3)["ref"]
    lo = hi.copy()
    for g in apm.B2M_HLA_ABC_GENES:
        lo[g] = 2 ** 1.0 - 1.0  # very low
    mat = pd.concat({"ref": hi, "kd": lo}, axis=1)
    floor = apm.score_b2m_hla_abc(mat, harmonize="none")
    c2 = apm.score_apm_class2(mat, harmonize="none")
    assert floor["ref"] > floor["kd"], "floor did not drop on HLA-ABC/B2M knockdown"
    assert abs(c2["ref"] - c2["kd"]) < 0.2, "class-II moved on class-I knockdown (should be ~flat)"
    print(f"  floor ref={floor['ref']:.2f} > kd={floor['kd']:.2f}; class-II ~flat "
          f"(Δ={abs(c2['ref']-c2['kd']):.3f})")


def test_exclude_gep_overlap():
    mat = _matrix({"s1": {}, "s2": {}})
    from analysis.baseline.gep_scores import _prep, _resolve_genes
    log = _prep(mat, AYERS_HOUSEKEEPING_GENES)
    full = _resolve_genes(log.index, apm.APM_CLASS1_GENES, "full")
    genes_excl = tuple(g for g in apm.APM_CLASS1_GENES if g not in set(apm.APM_CLASS1_GEP_OVERLAP))
    excl = _resolve_genes(log.index, genes_excl, "excl")
    assert set(full) - set(excl) == set(apm.APM_CLASS1_GEP_OVERLAP)
    # and the scorer runs with the flag
    s = apm.score_apm_class1(mat, harmonize="none", exclude_gep_overlap=True)
    assert s.notna().all()
    print(f"  exclude_gep_overlap drops {apm.APM_CLASS1_GEP_OVERLAP}: "
          f"{len(full)} -> {len(excl)} genes")


def test_score_all_contract():
    mat = _matrix({"r1": {}, "r2": {}, "r3": {}}, seed=4)
    sm = pd.DataFrame(
        {"run_accession": ["r1", "r2", "r3"], "cohort": ["A", "A", "B"]},
        index=["r1", "r2", "r3"],
    )
    out = apm.score_all(mat, sm)
    assert list(out.columns) == ["run_accession", "cohort", *apm.FEATURE_COLUMNS]
    assert len(out) == 3 and out[list(apm.FEATURE_COLUMNS)].notna().all().all()
    print("  score_all contract OK:", list(out.columns))


if __name__ == "__main__":
    print("APM / HLA-expression score module validation:")
    test_all_genes_resolve()
    test_high_vs_low_apm()
    test_harmonization_collapses_offset()
    test_floor_is_hla_abc_specific()
    test_exclude_gep_overlap()
    test_score_all_contract()
    print("\nALL APM SCORE TESTS PASSED")
