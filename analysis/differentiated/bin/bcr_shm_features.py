#!/usr/bin/env python
"""
bin/bcr_shm_features.py — per-sample BCR/SHM feature-row driver, called by the
BCR_SHM_FEATURES Nextflow process.

Reads ONE sample's TRUST4 output directory (containing <run>_cdr3.out and
optionally <run>_airr.tsv) and writes a single contract-shaped feature-row CSV
via analysis/differentiated/bcr_shm.py (build_bcr_features). Kept deliberately
thin — all repertoire logic lives in the unit-tested module.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parents[1]))   # analysis/differentiated
sys.path.insert(0, str(_HERE.parents[2]))   # analysis

import pandas as pd  # noqa: E402

import bcr_shm as bs  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--trust4-dir", required=True,
                    help="directory holding <run>_cdr3.out (TRUST4 --od/<run>)")
    ap.add_argument("--run-accession", required=True)
    ap.add_argument("--cohort", default=None)
    ap.add_argument("--min-clonotypes", type=int, default=3)
    ap.add_argument("--no-weight", action="store_true",
                    help="unweighted SHM mean for the cdr3 source (default: read-weighted)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    # trust4-dir may be the <run> dir itself or its parent; normalise so
    # build_bcr_features sees <root>/<run>/<run>_cdr3.out.
    tdir = Path(a.trust4_dir)
    if (tdir / f"{a.run_accession}_cdr3.out").exists():
        root = tdir.parent
    else:
        root = tdir
    idx = pd.DataFrame({"run_accession": [a.run_accession], "cohort": [a.cohort]})
    out = bs.build_bcr_features(str(root), idx,
                                weight=not a.no_weight,
                                min_clonotypes=a.min_clonotypes)
    out.to_csv(a.out, index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
