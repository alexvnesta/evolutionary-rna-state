# Clonal-trajectory analysis code

Pipeline for the clonal-evolutionary-state investigation. Each cohort follows the same
stages: build AnnData -> infercnvpy clone inference -> non-circular validation -> clone
pseudobulk + program scoring -> harmonized immune-hot label. Results feed the pooled
classifier and the real-bulk transfer / ICB tests. Run in the `scrna-cnv` conda env
(python 3.11, infercnvpy 0.6.1, scanpy 1.11.5); set `NUMBA_CACHE_DIR` before first import.

## Layout

Shared / cSCC (phase 0):
- `load_adata.py` — AnnData loader / normalization helpers.
- `run_infercnv.py` — infercnvpy clone inference (serial process_map monkeypatch, cnv_pca).
- `phase0_validate.py` — non-circular validation (silhouette vs within-patient permutation null).
- `phase0_effectsize.py` — within- vs between-patient subclone separation.
- `phase0_programs.py` — hallmark program scoring + k=2 recurrent classes.
- `phase0_dilution.py` — purity-gradient separability (stromal/immune admixture robustness).
- `phase0_compartment.py` — malignant-compartment-specific feature recovery at low purity.

DDLPS (cross-cohort replication + real matched bulk):
- `build_ddlps.py`, `run_ddlps.py` — build + clone pipeline.
- `ddlps_pseudobulk.py` — clone pseudobulk for the pooled classifier.
- `ddlps_compartment.py` — compartment-specific features scored on real DDLPS bulk.

AML (`aml/`) — van Galen GSE116256, genomically quiet:
- `build_adata.py`, `build_ref.py`, `run_infercnv.py`/`run_infercnv2.py`, `run_clones.py`,
  `validate.py`, `programs.py`, `characterize.py`, `coupling_test.py`, `add_strict_gate.py`.

LUAD (`luad/`) — Kim GSE131907, largest single cohort:
- `build_luad.py` — streaming sparse reader for the 208k-cell text matrix.
- `run_luad.py` — full clone pipeline.

Classifier and endpoint:
- `pooled_classifier.py` — pooled L2-logistic classifier, grouped-by-patient CV +
  leave-one-cohort-out + within-cohort permutation null.
- `icb/riaz_icb_test.py` — Riaz melanoma bulk ICB outcome association test.

## Reproducing

Per-cohort pseudobulk (`../tables/*_clone_pseudobulk.npy` + meta + genes) is committed so the
pooled classifier reproduces without re-running inferCNV. Full method: `../reports/methods.md`.
