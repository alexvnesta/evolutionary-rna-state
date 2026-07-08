"""Analysis package for the evolutionary-rna-state ICB project.

Per-sample, named, interpretable RNA feature modules for immune-checkpoint-blockade
response modeling, organized into subpackages:

- ``antigen_core`` — shared MHC-binding engine (MHCflurry) + arcasHLA HLA typing.
- ``baseline`` — proven-floor features (standardized TMB, SNV/indel neoantigen load,
  T-cell-inflamed GEP / IFN-gamma).
- ``differentiated`` — the project's differentiated features (splicing / TE / intron-
  retention / RNA-editing / fusion neoantigens).
- ``eval`` — batch-robustness and leave-one-cohort/patient-out evaluation vs the floor.

Top-level modules (``pilot_ingest``, ``pilot_crosswalk``, ``pilot_gep``,
``regulator_activity``, ``registry_update``) form the pilot ingest + orchestration
layer that turns a real pipeline TPM matrix into a scored, label-attached feature
matrix with a floor-referenced verdict.
"""
