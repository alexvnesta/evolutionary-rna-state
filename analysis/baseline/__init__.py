"""Baseline (proven-floor) IO-response features — reproduce, don't innovate.

FDA-approved and research-grade features every ICB predictor is measured against:

- ``tmb_standardized`` — batch-harmonized tumor mutational burden (FoCR/Merino).
- ``snv_indel_neoantigen`` — SNV/indel neoantigen load via the shared MHC engine.
- ``gep_scores`` — Ayers 18-gene T-cell-inflamed GEP, 6-gene IFN-gamma, and the
  Mariathasan Teff/TGF-beta balance, with within-cohort harmonization.
"""
