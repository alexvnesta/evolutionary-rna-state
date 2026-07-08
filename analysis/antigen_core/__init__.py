"""
antigen_core — the shared antigen-derivation keystone.

Every downstream antigen module (splicing, TE/ERV, intron-retention, fusion,
SNV/indel) imports the SAME peptide-scoring engine and HLA typing from here, so
the "*_neoantigen_burden" features are all defined identically and comparably.

Public surface
--------------
    from analysis.antigen_core.mhc_binding import (
        score_peptides, count_binders, binder_counts, best_per_peptide,
        STRONG_BINDER_RANK, WEAK_BINDER_RANK,
    )
    from analysis.antigen_core.hla_typing import (
        type_sample, summarize_genotype, build_hla_table,
        parse_arcashla_genotype, parse_optitype_result,
    )
"""
from . import mhc_binding, hla_typing  # noqa: F401

__all__ = ["mhc_binding", "hla_typing"]
