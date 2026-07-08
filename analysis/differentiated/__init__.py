"""
analysis/differentiated — the DIFFERENTIATED antigen modules.

These are the RNA-state-specific neoantigen modules that a WES/annotation
pipeline cannot see (splicing, TE/ERV, intron retention, fusion, editing).
Each derives candidate peptides from a de-novo phenotype matrix, then scores
them through the SHARED antigen core (analysis/antigen_core) so every
``*_neoantigen_burden`` is defined identically and comparably.

Public surface
--------------
    from analysis.differentiated.splicing_neoantigen import (
        splice_neoantigen_burden,          # per-sample burden int (the feature)
        peptides_from_neojunctions,        # neojunction -> 8-11mer peptides
        call_neojunctions,                 # SNAF tumor-specificity gate
        translate_junction,                # SNAF 3-frame junction translation
        parse_star_sj,                     # STAR SJ.out.tab adapter
        Junction,
    )
"""
