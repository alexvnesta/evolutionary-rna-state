//
// Subworkflow: TE_ANTIGEN
// Transposable-element / ERV-derived antigen burden.
// Consumes the pipeline session's Telescope per-locus TE counts + the rnaseq
// spine's HLA typing, derives ORF peptides from EXPRESSED TE/ERV loci, and
// scores them through the SHARED antigen_core MHCflurry engine to produce the
// per-sample te_antigen_burden family (locus- and family-resolved).
//
// Emits: te_antigen.parquet — one row per sample (run_accession, cohort) with
//   te_antigen_burden, te_antigen_burden_strong,
//   te_antigen_burden_{LINE,SINE,LTR,ERV}, and QC columns.
//
// This is a documented, correct STUB the pipeline wires into. Peptide
// derivation + the burden logic live in analysis/differentiated/te_antigen.py
// and are unit-tested in test_te_antigen.py; the shared peptide-scoring engine
// is analysis/antigen_core/mhc_binding.py. The stub{} blocks let
// `nextflow -stub-run` exercise the wiring without Telescope / a genome / the
// MHCflurry models present.
//
// BATCH ROBUSTNESS (project hard constraint): TE quantification is
// MULTIMAP-sensitive. Every sample's Telescope counts MUST come from the SAME
// STAR multimap settings (--outFilterMultimapNmax / --winAnchorMultimapNmax /
// --outSAMmultNmax) and the same Telescope run across the cohort, else the
// per-locus counts are not comparable. Activity is called on a WITHIN-sample
// CPM floor and the burden is an allele-calibrated MHCflurry percentile-rank
// count, both composition/depth-invariant. See te_antigen.build_batch_robustness_note().
//
// Channel contract (matches the hla_typing / intron_retention subworkflows):
//   ch_te_locus   : path te_locus.parquet   (cohort-level, rows=run_accession)
//   ch_annotation : path TE locus annotation (RepeatMasker/GENCODE-TE)
//   ch_genome     : path genome FASTA (+ .fai)   [ meta-free, value channel ]
//   ch_hla        : path hla_typing.parquet  (from HLA_TYPING subworkflow)
//

include { BUILD_TE_ANTIGEN_TABLE } from './modules/build_te_antigen_table.nf'

workflow TE_ANTIGEN {
    take:
    ch_te_locus     // channel: path te_locus.parquet
    ch_annotation   // channel: path TE annotation (locus_id, repeat_class, chrom, start, end, strand)
    ch_genome       // channel: path genome.fa (+ .fai alongside)
    ch_hla          // channel: path hla_typing.parquet

    main:
    ch_versions = Channel.empty()

    // one-shot per cohort: expressed loci -> ORF peptides -> shared MHC engine
    // -> te_antigen.parquet. (Per-sample fan-out is an internal loop in the
    // build script; the matrices are already cohort-level tidy-wide.)
    BUILD_TE_ANTIGEN_TABLE(
        ch_te_locus,
        ch_annotation,
        ch_genome,
        ch_hla
    )
    ch_versions = ch_versions.mix( BUILD_TE_ANTIGEN_TABLE.out.versions )

    emit:
    te_table = BUILD_TE_ANTIGEN_TABLE.out.table    // te_antigen.parquet
    versions = ch_versions
}

// ---------------------------------------------------------------------------
// Inline module definition (kept in one file for a self-contained stub; in the
// pipeline this lives under ./modules/build_te_antigen_table.nf like the
// sibling subworkflows).
// ---------------------------------------------------------------------------

process BUILD_TE_ANTIGEN_TABLE {
    label 'process_medium'
    publishDir "${params.outdir}/features", mode: 'copy'
    // needs the antigen env: mhcflurry + biopython + pysam + pandas/pyarrow
    conda "bioconda::mhcflurry conda-forge::biopython bioconda::pysam conda-forge::pandas conda-forge::pyarrow"

    input:
    path te_locus       // te_locus.parquet
    path annotation     // TE locus annotation
    path genome         // genome FASTA (.fai expected alongside)
    path hla            // hla_typing.parquet

    output:
    path "te_antigen.parquet", emit: table
    path "versions.yml",       emit: versions

    script:
    def min_reads = params.te_min_reads ?: 10
    def min_cpm   = params.te_min_cpm   ?: 1.0
    """
    build_te_antigen_table.py \\
        --te-locus ${te_locus} \\
        --annotation ${annotation} \\
        --genome ${genome} \\
        --hla ${hla} \\
        --min-reads ${min_reads} \\
        --min-cpm ${min_cpm} \\
        --out te_antigen.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcflurry: \$(python -c 'import mhcflurry; print(mhcflurry.__version__)')
        biopython: \$(python -c 'import Bio; print(Bio.__version__)')
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    """
    python -c "import pandas as pd; pd.DataFrame({'run_accession':['ERRTEST01'],'cohort':['gide2019'],'te_antigen_burden':[57],'te_antigen_burden_strong':[15],'te_antigen_burden_LINE':[9],'te_antigen_burden_SINE':[13],'te_antigen_burden_LTR':[35],'te_antigen_burden_ERV':[35],'te_antigen_n_expressed_loci':[3],'te_antigen_n_binder_loci':[3],'te_antigen_top_locus':['ERVK_locus_1']}).to_parquet('te_antigen.parquet')"
    echo '"${task.process}":' > versions.yml
    """
}
