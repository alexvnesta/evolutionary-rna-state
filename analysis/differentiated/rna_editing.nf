//
// Subworkflow: EDITING_ANTIGEN
// Derives the two NAMED per-sample editing features on top of the pipeline
// session's RNA-editing quantification (pipelines/rna_editing):
//
//   alu_editing_index          (float [0,1]) — the primary, batch-robust
//                                ADAR A-to-I "dial": AG_mismatches / A_coverage
//                                pooled over all Alu adenosines (Roth/Levanon
//                                2019). Read straight from the pipeline's
//                                per-sample *.aei.tsv (compute_aei.py) — no
//                                recomputation, so it stays identical to the
//                                quantification-layer rna_editing_aei column.
//   editing_neoantigen_burden  (int)         — MHCflurry binder count over the
//                                altered (I=G) recoding peptides from CDS
//                                editing sites (REDIportal recoding catalog
//                                intersected with the sample's *.editing_sites.tsv),
//                                scored through the SHARED antigen_core engine.
//
// This is a documented, correct STUB the pilot wires into. Feature logic +
// heterozygosity/recoding translation are unit-tested in
// analysis/differentiated/test_rna_editing.py (no tools required). The MHCflurry
// engine + REDIportal catalog build run on the pipeline host at pilot time; the
// stub{} blocks let `nextflow -stub-run` exercise the wiring without them.
//
// Channel contract (matches hla_typing / intron_retention / te_erv):
//   ch_aei   : [ val(meta), path(aei_tsv) ]     from RNA_EDITING.out.aei
//   ch_sites : [ val(meta), path(sites_tsv) ]   from RNA_EDITING.out.sites
//   meta = [id, cohort]. run_accession == meta.id.
//
// Upstream (produced by pipelines/rna_editing/main.nf):
//   RNA_EDITING.out.aei   -> per-sample <id>.aei.tsv
//   RNA_EDITING.out.sites -> per-sample <id>.editing_sites.tsv
//

include { BUILD_RECODING_CATALOG } from './modules/build_recoding_catalog.nf'
include { EDITING_FEATURES        } from './modules/editing_features.nf'
include { MERGE_EDITING_FEATURES  } from './modules/editing_features.nf'

workflow EDITING_ANTIGEN {
    take:
    ch_aei          // channel: [ val(meta), path(aei_tsv) ]
    ch_sites        // channel: [ val(meta), path(editing_sites_tsv) ]
    hla_table       // path: hla_typing.parquet (from HLA_TYPING subworkflow)
    reditportal     // path: REDIportal coding/recoding table (TABLE1_*.txt[.gz])
    gtf             // path: GENCODE GTF (CDS models for the catalog)
    fasta           // path: genome FASTA (to fill each site's reference CDS window)
    fai             // path: FASTA index

    main:
    ch_versions = Channel.empty()

    // 1) Build the fixed recoding-site catalog ONCE (site_id, chrom, pos,
    //    strand, gene, cds_window, edit_offset). Reused for every sample so
    //    recoding detection is not depth-driven (batch-robustness lever).
    BUILD_RECODING_CATALOG( reditportal, gtf, fasta, fai )
    ch_catalog  = BUILD_RECODING_CATALOG.out.catalog
    ch_versions = ch_versions.mix( BUILD_RECODING_CATALOG.out.versions )

    // 2) Per-sample editing features: alu_editing_index (from aei_tsv) +
    //    editing_neoantigen_burden (recoding catalog x editing_sites -> shared
    //    MHCflurry engine, with the sample's HLA-I alleles from hla_table).
    ch_feat_in = ch_aei.join( ch_sites )          // [ meta, aei_tsv, sites_tsv ]
    EDITING_FEATURES( ch_feat_in, ch_catalog, hla_table )
    ch_versions = ch_versions.mix( EDITING_FEATURES.out.versions.first() )

    // 3) Merge per-sample rows into the contract feature matrix.
    MERGE_EDITING_FEATURES(
        EDITING_FEATURES.out.feature.map { meta, tsv -> tsv }.collect()
    )
    ch_versions = ch_versions.mix( MERGE_EDITING_FEATURES.out.versions )

    emit:
    per_sample = EDITING_FEATURES.out.feature       // [ meta, <id>.editing_features.tsv ]
    features   = MERGE_EDITING_FEATURES.out.matrix  // rna_editing_features.parquet
    catalog    = ch_catalog                         // recoding_catalog.tsv
    versions   = ch_versions
}

// ---------------------------------------------------------------------------
// Inline module definitions (kept in one file for a self-contained stub; in
// the pipeline these live under ./modules/*.nf like the sibling subworkflows).
// ---------------------------------------------------------------------------

process BUILD_RECODING_CATALOG {
    label 'process_low'
    conda "conda-forge::pandas conda-forge::pysam bioconda::pyensembl"
    publishDir "${params.outdir}/editing_antigen", mode: 'copy'

    input:
    path reditportal    // REDIportal coding sites (recoding) table
    path gtf            // GENCODE GTF for CDS models
    path fasta          // genome FASTA
    path fai            // .fai

    output:
    path "recoding_catalog.tsv", emit: catalog
    path "versions.yml",         emit: versions

    script:
    """
    build_recoding_catalog.py \\
        --reditportal ${reditportal} \\
        --gtf ${gtf} \\
        --fasta ${fasta} \\
        --codon-flank ${params.editing_codon_flank ?: 10} \\
        --out recoding_catalog.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """

    stub:
    // A single placeholder catalog row for -stub-run wiring: the canonical
    // GRIA2 Q607R (Q->R) ADAR recoding site. Illustrative only (the cds_window
    // is elided); real windows are filled by build_recoding_catalog.py from the
    // genome FASTA at pilot time.
    """
    printf 'site_id\\tchrom\\tpos\\tstrand\\tgene\\tcds_window\\tedit_offset\\n' > recoding_catalog.tsv
    printf 'GRIA2_Q607R\\tchr4\\t157336723\\t+\\tGRIA2\\tGGTGGTGGT...CAA...GGTGGT\\t32\\n' >> recoding_catalog.tsv
    echo '"${task.process}":' > versions.yml
    """
}

process EDITING_FEATURES {
    tag { meta.id }
    label 'process_medium'
    // env with pandas + mhcflurry + the fetched models (antigen_core engine).
    conda "${moduleDir}/../environment.yml"
    publishDir "${params.outdir}/editing_antigen/per_sample", mode: 'copy'

    input:
    tuple val(meta), path(aei_tsv), path(sites_tsv)
    path catalog
    path hla_table

    output:
    tuple val(meta), path("${meta.id}.editing_features.tsv"), emit: feature
    path "versions.yml",                                      emit: versions

    script:
    def freq = params.editing_min_freq ?: 0.10
    """
    # Thin CLI over analysis/differentiated/rna_editing.py: computes
    # alu_editing_index from the AEI tsv and editing_neoantigen_burden from the
    # recoding catalog x this sample's editing sites, using the sample's HLA-I
    # alleles (looked up from hla_table by run_accession == meta.id).
    editing_features_cli.py \\
        --sample ${meta.id} \\
        --cohort ${meta.cohort} \\
        --aei ${aei_tsv} \\
        --sites ${sites_tsv} \\
        --recoding-catalog ${catalog} \\
        --hla-table ${hla_table} \\
        --freq-threshold ${freq} \\
        --out ${meta.id}.editing_features.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcflurry: \$(python -c "import mhcflurry; print(mhcflurry.__version__)")
    END_VERSIONS
    """

    stub:
    """
    printf 'run_accession\\tcohort\\talu_editing_index\\talu_editing_index_percent\\taei_signal_to_noise\\tediting_neoantigen_burden\\tn_recoding_sites_edited\\tn_candidate_peptides\\n' > ${meta.id}.editing_features.tsv
    printf '${meta.id}\\t${meta.cohort}\\t0.0150\\t1.5000\\t30.000\\t1\\t1\\t38\\n' >> ${meta.id}.editing_features.tsv
    echo '"${task.process}":' > versions.yml
    """
}

process MERGE_EDITING_FEATURES {
    label 'process_low'
    conda "conda-forge::pandas conda-forge::pyarrow"
    publishDir "${params.outdir}/../features", mode: 'copy'   // -> results/features/

    input:
    path feature_tsvs    // collected per-sample *.editing_features.tsv

    output:
    path "rna_editing_features.parquet", emit: matrix
    path "versions.yml",                 emit: versions

    script:
    """
    python - <<'PY'
    import glob, pandas as pd
    frames = [pd.read_csv(f, sep="\\t") for f in sorted(glob.glob("*.editing_features.tsv"))]
    df = pd.concat(frames, ignore_index=True).sort_values(["cohort", "run_accession"])
    # Contract: first two columns run_accession, cohort; keyed per-sample.
    lead = ["run_accession", "cohort"]
    df = df[lead + [c for c in df.columns if c not in lead]]
    df.to_parquet("rna_editing_features.parquet", index=False)
    print(df[["run_accession","cohort","alu_editing_index","editing_neoantigen_burden"]].to_string(index=False))
    PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    """
    python -c "import pandas as pd; pd.DataFrame({'run_accession':['ERRTEST01'],'cohort':['gide2019'],'alu_editing_index':[0.015],'editing_neoantigen_burden':[1]}).to_parquet('rna_editing_features.parquet')"
    echo '"${task.process}":' > versions.yml
    """
}
