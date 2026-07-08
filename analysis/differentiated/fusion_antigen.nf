//
// Subworkflow: FUSION_ANTIGEN
// Fusion-transcript-derived neoantigen burden from the rnafusion pipeline's
// per-sample fusion calls + the antigen_core HLA table.
//
// Emits, per sample: fusion_neoantigen_burden (unique MHC-I binders across
// in-frame fusion junctions, via the SHARED antigen_core MHCflurry engine),
// plus the named interpretable features n_fusions and n_inframe_fusions.
//
// FIXED-CALLER REQUIREMENT (batch robustness): a run uses ONE caller for every
// sample (params.fusion_caller = 'arriba' | 'starfusion'); the caller + version
// are recorded on every output row. Fusion detection is caller/version/depth
// sensitive, so the caller is held constant across a comparison. See
// BATCH_ROBUSTNESS_NOTE in analysis/differentiated/fusion_antigen.py.
//
// This is a documented, correct STUB that wires into the existing pattern
// (mirrors analysis/antigen_core/hla_typing.nf). The peptide-derivation and
// burden logic live in analysis/differentiated/fusion_antigen.py and are
// unit-tested in test_fusion_antigen.py. The stub{} blocks let
// `nextflow -stub-run` exercise the wiring without MHCflurry/arriba present.
//
// Channel contract (matches the sibling subworkflows):
//   ch_fusions : [ val(meta), path(fusion_tsv) ]   meta = [id, cohort]
//                fusion_tsv = Arriba fusions.tsv OR STAR-Fusion coding TSV
//   ch_hla     : path(hla_typing.parquet)          from HLA_TYPING subworkflow
//

include { FUSION_NEOANTIGEN_BURDEN } from './modules/fusion_neoantigen_burden.nf'
include { MERGE_FUSION_FEATURES    } from './modules/merge_fusion_features.nf'

workflow FUSION_ANTIGEN {
    take:
    ch_fusions      // channel: [ val(meta), path(fusion_tsv) ]
    ch_hla          // channel: path(hla_typing.parquet)   (single cohort table)

    main:
    ch_versions = Channel.empty()

    // 1) per-sample: derive junction peptides -> shared MHC engine -> burden row
    FUSION_NEOANTIGEN_BURDEN(
        ch_fusions,
        ch_hla,
    )
    ch_versions = ch_versions.mix( FUSION_NEOANTIGEN_BURDEN.out.versions.first() )

    // 2) merge per-sample feature rows -> one cohort fusion_features table
    MERGE_FUSION_FEATURES(
        FUSION_NEOANTIGEN_BURDEN.out.row.map { meta, csv -> csv }.collect()
    )
    ch_versions = ch_versions.mix( MERGE_FUSION_FEATURES.out.versions )

    emit:
    row      = FUSION_NEOANTIGEN_BURDEN.out.row     // per-sample feature CSV
    table    = MERGE_FUSION_FEATURES.out.table      // cohort fusion_features.parquet
    versions = ch_versions
}

// ---------------------------------------------------------------------------
// Inline module definitions (kept in one file for a self-contained stub; in
// the pipeline these would live under ./modules/*.nf like the sibling
// subworkflows).
// ---------------------------------------------------------------------------

process FUSION_NEOANTIGEN_BURDEN {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/fusion_antigen/per_sample", mode: 'copy'
    // the 'antigen' env: mhcflurry + biopython + pandas + pyarrow
    conda "bioconda::mhcflurry=2.2.0 conda-forge::biopython conda-forge::pandas conda-forge::pyarrow"

    input:
    tuple val(meta), path(fusion_tsv)
    path  hla_parquet

    output:
    tuple val(meta), path("${meta.id}.fusion_features.csv"), emit: row
    path  "versions.yml",                                     emit: versions

    script:
    // params.fusion_caller fixes the caller for the whole run (batch robustness).
    // The HLA-I alleles for this sample are pulled from the cohort HLA table by
    // (run_accession == meta.id). MHCFLURRY_DATA_DIR points at the repo models.
    def caller_flag = params.fusion_caller == 'starfusion' ? '--starfusion' : '--arriba'
    """
    export MHCFLURRY_DATA_DIR=${params.mhcflurry_models}
    fusion_burden_cli.py \\
        ${caller_flag} ${fusion_tsv} \\
        --hla-table ${hla_parquet} \\
        --run-accession ${meta.id} \\
        --cohort ${meta.cohort} \\
        --caller-version '${params.fusion_caller_version}' \\
        --out ${meta.id}.fusion_features.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcflurry: \$(python -c 'import mhcflurry; print(mhcflurry.__version__)')
        caller: ${params.fusion_caller} ${params.fusion_caller_version}
    END_VERSIONS
    """

    stub:
    // wiring check only: emit one contract-shaped row with an integer burden
    """
    cat <<-END_CSV > ${meta.id}.fusion_features.csv
    run_accession,cohort,n_fusions,n_inframe_fusions,fusion_neoantigen_burden,fusion_neoantigen_burden_strong,caller,caller_version
    ${meta.id},${meta.cohort},2,1,6,3,arriba,2.4.0
    END_CSV
    echo '"${task.process}":' > versions.yml
    """
}

process MERGE_FUSION_FEATURES {
    label 'process_low'
    publishDir "${params.outdir}/fusion_antigen", mode: 'copy'
    conda "conda-forge::pandas conda-forge::pyarrow"

    input:
    path rows    // collected per-sample *.fusion_features.csv

    output:
    path "fusion_features.parquet", emit: table
    path "versions.yml",            emit: versions

    script:
    """
    python - <<'PY'
    import glob, pandas as pd
    df = pd.concat([pd.read_csv(f) for f in glob.glob("*.fusion_features.csv")],
                   ignore_index=True)
    df.to_parquet("fusion_features.parquet", index=False)
    PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    """
    python -c "import pandas as pd; pd.DataFrame({'run_accession':['ERRTEST01'],'cohort':['gide2019'],'n_fusions':[2],'n_inframe_fusions':[1],'fusion_neoantigen_burden':[6],'fusion_neoantigen_burden_strong':[3],'caller':['arriba'],'caller_version':['2.4.0']}).to_parquet('fusion_features.parquet')"
    echo '"${task.process}":' > versions.yml
    """
}
