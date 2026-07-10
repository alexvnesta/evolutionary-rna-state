//
// Module: BCR_SHM_FEATURES
// Turns one sample's TRUST4 output dir into the NAMED BCR features
// (bcr_shm_rate, bcr_igg_fraction, bcr_switched_fraction, bcr_clonality,
// bcr_n_clonotypes, bcr_n_reads) by running
// analysis/differentiated/bin/bcr_shm_features.py, which imports the
// unit-tested logic in bcr_shm.py. Stub-runnable without pandas present.
//
process BCR_SHM_FEATURES {
    tag "${meta.id}"
    label 'process_single'
    conda "conda-forge::pandas conda-forge::numpy conda-forge::pyarrow"

    input:
    tuple val(meta), path(trust4_dir)

    output:
    tuple val(meta), path("${meta.id}_bcr_row.csv"), emit: row
    path "versions.yml",                             emit: versions

    script:
    def minc = params.min_clonotypes ?: 3
    def cohort = meta.cohort ? "--cohort ${meta.cohort}" : ""
    """
    bcr_shm_features.py \\
        --trust4-dir ${trust4_dir} \\
        --run-accession ${meta.id} ${cohort} \\
        --min-clonotypes ${minc} \\
        --out ${meta.id}_bcr_row.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf 'run_accession,cohort,bcr_shm_rate,bcr_igg_fraction,bcr_switched_fraction,bcr_clonality,bcr_n_clonotypes,bcr_n_reads,_shm_source,_shm_n\\n${meta.id},${meta.cohort ?: "NA"},0.05,0.8,0.9,0.1,10,40,cdr3_germline_similarity,10\\n' > ${meta.id}_bcr_row.csv
    echo '"${task.process}":' > versions.yml
    """
}
