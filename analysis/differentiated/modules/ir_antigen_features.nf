//
// Module: IR_ANTIGEN_FEATURES
// Turns intron_retention.parquet + introns.saf + genome + hla_typing.parquet
// into the two NAMED features (retained_intron_load, ir_neoantigen_burden) by
// running analysis/differentiated/bin/ir_antigen_features.py, which imports the
// unit-tested feature logic in intron_retention.py and the SHARED antigen_core
// MHCflurry engine. Stub-runnable without the genome / models present.
//
process IR_ANTIGEN_FEATURES {
    label 'process_medium'
    publishDir "${params.outdir}/features", mode: 'copy'
    conda "conda-forge::pandas conda-forge::pyarrow conda-forge::pysam bioconda::mhcflurry"

    input:
    path ir_matrix
    path intron_saf
    path genome
    path hla_table

    output:
    path "ir_antigen_features.parquet", emit: features
    path "versions.yml",                emit: versions

    script:
    def thr = params.ir_retained_threshold ?: 0.10
    """
    ir_antigen_features.py \\
        --ir-matrix ${ir_matrix} \\
        --intron-saf ${intron_saf} \\
        --genome ${genome} \\
        --hla-table ${hla_table} \\
        --threshold ${thr} \\
        --out ir_antigen_features.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        mhcflurry: \$(python3 -c 'import mhcflurry; print(mhcflurry.__version__)')
    END_VERSIONS
    """

    stub:
    """
    # dependency-free wiring stub: emit the declared output so `-stub-run`
    # exercises the channel graph without the genome / MHCflurry models present.
    touch ir_antigen_features.parquet
    echo '"${task.process}":' > versions.yml
    """
}
