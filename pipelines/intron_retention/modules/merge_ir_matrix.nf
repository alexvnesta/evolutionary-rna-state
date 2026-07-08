process MERGE_IR_MATRIX {
    tag "merge"
    label 'process_low'
    publishDir "${params.outdir}/features", mode: 'copy'
    conda "${moduleDir}/../environment.yml"

    input:
    path long_tables   // all *.ir_long.tsv
    path summaries     // all *.ir_summary.tsv

    output:
    path "intron_retention.parquet", emit: matrix
    path "intron_retention_summary.tsv", emit: summary
    path "versions.yml", emit: versions

    script:
    """
    merge_ir_matrix.py \\
        --long ${long_tables} \\
        --summaries ${summaries} \\
        --out-matrix intron_retention.parquet \\
        --out-summary intron_retention_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    """
    touch intron_retention.parquet
    printf 'run_accession\\tcohort\\tn_introns_evaluated\\n' > intron_retention_summary.tsv
    echo '"${task.process}":' > versions.yml
    """
}
