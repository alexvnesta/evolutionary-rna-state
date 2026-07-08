process COMPUTE_IR_RATIO {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/intron_retention", mode: 'copy'
    conda "${moduleDir}/../environment.yml"

    input:
    tuple val(meta), path(intron_counts), path(exon_counts)
    path intron_map

    output:
    tuple val(meta), path("${meta.id}.ir_long.tsv"),      emit: long
    tuple val(meta), path("${meta.id}.ir_ratio.parquet"), emit: wide
    tuple val(meta), path("${meta.id}.ir_summary.tsv"),   emit: summary
    path "versions.yml",                                  emit: versions

    script:
    def cohort = meta.cohort ?: 'NA'
    """
    compute_ir_ratio.py \\
        --intron-counts ${intron_counts} \\
        --exon-counts ${exon_counts} \\
        --map ${intron_map} \\
        --run-accession ${meta.id} \\
        --cohort ${cohort} \\
        --min-gene-exon-count ${params.ir_min_gene_exon_count} \\
        --high-ir-threshold ${params.ir_high_threshold} \\
        --out-long ${meta.id}.ir_long.tsv \\
        --out-wide ${meta.id}.ir_ratio.parquet \\
        --out-summary ${meta.id}.ir_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    """
    printf 'run_accession\\tcohort\\tintron_id\\tgene_id\\tIR_ratio\\n${meta.id}\\t${meta.cohort ?: "NA"}\\tENSG_test__intron_1\\tENSG_test\\t0.12\\n' > ${meta.id}.ir_long.tsv
    touch ${meta.id}.ir_ratio.parquet
    printf 'run_accession\\tcohort\\tn_introns_evaluated\\tmedian_IR\\n${meta.id}\\t${meta.cohort ?: "NA"}\\t1\\t0.12\\n' > ${meta.id}.ir_summary.tsv
    echo '"${task.process}":' > versions.yml
    """
}
