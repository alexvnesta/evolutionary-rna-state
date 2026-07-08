// Alu Editing Index (AEI) per sample: A>G mismatch rate pooled over all Alu
// elements genome-wide (Roth et al., Nat Methods 2019). Robust, cohort-level.
process ALU_EDITING_INDEX {
    tag { meta.id }
    label 'process_editing'
    conda "${moduleDir}/../environment.yml"
    publishDir "${params.outdir}/rna_editing/aei", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path alu_bed

    output:
    tuple val(meta), path("${meta.id}.aei.tsv"), emit: aei
    path "versions.yml",                         emit: versions

    script:
    def snp_arg = params.editing_snp_bed ? "--snp-bed ${params.editing_snp_bed}" : ""
    """
    compute_aei.py \\
        --bam ${bam} \\
        --fasta ${fasta} \\
        --alu ${alu_bed} \\
        --sample ${meta.id} \\
        --min-baseq ${params.editing_min_baseq} \\
        --min-mapq ${params.editing_min_mapq} \\
        ${snp_arg} \\
        --out ${meta.id}.aei.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """

    stub:
    """
    printf 'sample\\tAEI_percent\\tAG_mismatches\\tA_coverage\\tsignal_to_noise\\n${meta.id}\\t1.234\\t100\\t8100\\t12.3\\n' > ${meta.id}.aei.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}

// Merge per-sample AEI TSVs into one cohort table.
process MERGE_AEI {
    label 'process_light'
    conda "${moduleDir}/../environment.yml"
    publishDir "${params.outdir}/rna_editing/aei", mode: 'copy'

    input:
    path aei_tsvs

    output:
    path "cohort_aei.tsv", emit: cohort

    script:
    """
    merge_aei.py --out cohort_aei.tsv ${aei_tsvs}
    """

    stub:
    """
    touch cohort_aei.tsv
    """
}
