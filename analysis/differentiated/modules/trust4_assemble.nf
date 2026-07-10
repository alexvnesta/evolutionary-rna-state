//
// Module: TRUST4_ASSEMBLE
// Runs TRUST4 (built from source, tools/TRUST4-1.1.5) on one sample's paired
// reads to reconstruct BCR/TCR contigs. Emits the sample's TRUST4 output dir
// (<run>_cdr3.out / _report.tsv / _airr.tsv). Stub-runnable without TRUST4.
//
// NOTE: for the stream-from-ENA pilot arm, the streaming + this assembly are
// done together by pipelines/bcr_repertoire/run_trust4_pilot.sh (peak disk ~1
// sample). This process is the on-disk-reads path for a full Nextflow run.
//
process TRUST4_ASSEMBLE {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/bcr_repertoire", mode: 'copy'
    conda "bioconda::trust4=1.1.5"

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}"), emit: trust4_dir
    path "versions.yml",                 emit: versions

    script:
    def t4dir = params.trust4_dir ?: 'tools/TRUST4-1.1.5'
    """
    run-trust4 \\
        -1 ${r1} -2 ${r2} \\
        -f ${t4dir}/hg38_bcrtcr.fa \\
        --ref ${t4dir}/human_IMGT+C.fa \\
        -t ${task.cpus} \\
        -o ${meta.id} --od ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trust4: \$(run-trust4 2>&1 | head -1 | sed 's/.*v//;s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}
    printf 'c0\\t0\\tIGHV1-2\\t*\\tIGHJ4\\tIGHG1\\tAAA\\tBBB\\tCARWYFDVW\\t1.00\\t50\\t0.9000\\t1\\n' > ${meta.id}/${meta.id}_cdr3.out
    echo '"${task.process}":' > versions.yml
    """
}
