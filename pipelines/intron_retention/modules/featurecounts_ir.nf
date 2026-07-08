process FEATURECOUNTS_IR {
    tag "${meta.id}:${feature}"
    label 'process_medium'
    publishDir "${params.outdir}/intron_retention/counts", mode: 'copy'
    conda "${moduleDir}/../environment.yml"

    input:
    tuple val(meta), path(bam), path(bai)
    each path(saf)          // introns.saf or exons.saf
    val  feature            // "intron" or "exon" (used only for output naming)

    output:
    tuple val(meta), val(feature), path("${meta.id}.${feature}.featureCounts.txt"), emit: counts
    path "${meta.id}.${feature}.featureCounts.txt.summary",                         emit: summary
    path "versions.yml",                                                            emit: versions

    script:
    // -O + --fraction: intronic reads overlapping >1 meta-feature are fractionally
    //   assigned rather than discarded (retained introns can span nested features).
    // -p --countReadPairs: paired-end fragment counting when reads are paired.
    // Strandedness comes from params.ir_strandedness (0 unstranded, 1 fwd, 2 rev).
    def paired = meta.single_end ? '' : '-p --countReadPairs'
    def strand = params.ir_strandedness
    """
    featureCounts \\
        -a ${saf} -F SAF \\
        -f -O --fraction \\
        -s ${strand} \\
        ${paired} \\
        -T ${task.cpus} \\
        -o ${meta.id}.${feature}.featureCounts.txt \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        featurecounts: \$(featureCounts -v 2>&1 | sed -n 's/.*featureCounts v//p')
    END_VERSIONS
    """

    stub:
    """
    printf '# stub\\nGeneid\\tChr\\tStart\\tEnd\\tStrand\\tLength\\t${meta.id}.bam\\nENSG_test__intron_1\\tchr1\\t12228\\t12612\\t+\\t385\\t42\\n' > ${meta.id}.${feature}.featureCounts.txt
    touch ${meta.id}.${feature}.featureCounts.txt.summary
    echo '"${task.process}":' > versions.yml
    """
}
