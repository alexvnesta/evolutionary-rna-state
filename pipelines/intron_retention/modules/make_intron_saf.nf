process MAKE_INTRON_SAF {
    tag "derive_introns"
    label 'process_low'
    publishDir "${params.outdir}/intron_retention/reference", mode: 'copy'
    conda "${moduleDir}/../environment.yml"

    input:
    path gtf

    output:
    path "introns.saf",       emit: introns_saf
    path "exons.saf",         emit: exons_saf
    path "intron2gene.tsv",   emit: map
    path "versions.yml",      emit: versions

    script:
    """
    make_intron_saf.py \\
        --gtf ${gtf} \\
        --out-introns introns.saf \\
        --out-exons exons.saf \\
        --out-map intron2gene.tsv \\
        --min-intron-len ${params.ir_min_intron_len}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf 'GeneID\\tChr\\tStart\\tEnd\\tStrand\\nENSG_test__intron_1\\tchr1\\t12228\\t12612\\t+\\n' > introns.saf
    printf 'GeneID\\tChr\\tStart\\tEnd\\tStrand\\nENSG_test\\tchr1\\t11869\\t12227\\t+\\n' > exons.saf
    printf 'intron_id\\tgene_id\\tintron_length\\nENSG_test__intron_1\\tENSG_test\\t385\\n' > intron2gene.tsv
    echo '"${task.process}":' > versions.yml
    """
}
