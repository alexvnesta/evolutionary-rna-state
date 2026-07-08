//
// Subworkflow: HLA_TYPING
// RNA-based HLA class I genotyping from the rnaseq spine's STAR genome BAMs.
// Primary tool: arcasHLA (RabadanLab/arcasHLA). Fallback: OptiType.
//
// Emits, per sample: the 6 class-I alleles (A/B/C x2) and HLA_I_heterozygous
// (het at all three loci = the Chowell 2017 favorable feature), merged into
// one cohort table in the contract format (rows = run_accession, + cohort).
//
// This is a documented, correct STUB other modules wire into. It is not
// guaranteed to run end-to-end in the arm64 dev sandbox (arcasHLA needs a
// kallisto build + the git-lfs IMGT/HLA reference); the pilot runs it on the
// Linux pipeline host. The stub{} blocks let `nextflow -stub-run` exercise the
// wiring without the tools present. Parsing + heterozygosity logic live in
// analysis/antigen_core/hla_typing.py and are unit-tested there.
//
// Channel contract (matches the intron_retention / rna_editing subworkflows):
//   ch_bam : [ val(meta), path(bam), path(bai) ]   where meta = [id, cohort, single_end]
//

include { ARCASHLA_EXTRACT   } from './modules/arcashla_extract.nf'
include { ARCASHLA_GENOTYPE  } from './modules/arcashla_genotype.nf'
include { MERGE_HLA_TABLE    } from './modules/merge_hla_table.nf'

workflow HLA_TYPING {
    take:
    ch_bam          // channel: [ val(meta), path(bam), path(bai) ]

    main:
    ch_versions = Channel.empty()

    // 1) extract chr6 + unmapped HLA-region reads to FASTQ (arcasHLA extract)
    ARCASHLA_EXTRACT( ch_bam )
    ch_versions = ch_versions.mix( ARCASHLA_EXTRACT.out.versions.first() )

    // 2) genotype HLA-A/B/C per sample from the extracted reads
    ARCASHLA_GENOTYPE( ARCASHLA_EXTRACT.out.reads )
    ch_versions = ch_versions.mix( ARCASHLA_GENOTYPE.out.versions.first() )

    // 3) merge per-sample genotype JSONs -> one cohort tidy table
    //    (parse + heterozygosity via analysis/antigen_core/hla_typing.py)
    MERGE_HLA_TABLE(
        ARCASHLA_GENOTYPE.out.genotype.map { meta, gj -> gj }.collect()
    )
    ch_versions = ch_versions.mix( MERGE_HLA_TABLE.out.versions )

    emit:
    genotype   = ARCASHLA_GENOTYPE.out.genotype   // per-sample genotype.json
    hla_table  = MERGE_HLA_TABLE.out.table        // cohort hla_typing.parquet
    versions   = ch_versions
}

// ---------------------------------------------------------------------------
// Inline module definitions (kept in one file for a self-contained stub; in
// the pipeline these would live under ./modules/*.nf like the sibling
// subworkflows).
// ---------------------------------------------------------------------------

process ARCASHLA_EXTRACT {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/hla_typing/extracted", mode: 'copy'
    conda "bioconda::arcas-hla=0.6.0 bioconda::kallisto"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.extracted.*.fq.gz"), emit: reads
    path "versions.yml",                          emit: versions

    script:
    def single = meta.single_end ? '--single' : ''
    """
    arcasHLA extract ${bam} ${single} -t ${task.cpus} -o . -v

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arcasHLA: \$(arcasHLA --version 2>&1 | sed -n 's/.*version //p')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.extracted.1.fq.gz ${meta.id}.extracted.2.fq.gz
    echo '"${task.process}":' > versions.yml
    """
}

process ARCASHLA_GENOTYPE {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/hla_typing/genotype", mode: 'copy'
    conda "bioconda::arcas-hla=0.6.0 bioconda::kallisto"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.genotype.json"), emit: genotype
    path "versions.yml",                               emit: versions

    script:
    """
    arcasHLA genotype ${reads} -g A,B,C -t ${task.cpus} -o . -v

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arcasHLA: \$(arcasHLA --version 2>&1 | sed -n 's/.*version //p')
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_JSON > ${meta.id}.genotype.json
    {"A": ["A*02:01:01", "A*01:01:01"], "B": ["B*07:02:01", "B*08:01:01"], "C": ["C*07:01:01", "C*07:02:01"]}
    END_JSON
    echo '"${task.process}":' > versions.yml
    """
}

process MERGE_HLA_TABLE {
    label 'process_low'
    publishDir "${params.outdir}/hla_typing", mode: 'copy'
    conda "conda-forge::pandas conda-forge::pyarrow"

    input:
    path genotypes    // collected per-sample *.genotype.json

    output:
    path "hla_typing.parquet", emit: table
    path "versions.yml",       emit: versions

    script:
    // sample_map.csv (run_accession,cohort,genotype_json) is provided by the
    // pipeline from the samplesheet; the merge script keys JSON basename ->
    // run_accession/cohort, then calls hla_typing.summarize_genotype per sample.
    """
    merge_hla_table.py \\
        --genotypes ${genotypes} \\
        --sample-map ${params.hla_sample_map} \\
        --out hla_typing.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    """
    python -c "import pandas as pd; pd.DataFrame({'run_accession':['ERRTEST01'],'cohort':['gide2019'],'HLA_A_1':['A*02:01'],'HLA_A_2':['A*01:01'],'HLA_B_1':['B*07:02'],'HLA_B_2':['B*08:01'],'HLA_C_1':['C*07:01'],'HLA_C_2':['C*07:02'],'HLA_I_heterozygous':[True],'n_het_loci':[3],'tool':['arcasHLA'],'tool_version':['0.6.0']}).to_parquet('hla_typing.parquet')"
    echo '"${task.process}":' > versions.yml
    """
}
