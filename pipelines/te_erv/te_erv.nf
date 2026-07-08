/*
 * ============================================================================
 *  TE / ERV activation subworkflow  (DSL2)
 * ----------------------------------------------------------------------------
 *  Quantifies transposable-element and endogenous-retrovirus expression at
 *  LOCUS level (Telescope, EM reassignment of multimappers) and FAMILY level
 *  (TEtranscripts / TEcount) for the evolutionary-RNA-state melanoma-ICB study.
 *
 *  Two independent alignment inputs, by design:
 *    - Telescope    needs its OWN STAR run with high multimapping retention
 *                   (--outFilterMultimapNmax 100 --winAnchorMultimapNmax 200,
 *                    unsorted BAM keeping all alignments of a read together).
 *                   The standard nf-core/rnaseq BAM discards multimappers and
 *                   is therefore NOT usable for locus-level TE EM.
 *    - TEtranscripts consumes the coordinate-sorted genome BAM already produced
 *                   by the nf-core/rnaseq hisat2 spine (multi mode handles
 *                   the multimappers statistically from the standard BAM).
 *
 *  Resource directives are set inline (this Mac: 40 GB / maxForks 1 for STAR,
 *  <=10 GB for the light steps) so the module is self-contained regardless of
 *  the withName patterns in conf/mac_arm64.config.
 * ============================================================================
 */

nextflow.enable.dsl = 2

// ----------------------------------------------------------------------------
// PROCESS: STAR_ALIGN_MULTI
//   Multimap-aware STAR alignment for Telescope. Consumes raw FASTQs because
//   the required multimapping retention differs from the rnaseq spine.
// ----------------------------------------------------------------------------
process STAR_ALIGN_MULTI {
    tag   { meta.id }
    label 'process_star'
    publishDir "${params.outdir}/star_multi", mode: 'copy',
               saveAs: { fn -> fn.endsWith('.bam') ? fn : "logs/${fn}" }

    cpus     12
    memory   40.GB
    maxForks 1
    time     24.h

    input:
    tuple val(meta), path(reads)
    path  star_index

    output:
    tuple val(meta), path("${meta.id}.Aligned.out.bam"), emit: bam
    path  "${meta.id}.Log.final.out",                    emit: log
    path  "versions.yml",                                emit: versions

    script:
    def readsArg = meta.single_end ? "${reads}" : "${reads[0]} ${reads[1]}"
    def gzipped  = reads instanceof List ? reads[0].name.endsWith('.gz') : reads.name.endsWith('.gz')
    def readCmd  = gzipped ? '--readFilesCommand zcat' : ''
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${star_index} \\
        --readFilesIn ${readsArg} \\
        ${readCmd} \\
        --outFilterMultimapNmax 100 \\
        --winAnchorMultimapNmax 200 \\
        --outSAMtype BAM Unsorted \\
        --outSAMattributes NH HI AS nM \\
        --outSAMprimaryFlag AllBestScore \\
        --outFileNamePrefix ${meta.id}. \\
        --limitBAMsortRAM 0 \\
        ${params.star_extra_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.Aligned.out.bam ${meta.id}.Log.final.out
    echo '"${task.process}": {star: stub}' > versions.yml
    """
}

// ----------------------------------------------------------------------------
// PROCESS: TELESCOPE_ASSIGN
//   Locus-level EM reassignment of multimapped reads to TE/ERV loci.
//   pip-installed (telescope-ngs / mlbendall telescope), pure-Python + Cython.
// ----------------------------------------------------------------------------
process TELESCOPE_ASSIGN {
    tag   { meta.id }
    label 'process_medium'
    publishDir "${params.outdir}/telescope", mode: 'copy'

    cpus     4
    memory   10.GB
    time     12.h

    input:
    tuple val(meta), path(bam)
    path  te_gtf_locus

    output:
    tuple val(meta), path("${meta.id}-telescope_report.tsv"), emit: report
    path  "${meta.id}-*.log",                                 emit: log
    path  "versions.yml",                                     emit: versions

    script:
    """
    telescope assign \\
        --exp_tag ${meta.id} \\
        --theta_prior 200000 \\
        --max_iter 200 \\
        --outdir . \\
        ${params.telescope_extra_args} \\
        ${bam} \\
        ${te_gtf_locus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telescope: \$(telescope --version 2>&1 | sed 's/^.*telescope //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}-telescope_report.tsv ${meta.id}-telescope.log
    echo '"${task.process}": {telescope: stub}' > versions.yml
    """
}

// ----------------------------------------------------------------------------
// PROCESS: TETRANSCRIPTS_COUNT  (TEcount, family-level)
//   Consumes the coordinate-sorted genome BAM from the rnaseq spine.
//   noarch conda package.
// ----------------------------------------------------------------------------
process TETRANSCRIPTS_COUNT {
    tag   { meta.id }
    label 'process_medium'
    publishDir "${params.outdir}/tetranscripts", mode: 'copy'

    cpus     4
    memory   10.GB
    time     12.h

    input:
    tuple val(meta), path(bam)
    path  gene_gtf
    path  te_gtf_family

    output:
    tuple val(meta), path("${meta.id}.cntTable"), emit: counts
    path  "versions.yml",                          emit: versions

    script:
    def stranded = meta.strandedness ?: 'no'
    """
    TEcount \\
        --BAM ${bam} \\
        --GTF ${gene_gtf} \\
        --TE ${te_gtf_family} \\
        --mode multi \\
        --stranded ${stranded} \\
        --sortByPos \\
        --project ${meta.id} \\
        ${params.tetranscripts_extra_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tetranscripts: \$(TEcount --version 2>&1 | sed 's/^.*TEcount v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    printf 'gene/TE\\t${meta.id}\\n' > ${meta.id}.cntTable
    echo '"${task.process}": {tetranscripts: stub}' > versions.yml
    """
}

// ----------------------------------------------------------------------------
// PROCESS: MERGE_TELESCOPE   -> locus-level TE/ERV count matrix
// ----------------------------------------------------------------------------
process MERGE_TELESCOPE {
    label 'process_low'
    publishDir "${params.outdir}/matrices", mode: 'copy'

    cpus   2
    memory 8.GB

    input:
    path reports   // collect() of *-telescope_report.tsv

    output:
    path "telescope_locus_counts.tsv", emit: matrix

    script:
    """
    merge_matrices.py telescope ${reports} > telescope_locus_counts.tsv
    """

    stub:
    """
    echo -e "locus\\tsampleA" > telescope_locus_counts.tsv
    """
}

// ----------------------------------------------------------------------------
// PROCESS: MERGE_TETRANSCRIPTS  -> family-level count matrix
// ----------------------------------------------------------------------------
process MERGE_TETRANSCRIPTS {
    label 'process_low'
    publishDir "${params.outdir}/matrices", mode: 'copy'

    cpus   2
    memory 8.GB

    input:
    path cnt_tables   // collect() of *.cntTable

    output:
    path "tetranscripts_family_counts.tsv", emit: matrix

    script:
    """
    merge_matrices.py tetranscripts ${cnt_tables} > tetranscripts_family_counts.tsv
    """

    stub:
    """
    echo -e "feature\\tsampleA" > tetranscripts_family_counts.tsv
    """
}

// ============================================================================
//  SUBWORKFLOW
// ============================================================================
workflow TE_ERV {

    take:
    ch_reads          // channel: [ val(meta), [ fastq(s) ] ]  (for Telescope)
    ch_bam            // channel: [ val(meta), path(coord_sorted_bam) ] (rnaseq spine, for TEtranscripts)
    star_index        // path
    te_gtf_locus      // path : Telescope retro.hg38.v1 transcripts.gtf
    gene_gtf          // path : GENCODE genic GTF
    te_gtf_family     // path : TEtranscripts family TE GTF

    main:
    ch_versions = Channel.empty()

    // ---- Locus-level branch: dedicated multimap-aware STAR -> Telescope ----
    STAR_ALIGN_MULTI ( ch_reads, star_index )
    ch_versions = ch_versions.mix( STAR_ALIGN_MULTI.out.versions.first() )

    TELESCOPE_ASSIGN ( STAR_ALIGN_MULTI.out.bam, te_gtf_locus )
    ch_versions = ch_versions.mix( TELESCOPE_ASSIGN.out.versions.first() )

    MERGE_TELESCOPE ( TELESCOPE_ASSIGN.out.report.map { meta, tsv -> tsv }.collect() )

    // ---- Family-level branch: reuse spine BAM -> TEcount ----
    TETRANSCRIPTS_COUNT ( ch_bam, gene_gtf, te_gtf_family )
    ch_versions = ch_versions.mix( TETRANSCRIPTS_COUNT.out.versions.first() )

    MERGE_TETRANSCRIPTS ( TETRANSCRIPTS_COUNT.out.counts.map { meta, tsv -> tsv }.collect() )

    emit:
    telescope_reports = TELESCOPE_ASSIGN.out.report      // per-sample locus reports
    te_counts         = TETRANSCRIPTS_COUNT.out.counts   // per-sample family counts
    locus_matrix      = MERGE_TELESCOPE.out.matrix        // locus-level matrix
    family_matrix     = MERGE_TETRANSCRIPTS.out.matrix    // family-level matrix
    versions          = ch_versions
}
