/*
 * ============================================================================
 *  TE / ERV activation subworkflow  (DSL2)
 * ----------------------------------------------------------------------------
 *  Quantifies transposable-element and endogenous-retrovirus expression at
 *  LOCUS level (Telescope, EM reassignment of multimappers) and FAMILY level
 *  (TEtranscripts / TEcount) for the evolutionary-RNA-state melanoma-ICB study.
 *
 *  Two independent alignment inputs, by design:
 *    - Telescope    needs its OWN multimap-permissive alignment. Telescope's own
 *                   documentation and paper (Bendall et al. 2019, PLOS Comp Biol)
 *                   specify bowtie2 with multimapping retained:
 *                     bowtie2 -k 100 --very-sensitive-local --score-min L,0,1.6
 *                   We use exactly that. (The original design used a dedicated
 *                   multimap-aware STAR pass, but the only osx-arm64 conda build
 *                   of STAR 2.7.11b ingests 0 reads — see pipelines/scripts/
 *                   run_rnaseq.sh — so STAR is unusable on this machine. bowtie2
 *                   is arm64-native, is Telescope's canonical aligner, and is
 *                   what the Telescope paper validated multimapper EM against.)
 *                   The standard nf-core/rnaseq spine BAM discards multimappers
 *                   and is therefore NOT usable for locus-level TE EM.
 *    - TEtranscripts consumes the coordinate-sorted genome BAM already produced
 *                   by the nf-core/rnaseq hisat2 spine (multi mode handles
 *                   the multimappers statistically from the standard BAM).
 *
 *  Resource directives are set inline (this Mac: <=16 GB for bowtie2/Telescope,
 *  <=10 GB for the light steps) so the module is self-contained regardless of
 *  the withName patterns in conf/mac_arm64.config.
 * ============================================================================
 */

nextflow.enable.dsl = 2

// ----------------------------------------------------------------------------
// PROCESS: BOWTIE2_ALIGN_MULTI
//   Multimap-permissive bowtie2 alignment for Telescope (Bendall et al. params).
//   Consumes raw FASTQs because the required multimapping retention differs from
//   the rnaseq spine. bowtie2 is arm64-native (the osx-arm64 STAR build is broken)
//   and is Telescope's documented/validated aligner.
// ----------------------------------------------------------------------------
process BOWTIE2_ALIGN_MULTI {
    tag   { meta.id }
    label 'process_high'
    publishDir "${params.outdir}/bowtie2_multi", mode: 'copy',
               saveAs: { fn -> fn.endsWith('.bam') ? fn : "logs/${fn}" }

    cpus     12
    memory   16.GB
    maxForks 1
    time     24.h

    input:
    tuple val(meta), path(reads)
    path  bt2_index_dir

    output:
    tuple val(meta), path("${meta.id}.multi.bam"), emit: bam
    path  "${meta.id}.bowtie2.log",                emit: log
    path  "versions.yml",                          emit: versions

    script:
    // paired vs single-end read arguments
    def readsArg = meta.single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    // bowtie2 index basename inside the passed index dir (built by BOWTIE2_BUILD)
    """
    IDX=\$(ls ${bt2_index_dir}/*.1.bt2* 2>/dev/null | head -1 | sed -E 's/\\.1\\.bt2l?\$//')
    bowtie2 \\
        -x \$IDX \\
        ${readsArg} \\
        -k 100 --very-sensitive-local --score-min L,0,1.6 \\
        -p ${task.cpus} \\
        --no-unal \\
        ${params.bowtie2_extra_args} \\
        2> ${meta.id}.bowtie2.log \\
        | samtools view -bS - > ${meta.id}.multi.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version 2>&1 | head -1 | sed 's/^.*version //')
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/^samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.multi.bam ${meta.id}.bowtie2.log
    echo '"${task.process}": {bowtie2: stub}' > versions.yml
    """
}

// ----------------------------------------------------------------------------
// PROCESS: BOWTIE2_BUILD
//   Build a bowtie2 index from the genome FASTA (arm64-native). Cached +
//   published so a cohort run builds it once. ~4 GB RAM for GRCh38.
// ----------------------------------------------------------------------------
process BOWTIE2_BUILD {
    tag   { fasta.baseName }
    label 'process_high'
    storeDir "${params.outdir}/bowtie2_index"

    cpus   12
    memory 16.GB
    time   12.h

    input:
    path fasta

    output:
    path "bowtie2_index", emit: index

    script:
    """
    mkdir -p bowtie2_index
    bowtie2-build --threads ${task.cpus} ${fasta} bowtie2_index/${fasta.baseName}
    """

    stub:
    """
    mkdir -p bowtie2_index
    touch bowtie2_index/${fasta.baseName}.1.bt2
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
    genome_fasta      // path : GRCh38 genome FASTA (bowtie2 index built from it)
    te_gtf_locus      // path : Telescope retro.hg38.v1 transcripts.gtf
    gene_gtf          // path : GENCODE genic GTF
    te_gtf_family     // path : TEtranscripts family TE GTF

    main:
    ch_versions = Channel.empty()

    // ---- Locus-level branch: multimap-permissive bowtie2 -> Telescope ----
    // (Telescope's documented aligner; the osx-arm64 STAR build ingests 0 reads.)
    BOWTIE2_BUILD ( genome_fasta )
    BOWTIE2_ALIGN_MULTI ( ch_reads, BOWTIE2_BUILD.out.index )
    ch_versions = ch_versions.mix( BOWTIE2_ALIGN_MULTI.out.versions.first() )

    TELESCOPE_ASSIGN ( BOWTIE2_ALIGN_MULTI.out.bam, te_gtf_locus )
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
