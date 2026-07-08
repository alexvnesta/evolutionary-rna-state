//
// Subworkflow: INTRON_RETENTION
// featureCounts-based intron-retention (IR) quantification, arm64-native.
// Consumes the rnaseq spine's coordinate-sorted STAR genome BAMs.
//
include { MAKE_INTRON_SAF                          } from './modules/make_intron_saf.nf'
include { FEATURECOUNTS_IR                         } from './modules/featurecounts_ir.nf'
include { FEATURECOUNTS_IR as FEATURECOUNTS_EXON   } from './modules/featurecounts_ir.nf'
include { COMPUTE_IR_RATIO                         } from './modules/compute_ir_ratio.nf'
include { MERGE_IR_MATRIX                          } from './modules/merge_ir_matrix.nf'

workflow INTRON_RETENTION {
    take:
    ch_bam      // channel: [ val(meta), path(bam), path(bai) ]
    gtf         // path:    GENCODE GTF (annotation)

    main:
    ch_versions = Channel.empty()

    // 1) derive pure-intronic + exonic SAF intervals once from the GTF
    MAKE_INTRON_SAF( gtf )
    ch_versions = ch_versions.mix( MAKE_INTRON_SAF.out.versions )

    // 2) featureCounts on intron intervals, per sample
    FEATURECOUNTS_IR(
        ch_bam,
        MAKE_INTRON_SAF.out.introns_saf,
        'intron'
    )
    ch_intron_counts = FEATURECOUNTS_IR.out.counts
    ch_versions = ch_versions.mix( FEATURECOUNTS_IR.out.versions.first() )

    // NOTE: a single process instance can only be invoked once per workflow.
    // We run the exon pass through a second, aliased include below.
    FEATURECOUNTS_EXON(
        ch_bam,
        MAKE_INTRON_SAF.out.exons_saf,
        'exon'
    )
    ch_exon_counts = FEATURECOUNTS_EXON.out.counts

    // 3) join intron + exon counts per sample -> compute IR ratio
    ch_pair = ch_intron_counts
        .map { meta, feat, f -> [ meta.id, meta, f ] }
        .join( ch_exon_counts.map { meta, feat, f -> [ meta.id, f ] } )
        .map { id, meta, intron_f, exon_f -> [ meta, intron_f, exon_f ] }

    COMPUTE_IR_RATIO( ch_pair, MAKE_INTRON_SAF.out.map )
    ch_versions = ch_versions.mix( COMPUTE_IR_RATIO.out.versions.first() )

    // 4) merge per-sample wide parquets into one cohort matrix (contract format)
    MERGE_IR_MATRIX(
        COMPUTE_IR_RATIO.out.long.map { meta, f -> f }.collect(),
        COMPUTE_IR_RATIO.out.summary.map { meta, f -> f }.collect()
    )
    ch_versions = ch_versions.mix( MERGE_IR_MATRIX.out.versions )

    emit:
    ir_long      = COMPUTE_IR_RATIO.out.long        // per-sample per-intron long tables
    ir_wide      = COMPUTE_IR_RATIO.out.wide        // per-sample wide parquet
    ir_matrix    = MERGE_IR_MATRIX.out.matrix       // cohort intron_retention.parquet
    ir_summary   = MERGE_IR_MATRIX.out.summary      // cohort per-sample summary
    versions     = ch_versions
}
