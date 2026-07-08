// RNA_EDITING subworkflow
// Inputs : ch_bam  = channel of [ meta, bam, bai ]  (STAR coordinate-sorted genome BAM)
//          fasta   = genome FASTA        (reference/GRCh38/GRCh38.primary_assembly.genome.fa)
//          fai     = FASTA index (.fai)
//          rmsk    = UCSC hg38 rmsk.txt.gz (for Alu intervals)
// Outputs : per-site A-to-I editing tables (JACUSA2) + per-sample & cohort AEI.
include { PREPARE_ALU_BED     } from '../modules/prepare_alu_bed.nf'
include { JACUSA2_CALL1       } from '../modules/jacusa2_call1.nf'
include { ALU_EDITING_INDEX   } from '../modules/alu_editing_index.nf'
include { MERGE_AEI           } from '../modules/alu_editing_index.nf'

workflow RNA_EDITING {
    take:
    ch_bam        // [ meta, bam, bai ]
    fasta         // path
    fai           // path
    rmsk          // path (UCSC rmsk.txt.gz)

    main:
    ch_versions = Channel.empty()

    // 1) Alu intervals (built once, reused for every sample).
    PREPARE_ALU_BED( rmsk )
    ch_alu = PREPARE_ALU_BED.out.alu_bed

    // 2) Per-site A-to-I calling (JACUSA2 call-1) -- optional, heavier.
    ch_sites = Channel.empty()
    if ( params.editing_call_sites ) {
        JACUSA2_CALL1( ch_bam, fasta, fai )
        ch_sites    = JACUSA2_CALL1.out.sites
        ch_versions = ch_versions.mix( JACUSA2_CALL1.out.versions.first() )
    }

    // 3) Alu Editing Index (robust cohort summary) -- always run.
    ALU_EDITING_INDEX( ch_bam, fasta, fai, ch_alu )
    ch_versions = ch_versions.mix( ALU_EDITING_INDEX.out.versions.first() )

    // 4) Merge per-sample AEI into a cohort table.
    MERGE_AEI( ALU_EDITING_INDEX.out.aei.map { meta, tsv -> tsv }.collect() )

    emit:
    sites      = ch_sites               // [ meta, editing_sites.tsv ]
    aei        = ALU_EDITING_INDEX.out.aei
    cohort_aei = MERGE_AEI.out.cohort   // cohort_aei.tsv
    versions   = ch_versions
}
