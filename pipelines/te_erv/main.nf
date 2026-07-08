#!/usr/bin/env nextflow
/*
 * Standalone entrypoint for the TE / ERV activation subworkflow.
 *
 * Runs the dedicated multimap-aware STAR -> Telescope (locus level) and
 * TEcount (family level) branches. Consumes:
 *   --input      samplesheet CSV: sample,fastq_1,fastq_2,bam,strandedness
 *   --star_index prebuilt STAR index dir (shared repo index)
 *   --gene_gtf   GENCODE genic annotation
 *   --te_gtf_locus   Telescope retro.hg38.v1 transcripts.gtf   (locus level)
 *   --te_gtf_family  TEtranscripts family TE GTF               (family level)
 *
 * Example:
 *   nextflow run pipelines/te_erv/main.nf \
 *     -profile conda -c pipelines/conf/mac_arm64.config \
 *     --input   pipelines/te_erv/assets/samplesheet_test.csv \
 *     --outdir  results/te_erv \
 *     --star_index   reference/GRCh38/star_index \
 *     --gene_gtf     reference/GRCh38/gencode.v46.primary_assembly.annotation.gtf \
 *     --te_gtf_locus  reference/te/retro.hg38.v1.transcripts.gtf \
 *     --te_gtf_family reference/te/GRCh38_GENCODE_rmsk_TE.gtf
 */

nextflow.enable.dsl = 2

include { TE_ERV } from './te_erv.nf'

// ---- parameter defaults ----
params.input               = null
params.outdir              = 'results/te_erv'
params.star_index          = null
params.gene_gtf            = null
params.te_gtf_locus        = null
params.te_gtf_family       = null
params.star_extra_args         = ''
params.telescope_extra_args    = ''
params.tetranscripts_extra_args = ''
params.help                = false

def helpMessage() {
    log.info """
    TE / ERV activation subworkflow
    ===============================
    Required:
      --input            samplesheet CSV (cols: sample,fastq_1,fastq_2,bam,strandedness)
      --star_index       prebuilt STAR index directory
      --gene_gtf         GENCODE genic GTF
      --te_gtf_locus     Telescope locus-level TE/ERV GTF (retro.hg38.v1 transcripts.gtf)
      --te_gtf_family    TEtranscripts family-level TE GTF
    Optional:
      --outdir           output directory        [${params.outdir}]
      --star_extra_args / --telescope_extra_args / --tetranscripts_extra_args
    Run:
      nextflow run pipelines/te_erv/main.nf -profile conda \\
        -c pipelines/conf/mac_arm64.config --input <csv> --star_index <dir> ...
    """.stripIndent()
}

workflow {
    if (params.help) { helpMessage(); return }

    // ---- validate required params ----
    def missing = []
    ['input','star_index','gene_gtf','te_gtf_locus','te_gtf_family'].each { p ->
        if (!params[p]) missing << "--${p}"
    }
    if (missing) {
        log.error "Missing required parameter(s): ${missing.join(', ')}\nRun with --help for usage."
        System.exit(1)
    }

    // ---- parse samplesheet ----
    // Two channels derived from one sheet:
    //   ch_reads : rows with fastq_1 (Telescope needs its own multimap STAR run)
    //   ch_bam   : rows with a bam column (TEcount reuses the rnaseq spine BAM)
    ch_rows = Channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)

    ch_reads = ch_rows
        .filter { row -> row.fastq_1 && row.fastq_1.trim() }
        .map { row ->
            def meta = [ id: row.sample,
                         single_end: !(row.fastq_2 && row.fastq_2.trim()),
                         strandedness: row.strandedness ?: 'no' ]
            def reads = meta.single_end ? [ file(row.fastq_1, checkIfExists: true) ]
                                        : [ file(row.fastq_1, checkIfExists: true),
                                            file(row.fastq_2, checkIfExists: true) ]
            tuple(meta, reads)
        }

    ch_bam = ch_rows
        .filter { row -> row.bam && row.bam.trim() }
        .map { row ->
            def meta = [ id: row.sample, strandedness: row.strandedness ?: 'no' ]
            tuple(meta, file(row.bam, checkIfExists: true))
        }

    TE_ERV(
        ch_reads,
        ch_bam,
        file(params.star_index,    checkIfExists: true),
        file(params.gene_gtf,      checkIfExists: true),
        file(params.te_gtf_locus,  checkIfExists: true),
        file(params.te_gtf_family, checkIfExists: true)
    )
}
