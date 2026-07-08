#!/usr/bin/env nextflow
/*
 * RNA-editing subworkflow entrypoint (A-to-I sites + Alu Editing Index).
 *
 * Consumes the coordinate-sorted STAR genome BAMs produced by the nf-core/rnaseq
 * (hisat2 aligner) spine. Two input modes:
 *   --bam_glob '<...>/*.markdup.sorted.bam'   (auto-pairs .bai)
 *   --input samplesheet.csv                    (columns: sample,bam,bai)
 *
 * Example:
 *   source pipelines/env.sh
 *   nextflow run pipelines/rna_editing/main.nf \
 *     -profile conda -c pipelines/conf/mac_arm64.config \
 *     -c pipelines/rna_editing/conf/editing.config \
 *     --bam_glob 'results/rnaseq/hisat2/*.markdup.sorted.bam' \
 *     --fasta reference/GRCh38/GRCh38.primary_assembly.genome.fa \
 *     --rmsk  reference/GRCh38/repeats/rmsk.hg38.txt.gz \
 *     --outdir results/rna_editing
 */
nextflow.enable.dsl = 2

include { RNA_EDITING } from './subworkflows/rna_editing.nf'

// ---- parameters (defaults; override on CLI or in editing.config) ----
params.help         = false
params.input        = null      // samplesheet CSV: sample,bam,bai
params.bam_glob     = null      // alternative to --input
params.fasta        = null
params.rmsk         = null      // UCSC hg38 rmsk.txt.gz
params.outdir       = 'results/rna_editing'

// editing thresholds
params.editing_call_sites        = true    // run JACUSA2 per-site caller
params.editing_min_mapq          = 255     // STAR unique-mapper MAPQ
params.editing_min_baseq         = 25
params.editing_min_cov           = 10
params.editing_min_freq          = 0.10
params.editing_min_edit_reads    = 3
params.editing_snp_bed           = null    // bgzipped+tabixed dbSNP BED (optional)
params.editing_alu_standard_only = true    // Alu on chr1..22,X,Y,M only

def helpMessage() {
    log.info """
    RNA-editing subworkflow
    -----------------------
    Required:
      --fasta        genome FASTA (GRCh38 primary assembly, chr-prefixed)
      --rmsk         UCSC hg38 rmsk.txt.gz  (Alu source; see fetch_alu.sh)
      one of:
      --bam_glob     glob of STAR coord-sorted BAMs (auto-pairs .bai)
      --input        samplesheet CSV (columns: sample,bam,bai)

    Key options (defaults):
      --editing_call_sites      ${params.editing_call_sites}
      --editing_min_mapq        ${params.editing_min_mapq}
      --editing_min_baseq       ${params.editing_min_baseq}
      --editing_min_cov         ${params.editing_min_cov}
      --editing_min_freq        ${params.editing_min_freq}
      --editing_snp_bed         (optional dbSNP BED to mask known SNPs)
      --outdir                  ${params.outdir}
    """.stripIndent()
}

workflow {
    if ( params.help ) { helpMessage(); return }
    if ( !params.fasta ) exit 1, "ERROR: --fasta is required"
    if ( !params.rmsk )  exit 1, "ERROR: --rmsk is required (UCSC hg38 rmsk.txt.gz; run pipelines/rna_editing/fetch_alu.sh)"
    if ( !params.input && !params.bam_glob ) exit 1, "ERROR: provide --bam_glob or --input"

    // Build [ meta, bam, bai ] channel
    if ( params.input ) {
        ch_bam = Channel.fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def bai = row.bai ?: "${row.bam}.bai"
                tuple([id: row.sample], file(row.bam, checkIfExists: true), file(bai, checkIfExists: true))
            }
    } else {
        ch_bam = Channel.fromPath(params.bam_glob, checkIfExists: true)
            .map { bam ->
                def id = bam.getBaseName().replaceAll(/\.(markdup\.)?(sorted|Aligned.*|genome)?$/, '')
                def bai = file("${bam}.bai")
                if ( !bai.exists() ) bai = file("${bam.toString().replaceAll(/\.bam$/, '.bai')}")
                tuple([id: id], bam, bai)
            }
    }

    fasta = file(params.fasta, checkIfExists: true)
    fai   = file("${params.fasta}.fai")   // built by fetch_alu.sh / samtools faidx if missing
    rmsk  = file(params.rmsk, checkIfExists: true)

    RNA_EDITING( ch_bam, fasta, fai, rmsk )
}
