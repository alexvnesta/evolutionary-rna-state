#!/usr/bin/env nextflow
//
// Standalone entrypoint for the INTRON_RETENTION subworkflow.
//
// Consumes coordinate-sorted STAR genome BAMs (produced by the nf-core/rnaseq
// hisat2 spine, or star_salmon on amd64) and emits per-sample + cohort intron-retention (IR) ratio
// matrices. featureCounts-based, fully arm64-native (see README.md).
//
// Two ways to supply BAMs:
//   --input   samplesheet CSV: run_accession,cohort,bam,bai  (recommended)
//   --bam_glob '/path/to/*.bam'  (bai discovered next to each bam)
//
nextflow.enable.dsl = 2

include { INTRON_RETENTION } from './subworkflow.nf'

// ---- parameters (override on CLI or via -params-file) ----
params.help                    = false
params.input                   = null
params.bam_glob                = null
params.gtf                     = "${projectDir}/../../reference/GRCh38/gencode.v46.primary_assembly.annotation.gtf"
params.outdir                  = "${projectDir}/../../results/intron_retention"
params.ir_min_intron_len       = 50      // drop pure-intronic sub-intervals shorter than this
params.ir_min_gene_exon_count  = 20      // require host-gene exonic signal to trust an IR value
params.ir_high_threshold       = 0.1     // IR ratio above which an intron is 'retained' (summary)
params.ir_strandedness         = 2       // featureCounts -s : 0 unstranded, 1 fwd, 2 rev (Illumina dUTP)

def helpMessage() {
    log.info """
    ========================================================================
    INTRON_RETENTION  --  featureCounts-based IR ratio (arm64-native)
    ========================================================================
    Usage:
      source pipelines/env.sh
      nextflow run pipelines/intron_retention/main.nf \\
          -profile conda -c pipelines/conf/mac_arm64.config \\
          --input samplesheet.csv \\
          --gtf   reference/GRCh38/gencode.v46.primary_assembly.annotation.gtf \\
          --outdir results/intron_retention

    Samplesheet CSV columns: run_accession,cohort,bam,bai
    Key params:
      --ir_strandedness        featureCounts -s (0/1/2)  [${params.ir_strandedness}]
      --ir_min_intron_len      min pure-intron length bp [${params.ir_min_intron_len}]
      --ir_min_gene_exon_count host-gene exonic read floor [${params.ir_min_gene_exon_count}]
      --ir_high_threshold      'retained' IR cutoff for summary [${params.ir_high_threshold}]
    """.stripIndent()
}

workflow {
    if (params.help) { helpMessage(); return }

    // ---- build the BAM channel: [ meta, bam, bai ] ----
    if (params.input) {
        ch_bam = Channel
            .fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def meta = [ id: row.run_accession,
                             cohort: (row.cohort ?: 'NA'),
                             single_end: (row.single_end?.toString()?.toLowerCase() in ['true','1','yes']) ]
                def bam = file(row.bam, checkIfExists: true)
                def bai = row.bai ? file(row.bai, checkIfExists: true) : file("${row.bam}.bai")
                [ meta, bam, bai ]
            }
    } else if (params.bam_glob) {
        ch_bam = Channel
            .fromPath(params.bam_glob, checkIfExists: true)
            .map { bam ->
                def id = bam.simpleName.replaceAll(/\.Aligned.*$/, '')
                def bai = file("${bam}.bai")
                [ [ id: id, cohort: 'NA', single_end: false ], bam, bai ]
            }
    } else {
        error "Provide --input <samplesheet.csv> or --bam_glob '<glob>'. See --help."
    }

    ch_gtf = file(params.gtf, checkIfExists: true)

    INTRON_RETENTION( ch_bam, ch_gtf )
}
