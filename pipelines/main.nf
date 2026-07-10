#!/usr/bin/env nextflow
/*
 * ============================================================================
 *  evolutionary-rna-state — UNIFIED non-reference feature workflow (Apple Silicon)
 * ============================================================================
 *  One top-level DSL2 workflow: a BAM channel (the HISAT2/nf-core spine output)
 *  fans out to the validated per-feature subworkflows and each emits a cohort
 *  matrix. This is a WRAPPER over arms already validated on osx-arm64, not a
 *  reimplementation.
 *
 *  Arms wired here (all consume the SAME [meta,bam,bai] spine BAM):
 *    - RNA_EDITING       (Alu Editing Index; unique-mapper)          rna_editing/
 *    - INTRON_RETENTION  (featureCounts intron/exon -> IR ratio)     intron_retention/
 *    - TE_ERV            (TEtranscripts family-level; +optional        te_erv/
 *                         Telescope locus-level when --te_locus)
 *  Splicing + fusion arms (rnasplice, rnafusion) are wrapped nf-core pipelines
 *  run as separate entries (see pipelines/scripts/); their per-sample outputs
 *  are merged into the final matrix by analysis/build_nonref_matrix.py. They are
 *  intentionally NOT fanned here because they re-alist from FASTQ rather than
 *  consuming the spine BAM.
 *
 *  Usage (Apple Silicon):
 *    nextflow run pipelines/main.nf -profile apple_silicon \
 *      --bam_glob 'results/rnaseq_cohort/hisat2/*.markdup.sorted.bam' \
 *      --outdir results/nonref_unified
 *
 *  Toggles:
 *    --te_locus false   (default) family-level TE only (~1 min/sample)
 *    --te_locus true    also run Telescope locus-level (bowtie2 -k100, ~2h/sample; cloud/Modal)
 *    --editing_call_sites false (default) AEI only; true adds JACUSA2 per-site (heavier)
 * ============================================================================
 */
// Param defaults that subworkflows read at include-time must be set BEFORE the includes.
params.editing_call_sites  = false     // AEI-only by default; true adds JACUSA2 per-site (heavier)
params.te_locus            = false     // heavy Telescope locus path OFF by default (CPU cost)

include { RNA_EDITING }      from './rna_editing/subworkflows/rna_editing.nf'
include { INTRON_RETENTION } from './intron_retention/subworkflow.nf'
include { TE_ERV }           from './te_erv/te_erv.nf'

// ---- params (arm64-safe defaults) ----
params.bam_glob            = null
params.input               = null      // optional samplesheet: sample,bam[,bai]
params.outdir              = 'results/nonref_unified'
params.fasta               = "${projectDir}/../reference/GRCh38/GRCh38.primary_assembly.genome.fa"
params.rmsk                = "${projectDir}/../reference/GRCh38/repeats/rmsk.hg38.txt.gz"
params.gene_gtf            = "${projectDir}/../reference/GRCh38/gencode.v46.primary_assembly.annotation.gtf"
params.te_gtf_family       = "${projectDir}/../reference/te/GRCh38_rmsk_TE.gtf"
params.te_gtf_locus        = "${projectDir}/../reference/te/retro.hg38.v1.transcripts.gtf"
params.help                = false

def helpMessage() {
    log.info """
    evolutionary-rna-state unified non-reference workflow
    Required: --bam_glob '<glob>'  OR  --input <samplesheet.csv>
    See header of pipelines/main.nf for options.
    """.stripIndent()
}

workflow {
    if ( params.help ) { helpMessage(); return }
    if ( !params.bam_glob && !params.input )
        exit 1, "ERROR: provide --bam_glob '<glob>' or --input <samplesheet.csv>"

    // ---- one BAM channel, fanned to every arm: [ meta, bam, bai ] ----
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
                def id  = bam.getBaseName().replaceAll(/\.(markdup\.)?(sorted|Aligned.*|genome)?$/, '')
                def bai = file("${bam}.bai")
                if ( !bai.exists() ) bai = file("${bam.toString().replaceAll(/\.bam$/, '.bai')}")
                tuple([id: id], bam, bai)
            }
    }

    // shared references
    fasta        = file(params.fasta,  checkIfExists: true)
    fai          = file("${params.fasta}.fai")
    rmsk         = file(params.rmsk,   checkIfExists: true)
    gene_gtf     = file(params.gene_gtf, checkIfExists: true)
    te_fam       = file(params.te_gtf_family, checkIfExists: true)

    // ---- ARM 1: RNA editing (AEI) ----
    RNA_EDITING( ch_bam, fasta, fai, rmsk )

    // ---- ARM 2: intron retention ----
    INTRON_RETENTION( ch_bam, gene_gtf )

    // ---- ARM 3: TE/ERV ----
    // TE_ERV take: (ch_reads, ch_bam[meta,bam], genome_fasta, te_gtf_locus, gene_gtf, te_gtf_family)
    // Family-level (TEtranscripts) consumes the spine BAM directly and is what runs by default.
    // Locus-level (Telescope) needs a bowtie2 -k100 re-align FROM READS. This unified entry consumes
    // spine BAMs, not FASTQs, so a real read source is NOT available here. Rather than feed Telescope
    // empty read lists (which would silently produce zero-alignment locus output), the locus path is
    // hard-gated OFF: use the standalone te_erv/ entry with --input fastqs (or the cloud/Modal runner)
    // when locus-level resolution is actually required.
    if ( params.te_locus )
        exit 1, "ERROR: --te_locus (Telescope locus-level) is not supported from the unified BAM-only entry — " +
                "it requires a FASTQ read source for bowtie2 -k100. Run pipelines/te_erv/ with --input <fastqs> " +
                "or the cloud runner instead. Family-level TE runs here by default."
    // Family-level path only reaches TETRANSCRIPTS_COUNT (uses gene_gtf + te_gtf_family). The locus-GTF
    // slot is only consumed by the Telescope branch, which is unreachable here (hard-gated above), so pass
    // the real params.te_gtf_locus for signature correctness rather than duplicating the family GTF.
    ch_bam_te = ch_bam.map { meta, bam, bai -> tuple(meta, bam) }
    te_locus_gtf = file(params.te_gtf_locus, checkIfExists: true)
    TE_ERV( Channel.empty(), ch_bam_te, fasta, te_locus_gtf, gene_gtf, te_fam )

    // ---- cohort matrices (each arm merges internally) ----
    RNA_EDITING.out.cohort_aei.view      { "[unified] AEI cohort matrix -> ${it}" }
    TE_ERV.out.family_matrix.view        { "[unified] TE family matrix  -> ${it}" }
}
