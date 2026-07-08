// JACUSA2 call-1 : single-BAM per-site variant/mismatch caller.
// For RNA editing we look for A>G (A-to-I) mismatches. JACUSA2 is a noarch Java
// jar (bioconda), so it runs natively on Apple Silicon under -profile conda.
//
// call-1 emits a BED-like table of positions with per-base counts and a score;
// we post-filter to A>G (sense) / T>C (antisense) with an editing frequency and
// coverage threshold, and optionally exclude dbSNP positions.
process JACUSA2_CALL1 {
    tag { meta.id }
    label 'process_editing'
    conda "${moduleDir}/../environment.yml"
    publishDir "${params.outdir}/rna_editing/sites", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("${meta.id}.jacusa2.out"),      emit: raw
    tuple val(meta), path("${meta.id}.editing_sites.tsv"), emit: sites
    path "versions.yml",                                   emit: versions

    script:
    // Keep only uniquely-mapped reads via MAPQ floor. HISAT2 (the arm64 spine)
    // marks unique mappers with MAPQ 60; STAR uses 255. Default params.editing_min_mapq=60
    // matches HISAT2 — set it to 255 if you feed STAR BAMs (amd64/Docker).
    // -a D,Y : filter distance-to-read-end & homopolymer artifacts (JACUSA2 recs).
    def snp_arg = params.editing_snp_bed ? "--snp-bed ${params.editing_snp_bed}" : ""
    """
    JAVA_MEM=\$(( ${task.memory.toGiga()} - 1 ))
    jacusa2 call-1 \\
        -r ${meta.id}.jacusa2.out \\
        -p ${task.cpus} \\
        -m ${params.editing_min_mapq} \\
        -q ${params.editing_min_baseq} \\
        -a D,Y \\
        -R ${fasta} \\
        ${bam}

    # Post-filter to A-to-I: A>G on '+' library reads, T>C on '-'.
    # JACUSA2 call-1 output columns: contig start end name score strand
    #   bases11(A,C,G,T) ...  We parse counts to compute editing frequency.
    filter_editing_sites.py \\
        --jacusa ${meta.id}.jacusa2.out \\
        --min-cov ${params.editing_min_cov} \\
        --min-edit-freq ${params.editing_min_freq} \\
        --min-edit-reads ${params.editing_min_edit_reads} \\
        ${snp_arg} \\
        --out ${meta.id}.editing_sites.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jacusa2: \$(jacusa2 --version 2>&1 | head -1 | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.jacusa2.out ${meta.id}.editing_sites.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jacusa2: stub
    END_VERSIONS
    """
}
