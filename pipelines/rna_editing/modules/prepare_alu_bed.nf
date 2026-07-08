// Build the Alu-only BED6 from a UCSC hg38 RepeatMasker table (rmsk.txt.gz).
// The rmsk file is a parameter (params.editing_rmsk); if absent, the pipeline
// stops with a clear message pointing at fetch_alu.sh. This process does NOT
// download anything itself (keeps the run reproducible / offline-friendly).
process PREPARE_ALU_BED {
    label 'process_light'
    conda "${moduleDir}/../environment.yml"
    publishDir "${params.outdir}/rna_editing/refs", mode: 'copy'

    input:
    path rmsk

    output:
    path "alu.hg38.bed6", emit: alu_bed

    script:
    def std = params.editing_alu_standard_only ? "--keep-standard" : ""
    """
    make_alu_bed.py --rmsk ${rmsk} --out alu.hg38.bed6 ${std}
    """

    stub:
    """
    printf 'chr1\\t100000\\t100300\\tAluY\\t2000\\t+\\n' > alu.hg38.bed6
    """
}
