// Minimal harness to -stub-run the IR_NEOANTIGEN subworkflow wiring.
nextflow.enable.dsl=2
include { IR_NEOANTIGEN } from '../intron_retention.nf'

workflow {
    IR_NEOANTIGEN(
        file(params.ir_matrix),
        file(params.intron_saf),
        file(params.genome),
        file(params.hla_table)
    )
    IR_NEOANTIGEN.out.features.view { "features -> ${it}" }
}
