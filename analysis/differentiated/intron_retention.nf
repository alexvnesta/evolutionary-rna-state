//
// Subworkflow: IR_NEOANTIGEN
// DIFFERENTIATED BUCKET — retained-intron load + cryptic-ORF neoantigen burden.
//
// Sits DOWNSTREAM of the pipeline session's intron-retention quantification
// (pipelines/intron_retention → intron_retention.parquet + introns.saf +
// intron2gene.tsv) and the antigen-core HLA typing (hla_typing.parquet). It
// derives, per sample, two NAMED interpretable features:
//     retained_intron_load   (int)   — # introns with IR ratio >= threshold
//                                       (+ depth-robust fraction & cohort-z)
//     ir_neoantigen_burden   (int)   — MHCflurry binder count over peptides
//                                       translated from retained introns
//                                       (junction-spanning + intronic ORFs),
//                                       via the SHARED antigen_core engine.
//
// The feature logic lives (and is unit-tested on synthetic inputs) in
// analysis/differentiated/intron_retention.py; this .nf only wires the pipeline
// outputs into the CLI wrapper analysis/differentiated/bin/ir_antigen_features.py.
//
// This is a documented, correct STUB (mirrors analysis/antigen_core/hla_typing.nf):
// the `stub:` block lets `nextflow -stub-run` exercise the wiring without the
// genome FASTA / MHCflurry models present; the full run is the pilot.
//
// Channel / input contract:
//   ir_matrix   : path  intron_retention.parquet (contract tidy-wide)
//   intron_saf  : path  introns.saf              (make_intron_saf.py output)
//   genome      : path  GRCh38 primary_assembly FASTA (+ .fai)
//   hla_table   : path  hla_typing.parquet       (antigen_core HLA_TYPING output)
//

include { IR_ANTIGEN_FEATURES } from './modules/ir_antigen_features.nf'

workflow IR_NEOANTIGEN {
    take:
    ir_matrix       // path: intron_retention.parquet
    intron_saf      // path: introns.saf
    genome          // path: GRCh38 genome FASTA (indexed)
    hla_table       // path: hla_typing.parquet

    main:
    ch_versions = Channel.empty()

    IR_ANTIGEN_FEATURES( ir_matrix, intron_saf, genome, hla_table )
    ch_versions = ch_versions.mix( IR_ANTIGEN_FEATURES.out.versions )

    emit:
    features  = IR_ANTIGEN_FEATURES.out.features   // ir_antigen_features.parquet
    versions  = ch_versions
}
