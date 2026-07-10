//
// Subworkflow: BCR_SHM
// B-cell-receptor repertoire, somatic hypermutation (SHM), isotype
// class-switch and clonality from bulk RNA-seq, via TRUST4.
//
// Emits, per sample (keyed run_accession, cohort):
//   bcr_shm_rate, bcr_igg_fraction, bcr_switched_fraction, bcr_clonality,
//   bcr_n_clonotypes, bcr_n_reads   (+ audit cols _shm_source, _shm_n)
//
// This is the repertoire-level read-out of B-cell affinity maturation — the
// same biology as the population-level baseline TLS/B-cell expression score
// (analysis/baseline/tls_bcell_scores.py), measured on reconstructed Ig
// sequences rather than aggregate expression.
//
// DEPTH REQUIREMENT (batch robustness): BCR contig yield, and therefore SHM /
// clonality stability, scale with sequencing depth and B-cell content. This
// arm streams a LARGER read subsample than the gene-TPM arm (params.trust4_nreads,
// default 15e6 pairs) and records bcr_n_clonotypes / bcr_n_reads on every row so
// per-sample reliability is auditable. Samples below params.min_clonotypes
// (default 3) IGH contigs return NaN features — never imputed. Hold the
// subsample size CONSTANT across a comparison (same rationale as the
// fixed-caller requirement in fusion_antigen.nf).
//
// This is a documented, correct wiring that mirrors the sibling subworkflows
// (fusion_antigen.nf / hla_typing.nf). The streaming + TRUST4 run is done by
// pipelines/bcr_repertoire/run_trust4_pilot.sh; the feature logic lives in
// analysis/differentiated/bcr_shm.py and is unit-tested in test_bcr_shm.py.
// The stub{} blocks let `nextflow -stub-run` exercise the wiring without
// TRUST4 / network present.
//
// Channel contract (matches the sibling subworkflows):
//   ch_reads : [ val(meta), path(r1), path(r2) ]   meta = [id, cohort]
//              OR, for the stream-from-ENA pilot arm:
//   ch_manifest : path(manifest_csv)   (run_accession, cohort, fastq_ftp)
//

include { TRUST4_ASSEMBLE      } from './modules/trust4_assemble.nf'
include { BCR_SHM_FEATURES     } from './modules/bcr_shm_features.nf'
include { MERGE_BCR_FEATURES   } from './modules/merge_bcr_features.nf'

workflow BCR_SHM {
    take:
    ch_reads        // channel: [ val(meta), path(r1), path(r2) ]

    main:
    ch_versions = Channel.empty()

    // 1) per-sample: TRUST4 reconstruction -> _cdr3.out/_report.tsv/_airr.tsv
    TRUST4_ASSEMBLE( ch_reads )
    ch_versions = ch_versions.mix( TRUST4_ASSEMBLE.out.versions.first() )

    // 2) per-sample: parse TRUST4 output -> one feature row (bcr_shm.py)
    BCR_SHM_FEATURES( TRUST4_ASSEMBLE.out.trust4_dir )
    ch_versions = ch_versions.mix( BCR_SHM_FEATURES.out.versions.first() )

    // 3) merge per-sample rows -> one cohort bcr_features table
    MERGE_BCR_FEATURES(
        BCR_SHM_FEATURES.out.row.map { meta, csv -> csv }.collect()
    )
    ch_versions = ch_versions.mix( MERGE_BCR_FEATURES.out.versions )

    emit:
    features = MERGE_BCR_FEATURES.out.features   // path: bcr_features.parquet
    versions = ch_versions
}
