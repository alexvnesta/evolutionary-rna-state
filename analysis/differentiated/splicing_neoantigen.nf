//
// Subworkflow: SPLICING_NEOANTIGEN
// Splicing-derived neoantigen burden (SNAF algorithm, Li 2024
// doi:10.1126/scitranslmed.ade2886), wired to the shared antigen core.
//
// Per sample:  STAR SJ.out.tab (junction counts)  +  HLA-I alleles
//   -> tumor-specific NEOJUNCTIONS   (SNAF count gate: count-normal_mean>=t_min
//                                      AND normal_mean<n_max)
//   -> junction-spanning 8-11mer peptides  (SNAF 3-frame in-silico translation)
//   -> count_binders() via analysis/antigen_core/mhc_binding.py (MHCflurry)
//   -> splice_neoantigen_burden  (int, per run_accession+cohort)
//
// SNAF itself is NOT pip-installable on the arm64 dev box (it hard-pins
// tensorflow==2.3.0, no arm64/py3.11 wheel), so the SNAF ALGORITHM is
// reimplemented faithfully in analysis/differentiated/splicing_neoantigen.py
// and validated in test_splicing_neoantigen.py. On the Linux pipeline host
// SNAF *is* installable; if the pilot swaps in upstream SNAF, only the
// SPLICE_TRANSLATE process body changes — the channel contract is unchanged.
//
// This is a documented STUB other modules wire into (like hla_typing.nf). The
// stub{} blocks let `nextflow -stub-run` exercise the wiring without genome
// FASTA / MHCflurry present. Real translation + binding logic lives in the
// Python module and are unit-tested there.
//
// Channel contract (matches hla_typing / intron_retention subworkflows):
//   ch_sj       : [ val(meta), path(sj_out_tab) ]   meta = [id, cohort, single_end]
//   ch_hla      : path(hla_typing.parquet)          (from HLA_TYPING.out.hla_table)
//   ch_fasta    : path(genome.fa)                   (rnaseq spine GRCh38 FASTA)
//   ch_gtf      : path(annotation.gtf)              (GENCODE exon annotation)
//

include { SPLICE_CALL_NEOJUNCTIONS } from './modules/splice_call_neojunctions.nf'
include { SPLICE_TRANSLATE         } from './modules/splice_translate.nf'
include { SPLICE_SCORE_BURDEN      } from './modules/splice_score_burden.nf'
include { MERGE_SPLICE_BURDEN      } from './modules/merge_splice_burden.nf'

workflow SPLICING_NEOANTIGEN {
    take:
    ch_sj           // channel: [ val(meta), path(SJ.out.tab) ]
    ch_hla          // channel: path(hla_typing.parquet)
    ch_fasta        // channel: path(genome.fa)
    ch_gtf          // channel: path(annotation.gtf)

    main:
    ch_versions = Channel.empty()

    // 1) count gate: SJ.out.tab -> tumor-specific neojunctions (BED-like)
    SPLICE_CALL_NEOJUNCTIONS( ch_sj )
    ch_versions = ch_versions.mix( SPLICE_CALL_NEOJUNCTIONS.out.versions.first() )

    // 2) SNAF 3-frame translation across each neojunction -> candidate peptides
    SPLICE_TRANSLATE(
        SPLICE_CALL_NEOJUNCTIONS.out.neojunctions
            .combine( ch_fasta )
            .combine( ch_gtf )
    )
    ch_versions = ch_versions.mix( SPLICE_TRANSLATE.out.versions.first() )

    // 3) per-sample MHC-I binding via the SHARED engine -> burden int.
    //    HLA alleles are looked up from the cohort HLA table by run_accession.
    SPLICE_SCORE_BURDEN(
        SPLICE_TRANSLATE.out.peptides.combine( ch_hla )
    )
    ch_versions = ch_versions.mix( SPLICE_SCORE_BURDEN.out.versions.first() )

    // 4) merge per-sample burdens -> one tidy feature table (contract format)
    MERGE_SPLICE_BURDEN(
        SPLICE_SCORE_BURDEN.out.burden.map { meta, tsv -> tsv }.collect()
    )
    ch_versions = ch_versions.mix( MERGE_SPLICE_BURDEN.out.versions )

    emit:
    neojunctions = SPLICE_CALL_NEOJUNCTIONS.out.neojunctions  // per-sample
    peptides     = SPLICE_TRANSLATE.out.peptides              // per-sample fasta
    burden_table = MERGE_SPLICE_BURDEN.out.table              // splice_neoantigen_burden.parquet
    versions     = ch_versions
}

// ---------------------------------------------------------------------------
// Inline module definitions (kept in one file for a self-contained stub; in
// the pipeline these would live under ./modules/*.nf like the siblings). Each
// wraps the corresponding analysis/differentiated/splicing_neoantigen.py entry
// point through a thin CLI shim (splice_neoantigen_cli.py) so the algorithm is
// single-sourced with the unit-tested Python.
// ---------------------------------------------------------------------------

process SPLICE_CALL_NEOJUNCTIONS {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/splicing_neoantigen/neojunctions", mode: 'copy'
    conda "conda-forge::pandas conda-forge::pyarrow"

    input:
    tuple val(meta), path(sj)

    output:
    tuple val(meta), path("${meta.id}.neojunctions.tsv"), emit: neojunctions
    path "versions.yml",                                   emit: versions

    script:
    // SNAF gate params: t_min (min tumor-over-normal), n_max (max normal mean).
    // With a GTEx-style normal reference supplied via params.splice_normal_ref,
    // the CLI subtracts per-junction normal means; without it the gate is
    // count >= t_min (documented degenerate behaviour).
    def normal_ref = params.splice_normal_ref ? "--normal-ref ${params.splice_normal_ref}" : ''
    """
    splice_neoantigen_cli.py call-neojunctions \\
        --sj ${sj} \\
        --t-min ${params.splice_t_min ?: 20} \\
        --n-max ${params.splice_n_max ?: 3} \\
        ${normal_ref} \\
        --out ${meta.id}.neojunctions.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        splice_neoantigen: \$(splice_neoantigen_cli.py --version)
    END_VERSIONS
    """

    stub:
    """
    printf 'chrom\\tstart\\tend\\tstrand\\tcount\\n1\\t100\\t200\\t+\\t120\\n' > ${meta.id}.neojunctions.tsv
    echo '"${task.process}":' > versions.yml
    """
}

process SPLICE_TRANSLATE {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/splicing_neoantigen/peptides", mode: 'copy'
    conda "bioconda::pysam conda-forge::biopython conda-forge::pandas"

    input:
    tuple val(meta), path(neojunctions), path(fasta), path(gtf)

    output:
    tuple val(meta), path("${meta.id}.peptides.txt"), emit: peptides
    path "versions.yml",                              emit: versions

    script:
    """
    splice_neoantigen_cli.py translate \\
        --neojunctions ${neojunctions} \\
        --fasta ${fasta} \\
        --gtf ${gtf} \\
        --ks 8,9,10,11 \\
        --flank ${params.splice_flank ?: 60} \\
        --out ${meta.id}.peptides.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        splice_neoantigen: \$(splice_neoantigen_cli.py --version)
    END_VERSIONS
    """

    stub:
    """
    printf 'GILGFVFTL\\nNLVPMVATV\\n' > ${meta.id}.peptides.txt
    echo '"${task.process}":' > versions.yml
    """
}

process SPLICE_SCORE_BURDEN {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/splicing_neoantigen/burden", mode: 'copy'
    // MHCflurry engine (shared antigen_core). MHCflurry is pip-only; the pilot
    // supplies it via the process container/env — pinned here for provenance.
    conda "bioconda::mhcflurry=2.* conda-forge::pandas conda-forge::pyarrow"

    input:
    tuple val(meta), path(peptides), path(hla_table)

    output:
    tuple val(meta), path("${meta.id}.burden.tsv"), emit: burden
    path "versions.yml",                            emit: versions

    script:
    """
    export MHCFLURRY_DATA_DIR=${params.mhcflurry_data_dir ?: "${projectDir}/reference/mhcflurry_models"}
    splice_neoantigen_cli.py score-burden \\
        --peptides ${peptides} \\
        --hla-table ${hla_table} \\
        --run-accession ${meta.id} \\
        --cohort ${meta.cohort} \\
        --rank ${params.splice_binder_rank ?: 2.0} \\
        --out ${meta.id}.burden.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcflurry: \$(python -c 'import mhcflurry; print(mhcflurry.__version__)')
    END_VERSIONS
    """

    stub:
    """
    printf 'run_accession\\tcohort\\tsplice_neoantigen_burden\\n${meta.id}\\t${meta.cohort}\\t2\\n' > ${meta.id}.burden.tsv
    echo '"${task.process}":' > versions.yml
    """
}

process MERGE_SPLICE_BURDEN {
    label 'process_low'
    publishDir "${params.outdir}/features", mode: 'copy'
    conda "conda-forge::pandas conda-forge::pyarrow"

    input:
    path burdens    // collected per-sample *.burden.tsv

    output:
    path "splice_neoantigen_burden.parquet", emit: table
    path "versions.yml",                     emit: versions

    script:
    """
    splice_neoantigen_cli.py merge-burden \\
        --burdens ${burdens} \\
        --out splice_neoantigen_burden.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    """
    python -c "import pandas as pd; pd.DataFrame({'run_accession':['ERRTEST01'],'cohort':['gide2019'],'splice_neoantigen_burden':[2]}).to_parquet('splice_neoantigen_burden.parquet')"
    echo '"${task.process}":' > versions.yml
    """
}
