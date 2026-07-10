//
// Module: MERGE_BCR_FEATURES
// Concatenates the per-sample BCR feature-row CSVs into one cohort
// bcr_features.parquet (keyed run_accession, cohort). Stub-runnable.
//
process MERGE_BCR_FEATURES {
    label 'process_single'
    publishDir "${params.outdir}/features", mode: 'copy'
    conda "conda-forge::pandas conda-forge::pyarrow"

    input:
    path rows

    output:
    path "bcr_features.parquet", emit: features
    path "versions.yml",         emit: versions

    script:
    """
    python3 - <<'PY'
    import glob, pandas as pd
    frames = [pd.read_csv(f) for f in sorted(glob.glob("*_bcr_row.csv"))]
    out = pd.concat(frames, ignore_index=True)
    out.to_parquet("bcr_features.parquet", index=False)
    print(f"merged {len(out)} BCR feature rows")
    PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch bcr_features.parquet
    echo '"${task.process}":' > versions.yml
    """
}
