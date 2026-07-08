#!/usr/bin/env python
"""
pilot_ingest.py — one-shot orchestrator: run every feature module on the
pipeline session's REAL per-sample matrices as they land, join into one
per-sample feature matrix, then run batch-robustness + LOCO/LOPO evaluation.

Reads the handoff contract's directory (results/features/) and DEGRADES
GRACEFULLY: any phenotype matrix that is absent is skipped and its features
left as NA; the evaluation harness scores whatever is present. This is the
'auto-fill the PENDING slots' path referenced in FEATURE_EVALUATION_REPORT.md.

Usage:
    MHCFLURRY_DATA_DIR=reference/mhcflurry_models \
    python analysis/pilot_ingest.py --features-dir results/features \
        --hla results/features/hla_typing.parquet \
        --clinical results/analysis_frame.parquet --outdir results/pilot

Nothing here fabricates values: if an input matrix is missing, the feature is
NA, not imputed.
"""
import os, sys, argparse, json
import pandas as pd

os.environ.setdefault("MHCFLURRY_DATA_DIR", os.path.abspath("reference/mhcflurry_models"))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

KEY = ["run_accession", "cohort"]

# phenotype matrix -> (filename in results/features, builder callable name, module)
INPUTS = {
    "quant_gene_tpm":     "quant_gene_tpm.parquet",
    "splicing_junctions": "splicing_junctions.parquet",
    "intron_retention":   "intron_retention.parquet",
    "rna_editing_aei":    "rna_editing_aei.parquet",
    "rna_editing_sites":  "rna_editing_sites.parquet",
    "te_locus":           "te_locus.parquet",
    "fusion_calls":       "fusion_calls.parquet",
    "variants_maf":       "variants.maf",
}

def _exists(features_dir, fn):
    p = os.path.join(features_dir, fn)
    return p if os.path.exists(p) else None

def build_all(features_dir, hla_path, clinical_path, cohort_map=None):
    """cohort_map: optional dict run_accession -> cohort, used for within-cohort
    GEP harmonization when the TPM matrix carries no cohort column. Without it,
    GEP within-batch z-scoring collapses to NaN."""
    frames = []
    log = {}

    # HLA table (per-sample 6 alleles + heterozygosity) — required for antigen modules
    hla = pd.read_parquet(hla_path) if hla_path and os.path.exists(hla_path) else None
    if hla is not None:
        frames.append(hla[[c for c in hla.columns if c in KEY or c == "HLA_I_heterozygous"]])
        from analysis.differentiated.intron_retention import hla_map_from_table
        hla_by_sample = hla_map_from_table(hla)
        log["hla"] = f"{len(hla)} samples typed"
    else:
        hla_by_sample = {}
        log["hla"] = "ABSENT — antigen-burden features will be NA"

    # ---- BASELINE: GEP (needs gene TPM) ----
    tpm_p = _exists(features_dir, INPUTS["quant_gene_tpm"])
    if tpm_p:
        from analysis.baseline.gep_scores import score_all
        from analysis.pilot_gep import to_symbol_gene_matrix
        raw = pd.read_parquet(tpm_p)
        # Normalize to genes(symbol) x samples regardless of input orientation
        # / gene-id namespace (real pilot TPM is samples x ENSG columns).
        mat, sm = to_symbol_gene_matrix(raw)
        # GEP uses within-cohort harmonization; supply cohort if the matrix lacks
        # it (else within-batch z-scoring yields NaN).
        if cohort_map is not None and sm["cohort"].isna().all():
            sm = sm.copy()
            sm["cohort"] = sm["run_accession"].map(cohort_map)
        if sm["cohort"].isna().all():
            # no cohort available at all — score without harmonization (single batch)
            sm = sm.copy(); sm["cohort"] = "ALL"
            log["gep_note"] = "no cohort map — scored as a single batch (no harmonization)"
        from analysis.baseline.gep_scores import GEP_TCELL_INFLAMED_GENES
        n_sig = sum(g in mat.index for g in GEP_TCELL_INFLAMED_GENES)
        frames.append(score_all(mat, sm))
        log["gep"] = f"scored {mat.shape[1]} samples ({n_sig}/{len(GEP_TCELL_INFLAMED_GENES)} GEP genes present)"
    else:
        log["gep"] = "quant_gene_tpm.parquet ABSENT — GEP features NA"

    # ---- DIFFERENTIATED (expression-derived): RNA-regulator activity ----
    # Computable NOW from the same gene TPM: splicing-factor / broad-RBP / ADAR
    # activity (rbp-activity-scorer). Named, interpretable proxies for the
    # differentiated RNA state; NOT the de-novo antigen burden. Guarded on the
    # skill being importable so a standalone run without it still degrades.
    if tpm_p:
        try:
            from analysis.regulator_activity import (
                build_regulator_activity, gene_map_from_matrix, ACTIVITY_COLUMNS, regulator_sets)
            raw2 = pd.read_parquet(tpm_p)
            if "gene_name" not in raw2.columns:
                log["regulator_activity"] = "TPM lacks gene_name column — skipped"
            else:
                reg_syms = sorted({s for v in regulator_sets().values() for s in v})
                gmap = gene_map_from_matrix(raw2, reg_syms)
                scols = [c for c in raw2.columns if c != "gene_name" and not str(c).startswith("ENSG")]
                ens_idx = [str(i).split(".")[0] for i in raw2.index]
                tdf = raw2.assign(_e=ens_idx).set_index("_e")[scols].T
                ra = build_regulator_activity(tdf, gene_symbol_index=gmap)  # run_accession + 3 cols
                frames.append(ra[["run_accession"] + [c for c in ACTIVITY_COLUMNS if c in ra.columns]])
                log["regulator_activity"] = f"scored {len(ra)} samples ({len(gmap)} regulator genes)"
        except Exception as e:
            log["regulator_activity"] = f"skipped: {e}"

    # ---- DIFFERENTIATED: intron-retention load (needs IR ratios; no genome) ----
    ir_p = _exists(features_dir, INPUTS["intron_retention"])
    if ir_p:
        from analysis.differentiated.intron_retention import compute_retained_intron_load
        frames.append(compute_retained_intron_load(pd.read_parquet(ir_p)))
        log["intron_retention"] = "retained_intron_load computed"
    else:
        log["intron_retention"] = "intron_retention.parquet ABSENT"

    # ---- DIFFERENTIATED: RNA editing AEI (needs aei table; no genome) ----
    aei_p = _exists(features_dir, INPUTS["rna_editing_aei"])
    if aei_p:
        from analysis.differentiated.rna_editing import compute_alu_editing_index
        aei = pd.read_parquet(aei_p)
        frames.append(compute_alu_editing_index(aei, cohort=None))
        log["rna_editing_aei"] = "alu_editing_index computed"
    else:
        log["rna_editing_aei"] = "rna_editing_aei.parquet ABSENT"

    # ---- Antigen-burden modules (need HLA + sequence refs) ----
    # These require genome/annotation refs staged on the pilot host; each module
    # exposes a build_*_features(...) entrypoint. Wired but guarded on inputs.
    for name, note in [
        ("te_locus", "analysis.differentiated.te_antigen.build_te_antigen_table (needs genome FASTA + RepeatMasker GTF)"),
        ("splicing_junctions", "analysis.differentiated.splicing_neoantigen.build_feature_table (needs genome FASTA + exon index)"),
        ("fusion_calls", "analysis.differentiated.fusion_antigen.build_fusion_feature_table (needs Arriba/STAR-Fusion TSV)"),
        ("variants_maf", "analysis.baseline.snv_indel_neoantigen.variants_from_maf + snv_indel_neoantigen_burden (needs proteome FASTA)"),
    ]:
        p = _exists(features_dir, INPUTS[name])
        log[name] = (f"present ({note}) — run module with staged refs" if p
                     else f"{INPUTS[name]} ABSENT — {name} burden NA")

    # join everything on the contract key. Some builders emit run_accession only
    # (cohort is attached from the crosswalk downstream); merge on whichever key
    # columns are shared so a cohort-less frame still joins cleanly.
    M = None
    for f in frames:
        f = f.loc[:, ~f.columns.duplicated()]
        if M is None:
            M = f
        else:
            on = [k for k in KEY if k in M.columns and k in f.columns] or ["run_accession"]
            M = M.merge(f, on=on, how="outer")
    return M, log

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--features-dir", default="results/features")
    ap.add_argument("--hla", default="results/features/hla_typing.parquet")
    ap.add_argument("--clinical", default="results/analysis_frame.parquet")
    ap.add_argument("--run-catalog", dest="run_catalog", default="results/features/run_catalog.csv",
                    help="pipeline run_catalog.csv (run_accession,cohort,patient_id,timepoint) — "
                         "preferred crosswalk source; covers all cohorts")
    ap.add_argument("--outdir", default="results/pilot")
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)

    # cohort map (run_accession -> cohort) for GEP harmonization, from run_catalog
    cohort_map = None
    if a.run_catalog and os.path.exists(a.run_catalog):
        _cat = pd.read_csv(a.run_catalog)
        if {"run_accession", "cohort"}.issubset(_cat.columns):
            cohort_map = dict(zip(_cat["run_accession"], _cat["cohort"]))

    M, log = build_all(a.features_dir, a.hla, a.clinical, cohort_map=cohort_map)
    print("INGEST LOG:"); print(json.dumps(log, indent=1))
    if M is None or M.empty:
        print("No feature matrices present yet — nothing to score. Re-run when results/features/ populates.")
        return 0
    out_mat = os.path.join(a.outdir, "pilot_feature_matrix.parquet")
    M.to_parquet(out_mat)
    print(f"wrote {out_mat}: {M.shape[0]} samples x {M.shape[1]} cols")

    # attach labels + run evaluation on whatever landed.
    # The pipeline keys on run_accession; the clinical frame keys on study
    # sampleId — bridge them with the validated ENA crosswalk, never a naive
    # run_accession merge (which silently yields all-NA labels).
    clin = pd.read_parquet(a.clinical)
    from analysis.pilot_crosswalk import (
        fetch_ena_metadata, build_crosswalk, build_crosswalk_from_catalog, attach_labels)
    if "run_accession" in clin.columns and clin["run_accession"].isin(M["run_accession"]).any():
        ev = M.merge(clin, on=KEY, how="left")
    elif a.run_catalog and os.path.exists(a.run_catalog):
        # Preferred: the pipeline's run_catalog covers all cohorts (incl. hugo)
        cat = pd.read_csv(a.run_catalog)
        cat = cat[cat["run_accession"].isin(M["run_accession"])]
        xwalk = build_crosswalk_from_catalog(cat, clin)
        ev = attach_labels(M, xwalk, clin)
        print(f"crosswalk (run_catalog): mapped {len(xwalk)}/{len(M)} runs; "
              f"RESPONDER labels attached to {ev['RESPONDER'].notna().sum()}")
    else:
        # Fallback: ENA metadata (gide2019 + riaz2017 only)
        ena_meta = fetch_ena_metadata(cache_path=os.path.join(a.outdir, "ena_run_metadata.parquet"))
        ena_meta = ena_meta[ena_meta["run_accession"].isin(M["run_accession"])]
        xwalk = build_crosswalk(ena_meta, clin)
        ev = attach_labels(M, xwalk, clin)
        print(f"crosswalk (ENA): mapped {len(xwalk)}/{len(M)} runs to study sampleIds; "
              f"RESPONDER labels attached to {ev['RESPONDER'].notna().sum()}")
    feats = [c for c in M.columns if c not in KEY]
    try:
        from analysis.eval.batch_robustness import run_batch_robustness
        run_batch_robustness(ev, [f for f in feats if ev[f].dtype.kind in "fi"],
                             batch_col="cohort", outcome_col="RESPONDER", outdir=a.outdir)
        print("batch-robustness report written")
    except Exception as e:
        print("batch-robustness skipped:", e)
    return 0

if __name__ == "__main__":
    sys.exit(main())
