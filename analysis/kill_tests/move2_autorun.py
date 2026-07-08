#!/usr/bin/env python
"""
move2_autorun.py — watch-and-execute trigger for Move 2 (regulator-activity tests).

WHAT IT DOES (idempotent, safe to re-run):
  1. Finds the latest quant_gene_tpm.parquet in the artifact store.
  2. Checks a state file (.move2_state.json) for the version_id last processed and
     the sample count. If unchanged since last run -> exits 0 with 'no change'.
  3. If new/grown: scores regulator activity (rbp_activity_scorer), builds the
     per-cohort run->iAtlas crosswalk (gide: reuse the pre-validated
     gide2019_id_crosswalk.csv artifact; riaz: derive this run from ENA PRJNA356761
     sample titles, validated 10/10 -> iAtlas), joins to analysis_frame, and runs
     the three tests (run_activity_response_tests): mechanistic / shared-state /
     LOCO-AUROC-over-TMB.
  4. Writes results_move2_<n>samp_* and updates the state file.

GATING PREDICATE ("has the full matrix landed?"):
  The pipeline is filling cohorts incrementally. This runner processes WHATEVER
  is present each time and records progress, so it produces the best-available
  result now (40 samples: gide+riaz) and re-fires as liu/hugo/dfci land. The
  headline LOCO-over-TMB test only becomes meaningful once >=2 cohorts carry TMB
  (gide has none); the runner reports readiness per test rather than blocking.

USAGE
  # one-shot (any session): run it; it no-ops if nothing changed
  python move2_autorun.py
  # force re-run even if unchanged
  python move2_autorun.py --force
  # in-kernel (live watch): from move2_autorun import run_if_changed; run_if_changed(host)

Requires (same dir, all in project store): rbp_activity_scorer.py,
run_activity_response_tests.py, analysis_frame.parquet, gide2019_id_crosswalk.csv.
"""
from __future__ import annotations
import json, os, sys, io, re, time, urllib.request
import numpy as np
import pandas as pd

STATE = ".move2_state.json"
REG_GENES_CACHE = ".regulator_ensembl_map.json"


# ---------------------------------------------------------------- helpers
def _ena_run_table(accession, fields):
    url = (f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={accession}"
           f"&result=read_run&fields={fields}&format=tsv")
    raw = urllib.request.urlopen(url, timeout=90).read().decode()
    return pd.read_csv(io.StringIO(raw), sep="\t")


def regulator_ensembl_map(symbols):
    """HGNC symbol -> ENSG via Ensembl REST (cached to disk)."""
    if os.path.exists(REG_GENES_CACHE):
        cached = json.load(open(REG_GENES_CACHE))
        if all(s in cached for s in symbols):
            return {s: cached[s] for s in symbols}
    url = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"
    req = urllib.request.Request(url, data=json.dumps({"symbols": symbols}).encode(),
                                 headers={"Content-Type": "application/json",
                                          "Accept": "application/json"})
    r = json.load(urllib.request.urlopen(req, timeout=90))
    out = {s: info["id"] for s, info in r.items()
           if isinstance(info, dict) and str(info.get("id", "")).startswith("ENSG")}
    json.dump(out, open(REG_GENES_CACHE, "w"))
    return out


def build_crosswalk(cohorts_present, gide_xw_path="gide2019_id_crosswalk.csv"):
    """Return unified run_accession -> (iatlas_patientId, cohort, timepoint,
    iatlas_burden_available) for whichever cohorts are in the matrix.
    gide: reuse the validated crosswalk artifact. riaz: derive from ENA titles."""
    parts = []
    if "gide2019" in cohorts_present and os.path.exists(gide_xw_path):
        g = pd.read_csv(gide_xw_path)
        gg = g[["run_accession", "iatlas_patientId", "timepoint",
                "iatlas_burden_available"]].copy()
        gg["cohort"] = "gide2019"
        parts.append(gg)
    if "riaz2017" in cohorts_present:
        rena = _ena_run_table("PRJNA356761",
                              "run_accession,sample_title")
        pt = rena["sample_title"].str.extract(r"^(Pt\d+)_(Pre|On)_")
        rena["iatlas_patientId"] = pt[0]
        rena["timepoint"] = pt[1].str.upper()
        rr = rena.dropna(subset=["iatlas_patientId"]).copy()
        rr["cohort"] = "riaz2017"
        rr["iatlas_burden_available"] = rr["timepoint"].eq("PRE")
        parts.append(rr[["run_accession", "iatlas_patientId", "timepoint",
                         "cohort", "iatlas_burden_available"]])
    if "hugo2016" in cohorts_present:
        # hugo expression = ENA/SRA PRJNA312948 (GEO GSE78220); sample_title IS the iAtlas
        # patientId ('Pt2'), with a couple of split samples 'Pt27A'/'Pt27B' -> 'Pt27'.
        # Validated: 27/28 matrix runs map into iAtlas hugo (Pt16 absent from response set).
        hena = _ena_run_table("PRJNA312948", "run_accession,sample_title")
        hena["iatlas_patientId"] = (hena["sample_title"].astype(str).str.strip()
                                    .str.replace(r"(Pt\d+)[AB]$", r"\1", regex=True))
        hh = hena.dropna(subset=["iatlas_patientId"]).drop_duplicates("iatlas_patientId").copy()
        hh["cohort"] = "hugo2016"; hh["timepoint"] = "PRE"
        hh["iatlas_burden_available"] = True
        parts.append(hh[["run_accession", "iatlas_patientId", "timepoint",
                         "cohort", "iatlas_burden_available"]])
    # liu expression:
    #  - liu2019 arrives "via DFCI matrix (processed)" (cohorts.csv): its run_accession
    #    column is expected to already carry iAtlas-compatible sample/patient IDs, so the
    #    safe rule is PASSTHROUGH — the caller validates against the iAtlas frame and the
    #    runner FLAGS any cohort whose IDs don't map (rather than silently dropping them).
    #  - hugo2016 raw-read accession in cohorts.csv (PRJNA356839) resolves to an unrelated
    #    submission on inspection, so no title rule is hardcoded; passthrough + flag applies.
    # These cohorts are handled generically in run_if_changed via _passthrough_cohorts.
    xw = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    return xw


def _passthrough_crosswalk(runs, cohort_name, frame=None):
    """For cohorts whose expression matrix already carries iAtlas-style IDs (liu processed
    matrix; hugo fallback): map run_accession to the iAtlas patientId JOIN KEY.

    The analysis_frame join keys on patientId, but a processed matrix may carry the SAMPLE
    id instead (e.g. liu2019: patientId='Patient4' but sampleId='Liu_Sample4'). So when a
    frame is provided we resolve each run_accession against BOTH the cohort's patientId and
    sampleId columns and emit the patientId. IDs that match neither are still emitted as-is
    (so the caller's diagnostic flags them as unmapped rather than silently dropping)."""
    runs = list(runs)
    id2patient = {}
    if frame is not None:
        f = frame[frame["cohort"] == cohort_name]
        for pid in f["patientId"].astype(str):
            id2patient[pid] = pid
        if "sampleId" in f.columns:
            for sid, pid in zip(f["sampleId"].astype(str), f["patientId"].astype(str)):
                id2patient[sid] = pid
    mapped = [id2patient.get(str(r), str(r)) for r in runs]
    return pd.DataFrame({"run_accession": runs, "iatlas_patientId": mapped,
                         "timepoint": "PRE", "cohort": cohort_name,
                         "iatlas_burden_available": True})


# ---------------------------------------------------------------- core
def run_if_changed(host, force=False, n_perm=5000):
    """Detect matrix, and if changed since last run, execute Move 2. Returns a
    dict summary. `host` is the kernel host singleton (for artifact resolution)."""
    import importlib.util

    # locate matrix
    arts = host.artifacts(filename="quant_gene_tpm.parquet", exact=True)["artifacts"]
    if not arts:
        return {"status": "no_matrix", "msg": "quant_gene_tpm.parquet not in store"}
    art = arts[0]
    vid = art["latest_version_id"]

    prev = json.load(open(STATE)) if os.path.exists(STATE) else {}
    # cheap short-circuit: if the version_id is unchanged, don't even read the parquet
    if not force and prev.get("version_id") == vid:
        return {"status": "no_change", "version_id": vid,
                "n_samples": prev.get("n_samples"), "cohorts": prev.get("cohorts"),
                "msg": "matrix version unchanged since last run"}
    q = pd.read_parquet(host.artifact_path(vid))
    cohorts = sorted(q["cohort"].unique().tolist()) if "cohort" in q else []
    fingerprint = {"version_id": vid, "n_samples": int(len(q)), "cohorts": cohorts}

    if (not force and prev.get("version_id") == vid
            and prev.get("n_samples") == len(q)):
        return {"status": "no_change", **fingerprint,
                "msg": "matrix unchanged since last run"}

    # --- score activity ---
    spec = importlib.util.spec_from_file_location("rbp", "rbp_activity_scorer.py")
    rbp = importlib.util.module_from_spec(spec); spec.loader.exec_module(rbp)
    allsyms = [g for s in rbp.REGULATOR_SETS.values() for g in s]
    gmap = regulator_ensembl_map(allsyms)
    ensg = [c for c in q.columns if c.startswith("ENSG")]
    tpm = q.set_index("run_accession")[ensg]
    S = rbp.score_regulator_activity(tpm, gmap).reset_index()
    S = S.rename(columns={S.columns[0]: "run_accession"})
    act_path = f"activity_scores_{len(q)}samp.csv"
    S.to_csv(act_path, index=False)

    # --- crosswalk (validated rules) + passthrough for processed-matrix cohorts ---
    xw = build_crosswalk(cohorts)
    covered = set(xw["cohort"].unique()) if len(xw) else set()
    frame_for_ids = pd.read_parquet(host.artifact_path(
        host.artifacts(filename="analysis_frame.parquet", exact=True)["artifacts"][0]["latest_version_id"]))
    id_diag = {}
    extra = []
    for coh in cohorts:
        if coh in covered:
            # report mapping rate for validated cohorts too
            runs = set(q[q["cohort"] == coh]["run_accession"])
            mapped = xw[(xw.cohort == coh) & (xw.run_accession.isin(runs))]
            ia = set(frame_for_ids[frame_for_ids.cohort == coh]["patientId"].astype(str))
            id_diag[coh] = {"runs": len(runs),
                            "mapped_to_iatlas": int(mapped.iatlas_patientId.isin(ia).sum()),
                            "rule": "validated"}
            continue
        # processed-matrix / unmapped cohort -> passthrough, then validate
        runs = q[q["cohort"] == coh]["run_accession"]
        pt = _passthrough_crosswalk(runs, coh, frame=frame_for_ids)
        ia = set(frame_for_ids[frame_for_ids.cohort == coh]["patientId"].astype(str))
        n_ok = int(pt.iatlas_patientId.isin(ia).sum())
        id_diag[coh] = {"runs": len(runs), "mapped_to_iatlas": n_ok, "rule": "passthrough"}
        if n_ok > 0:
            extra.append(pt)  # only include if passthrough IDs actually resolve
    if extra:
        xw = pd.concat([xw] + extra, ignore_index=True)
    xw_path = "combined_id_crosswalk.csv"
    xw.to_csv(xw_path, index=False)

    spec = importlib.util.spec_from_file_location("h", "run_activity_response_tests.py")
    H = importlib.util.module_from_spec(spec); spec.loader.exec_module(H)
    frame_path = host.artifact_path(
        host.artifacts(filename="analysis_frame.parquet", exact=True)["artifacts"][0]["latest_version_id"])
    m, setcols, meta = H.load_and_join(act_path, xw_path, frame_path)
    t1 = H.test_mechanistic(m, setcols)
    t2 = H.test_shared_state(m, setcols, n_perm=n_perm)
    t3 = H.test_prediction(m, setcols)

    # --- reframing tests: immune-floor coupling + incremental value (P1-P3) ---
    reframe = None
    try:
        spec = importlib.util.spec_from_file_location("rf", "reframing_tests.py")
        RF = importlib.util.module_from_spec(spec); spec.loader.exec_module(RF)
        reg_ensg = [e for e in gmap.values() if e in tpm.columns]
        resp_df = xw[["run_accession", "iatlas_patientId", "cohort"]].merge(
            frame_for_ids[["patientId", "cohort", "RESPONDER"]].rename(
                columns={"patientId": "iatlas_patientId"}),
            on=["iatlas_patientId", "cohort"], how="left")
        reframe, _ = RF.run(tpm, reg_ensg, resp_df, n_perm=min(n_perm, 2000))
    except Exception as e:
        reframe = {"error": f"{type(e).__name__}: {e}"}

    # --- untested-core: learned-representation vs immune-floor head-to-head (LOCO perm) ---
    headtohead = None
    try:
        spec = importlib.util.spec_from_file_location("hh", "headtohead_repr.py")
        HH = importlib.util.module_from_spec(spec); spec.loader.exec_module(HH)
        rmap = xw[["run_accession", "iatlas_patientId", "cohort"]].merge(
            frame_for_ids[["patientId", "cohort", "RESPONDER"]].rename(
                columns={"patientId": "iatlas_patientId"}),
            on=["iatlas_patientId", "cohort"], how="left").dropna(subset=["RESPONDER"])
        rmap = rmap[rmap.run_accession.isin(tpm.index)].drop_duplicates("run_accession")
        logX_hh = np.log2(tpm.loc[rmap.run_accession] + 1)
        headtohead = HH.run(logX_hh, rmap.RESPONDER.astype(int).values,
                            rmap.cohort.values, n_perm=min(n_perm, 300))
    except Exception as e:
        headtohead = {"error": f"{type(e).__name__}: {e}"}

    out = f"results_move2_{len(q)}samp"
    t1.to_csv(f"{out}_T1_mechanistic.csv", index=False)
    report = {"fingerprint": fingerprint, "join_meta": meta,
              "id_mapping_diagnostics": id_diag,
              "T2_shared_state": t2, "T3_prediction": t3,
              "reframing_P1P3": reframe,
              "headtohead_learned_rep": headtohead,
              "n_cohorts_with_TMB": int(m.groupby("cohort")["TMB_NONSYNONYMOUS"]
                                        .apply(lambda s: s.notna().any()).sum())
              if "cohort" in m else None}
    json.dump(report, open(f"{out}_report.json", "w"), indent=2, default=str)
    json.dump({**fingerprint, "ran_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
               "out_prefix": out}, open(STATE, "w"), indent=2)

    hh_loco = (headtohead.get("loco", {}) if isinstance(headtohead, dict) else {}) or {}
    hh_p = (headtohead.get("loco_perm_p", {}) if isinstance(headtohead, dict) else {}) or {}
    return {"status": "ran", **fingerprint, "out_prefix": out,
            "join_n": meta.get("n_joined"),
            "id_mapping": id_diag,
            "T2_p": t2.get("p_value"),
            "T3_activity_auroc": t3.get("activity_alone", {}).get("auroc"),
            "headtohead_loco": {k: hh_loco.get(k) for k in ("immune_floor", "learned_rep", "floor+rep")},
            "headtohead_loco_perm_p": {k: hh_p.get(k) for k in ("immune_floor", "learned_rep", "floor+rep")}}


if __name__ == "__main__":
    # CLI path needs a host; only works inside the kernel. Guard clearly.
    force = "--force" in sys.argv
    try:
        host  # noqa: F821  (injected in kernel)
    except NameError:
        print("move2_autorun: run inside the analysis kernel "
              "(from move2_autorun import run_if_changed; run_if_changed(host))",
              file=sys.stderr)
        sys.exit(2)
    print(json.dumps(run_if_changed(host, force=force), indent=2, default=str))
