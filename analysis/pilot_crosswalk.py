#!/usr/bin/env python
"""
pilot_crosswalk.py — map pipeline RNA run accessions to the analysis frame's
study sample IDs, and attach response + floor labels.

The pipeline emits per-sample matrices keyed on ENA/SRA ``run_accession``
(e.g. ``ERR2208902``, ``SRR5088861``). The analysis frame keys on study
``sampleId`` (e.g. ``PD02_Pre``, ``Pt72_pre``). No shared key exists, so this
module builds the bridge from ENA run metadata (``sample_alias`` / study patient
id) via deterministic, per-study rules. Every mapped run is validated to land on
a real frame row before labels are attached; unmatched runs are dropped, never
guessed.

Validated on the first pilot (gide2019 n=30 + riaz2017 n=10): all 40 runs mapped,
responder balance consistent with the frame (gide 0.53 vs frame-wide 0.55; riaz
deliberately response-balanced at 0.50).

Rules
-----
gide2019 (ENA PRJEB23709): sample_alias ``PD1_<n>_PRE`` -> ``PD<nn>_Pre``
    (anti-PD1 mono); ``ipiPD1_<n>_PRE`` -> ``iP<nn>_Pre`` (combo). Zero-padded to 2.
riaz2017 (SRA PRJNA356761): sample_title ``Pt<n>_<Pre|On>_...`` -> patientId
    ``Pt<n>``, joined to frame on patientId (pilot is Pre-treatment only).
"""
import os, re
import pandas as pd

KEY = ["run_accession", "cohort"]
LABEL_COLS = ["RESPONDER", "TMB_NONSYNONYMOUS", "TIDE_RESPONDER"]


def _gide_alias_to_sampleid(alias):
    m = re.match(r"(ipiPD1|PD1)_(\d+)_PRE", str(alias))
    if not m:
        return None
    arm, n = m.group(1), int(m.group(2))
    prefix = "iP" if arm == "ipiPD1" else "PD"
    return f"{prefix}{n:02d}_Pre"


def _riaz_title_to_patientid(title):
    m = re.match(r"(Pt\d+)_", str(title))
    return m.group(1) if m else None


def build_crosswalk(ena_meta, analysis_frame):
    """Return a DataFrame [run_accession, cohort, sampleId] for runs that map to
    a real frame row. ``ena_meta`` needs columns run_accession, project,
    sample_alias, sample_title. ``analysis_frame`` needs sampleId, cohort,
    patientId.
    """
    af = analysis_frame
    rows = []

    # gide2019 via sample_alias -> sampleId
    gide = ena_meta[ena_meta["project"] == "PRJEB23709"].copy()
    gide["sampleId"] = gide["sample_alias"].map(_gide_alias_to_sampleid)
    gide["cohort"] = "gide2019"
    frame_ids = set(af.loc[af["cohort"] == "gide2019", "sampleId"].astype(str))
    gide = gide[gide["sampleId"].isin(frame_ids)]
    rows.append(gide[["run_accession", "cohort", "sampleId"]])

    # riaz2017 via patientId
    riaz = ena_meta[ena_meta["project"] == "PRJNA356761"].copy()
    riaz["pt"] = riaz["sample_title"].map(_riaz_title_to_patientid)
    af_riaz = af[af["cohort"] == "riaz2017"].copy()
    af_riaz["pt"] = af_riaz["patientId"].astype(str)
    riaz = riaz.merge(af_riaz[["pt", "sampleId"]], on="pt", how="inner")
    riaz["cohort"] = "riaz2017"
    rows.append(riaz[["run_accession", "cohort", "sampleId"]])

    xwalk = pd.concat(rows, ignore_index=True)
    return xwalk


def build_crosswalk_from_catalog(run_catalog, analysis_frame):
    """Crosswalk from a pipeline ``run_catalog.csv`` (run_accession, cohort,
    patient_id, timepoint, ...) to frame sampleId — covers all cohorts including
    hugo2016. gide patient_ids are in ENA-alias form (PD1_n / ipiPD1_n) and are
    normalized with the same rule as build_crosswalk; hugo/riaz map by patientId.
    Only PRE-treatment runs are kept (baseline).
    """
    af = analysis_frame
    cat = run_catalog.copy()
    if "timepoint" in cat.columns:
        cat = cat[cat["timepoint"].astype(str).str.upper().str.startswith("PRE")]
    rows = []

    # gide: patient_id like 'PD1_13' / 'ipiPD1_24' -> sampleId via alias rule
    gide = cat[cat["cohort"] == "gide2019"].copy()
    gide["sampleId"] = (gide["patient_id"].astype(str) + "_PRE").map(_gide_alias_to_sampleid)
    frame_ids = set(af.loc[af["cohort"] == "gide2019", "sampleId"].astype(str))
    gide = gide[gide["sampleId"].isin(frame_ids)]
    rows.append(gide[["run_accession", "cohort", "sampleId"]])

    # hugo + riaz: patient_id == frame sampleId (hugo) or patientId (riaz)
    for coh in ("hugo2016", "riaz2017"):
        sub = cat[cat["cohort"] == coh].copy()
        af_c = af[af["cohort"] == coh].copy()
        af_c["patientId_str"] = af_c["patientId"].astype(str)
        sub = sub.merge(af_c[["patientId_str", "sampleId"]],
                        left_on="patient_id", right_on="patientId_str", how="inner")
        rows.append(sub[["run_accession", "cohort", "sampleId"]])

    return pd.concat(rows, ignore_index=True).drop_duplicates("run_accession")


def attach_labels(feature_matrix, xwalk, analysis_frame, label_cols=LABEL_COLS):
    """Left-join response + floor labels onto a per-run feature matrix via the
    crosswalk. Runs with no crosswalk entry keep NA labels (never imputed).

    Joins on run_accession only: feature frames often carry a null/placeholder
    ``cohort`` (builders that don't know it set it None), so the authoritative
    cohort comes from the crosswalk, not the feature matrix. Any ``cohort``
    column already on the feature matrix is dropped and replaced.
    """
    lab = xwalk.merge(
        analysis_frame[["sampleId", "cohort"] + [c for c in label_cols if c in analysis_frame.columns]],
        on=["sampleId", "cohort"], how="left",
    )
    keep = ["run_accession", "cohort"] + [c for c in label_cols if c in lab.columns]
    fm = feature_matrix.drop(columns=[c for c in ["cohort"] if c in feature_matrix.columns])
    return fm.merge(lab[keep], on="run_accession", how="left")


def fetch_ena_metadata(projects=("PRJEB23709", "PRJNA356761"), cache_path=None, contact_email=None):
    """Fetch run->sample metadata from ENA for the given study accessions.
    Returns a DataFrame; caches to ``cache_path`` if given.
    """
    import urllib.request, urllib.parse, json
    if cache_path and os.path.exists(cache_path):
        return pd.read_parquet(cache_path)
    fields = "run_accession,sample_accession,sample_alias,sample_title,experiment_title,library_name"
    out = []
    for proj in projects:
        q = {"accession": proj, "result": "read_run", "fields": fields, "format": "json"}
        url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + urllib.parse.urlencode(q)
        with urllib.request.urlopen(url, timeout=60) as r:
            out += [dict(row, project=proj) for row in json.load(r)]
    df = pd.DataFrame(out)
    if cache_path:
        df.to_parquet(cache_path)
    return df
