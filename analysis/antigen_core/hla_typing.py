"""
analysis/antigen_core/hla_typing.py

SHARED ANTIGEN CORE — RNA-based HLA class I genotyping.

Produces, per sample, the 6 HLA-I alleles (A/B/C x2) and the boolean
``HLA_I_heterozygous`` — heterozygous at ALL THREE class-I loci — which is the
Chowell 2017 favorable checkpoint-response feature (doi:10.1126/science.aao4572).

Tool decision (arcasHLA vs OptiType)
------------------------------------
Primary tool: **arcasHLA** (RabadanLab/arcasHLA). It genotypes HLA directly
from an RNA-seq STAR-aligned BAM by extracting chromosome-6 + unmapped reads,
re-quantifying against the IPD-IMGT/HLA reference with kallisto, and calling a
genotype. It is the standard RNA-seq HLA typer and integrates cleanly with the
pipeline session's STAR BAMs (no re-alignment of FASTQ needed).

Portability note (arm64 mac dev box): arcasHLA depends on kallisto and a
git-lfs IMGT/HLA reference. When it cannot be built/run in-sandbox (no osx-arm64
kallisto build, or the reference cannot be fetched offline), the module still
provides the full wrapper + parser, whose parsing and heterozygosity logic are
validated on synthetic genotype fixtures in ``test_hla_typing.py`` (run it with
``PYTHONPATH=. python test_hla_typing.py``). **OptiType** (FRED-2/OptiType) is
the documented fallback: it types HLA-I
from reads via ILP on the IMGT reference and is the nf-core/hlatyping default —
also BAM/FASTQ-driven. Whichever ran for the pilot is recorded in the emitted
manifest (``tool``/``tool_version`` columns) so provenance is explicit; we never
fabricate allele calls for real samples.

Both paths converge on the same tidy output schema below, so downstream code is
tool-agnostic.

Output schema (tidy, one row per sample)
----------------------------------------
    run_accession, cohort,
    HLA_A_1, HLA_A_2, HLA_B_1, HLA_B_2, HLA_C_1, HLA_C_2,
    HLA_I_heterozygous,           # bool: het at A AND B AND C
    n_het_loci,                   # 0-3, for QC
    tool, tool_version            # provenance
"""
from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path
from typing import Iterable

import pandas as pd

LOCI = ("A", "B", "C")
ALLELE_COLS = ["HLA_A_1", "HLA_A_2", "HLA_B_1", "HLA_B_2", "HLA_C_1", "HLA_C_2"]
OUTPUT_COLS = (["run_accession", "cohort"] + ALLELE_COLS
               + ["HLA_I_heterozygous", "n_het_loci", "tool", "tool_version"])


# ---------------------------------------------------------------------------
# Allele normalization + heterozygosity logic (the part that must be correct
# and is validated on synthetic data below)
# ---------------------------------------------------------------------------
def normalize_allele(allele: str, resolution: int = 2) -> str:
    """Normalize an HLA allele call to a comparable field-resolution string.

    arcasHLA emits e.g. 'A*02:01:01' / 'A*02:01'; OptiType emits 'A*02:01'.
    We compare at 2-field ('A*02:01') resolution by default — the resolution
    at which zygosity is biologically defined for antigen presentation, and
    the common denominator across tools/reference versions.

    Handles: 'HLA-A*02:01', 'A*02:01:01:02', 'A_02_01' -> 'A*02:01'.
    Returns '' for empty/NA input.
    """
    if allele is None:
        return ""
    a = str(allele).strip()
    if a == "" or a.lower() in {"nan", "none", "na"}:
        return ""
    a = a.replace("HLA-", "").replace("_", ":").replace("*", "*")
    # split locus from fields: locus letter(s) then '*' then colon-fields
    if "*" in a:
        locus, fields = a.split("*", 1)
    else:
        # e.g. 'A:02:01' -> locus 'A', rest '02:01'
        parts = a.split(":")
        locus, fields = parts[0], ":".join(parts[1:])
    field_parts = [f for f in fields.split(":") if f != ""]
    field_parts = field_parts[:resolution]
    if not field_parts:
        return locus.strip()
    return f"{locus.strip()}*{':'.join(field_parts)}"


def is_heterozygous_locus(allele1: str, allele2: str) -> bool:
    """Two alleles at a locus are heterozygous if they differ at the compared
    resolution and both are present."""
    a1 = normalize_allele(allele1)
    a2 = normalize_allele(allele2)
    if not a1 or not a2:
        return False
    return a1 != a2


def summarize_genotype(genotype: dict[str, list[str]],
                       run_accession: str,
                       cohort: str,
                       tool: str = "",
                       tool_version: str = "") -> dict:
    """Turn a per-locus genotype dict into one tidy output row.

    Parameters
    ----------
    genotype : {"A": [a1, a2], "B": [b1, b2], "C": [c1, c2]}. Missing/homozygous
        loci may provide one allele or two identical alleles; both are handled.

    Returns a dict with the OUTPUT_COLS keys. ``HLA_I_heterozygous`` is True iff
    the sample is heterozygous at A AND B AND C (Chowell 2017 definition).
    """
    row: dict = {"run_accession": run_accession, "cohort": cohort}
    het_flags = []
    for locus in LOCI:
        calls = [normalize_allele(x) for x in genotype.get(locus, [])]
        calls = [c for c in calls if c]
        # pad/truncate to exactly 2 slots; homozygous -> duplicate
        if len(calls) == 0:
            a1, a2 = "", ""
        elif len(calls) == 1:
            a1, a2 = calls[0], calls[0]
        else:
            a1, a2 = calls[0], calls[1]
        row[f"HLA_{locus}_1"] = a1
        row[f"HLA_{locus}_2"] = a2
        het_flags.append(is_heterozygous_locus(a1, a2))
    row["n_het_loci"] = int(sum(het_flags))
    # heterozygous at ALL THREE loci (Chowell favorable feature)
    row["HLA_I_heterozygous"] = bool(all(het_flags))
    row["tool"] = tool
    row["tool_version"] = tool_version
    return row


# ---------------------------------------------------------------------------
# arcasHLA output parser
# ---------------------------------------------------------------------------
def parse_arcashla_genotype(json_path: str | Path) -> dict[str, list[str]]:
    """Parse an arcasHLA ``*.genotype.json`` file into a per-locus dict.

    arcasHLA writes {"A": ["A*02:01:01", "A*01:01:01"], "B": [...], ...} —
    possibly including non-class-I loci (DRB1, DQB1...). We keep only A/B/C.
    A homozygous locus is reported by arcasHLA as a single-element list.
    """
    with open(json_path) as fh:
        data = json.load(fh)
    return {locus: list(data.get(locus, [])) for locus in LOCI}


def parse_optitype_result(tsv_path: str | Path) -> dict[str, list[str]]:
    """Parse an OptiType ``*_result.tsv`` into a per-locus dict.

    OptiType writes a one-row TSV with columns A1 A2 B1 B2 C1 C2 (+ Reads,
    Objective). Values are 2-field alleles like 'A*02:01'.
    """
    df = pd.read_csv(tsv_path, sep="\t")
    r = df.iloc[0]
    return {
        "A": [str(r.get("A1", "")), str(r.get("A2", ""))],
        "B": [str(r.get("B1", "")), str(r.get("B2", ""))],
        "C": [str(r.get("C1", "")), str(r.get("C2", ""))],
    }


# ---------------------------------------------------------------------------
# arcasHLA CLI wrapper (runs when the tool + reference are available)
# ---------------------------------------------------------------------------
def _tool_version(exe: str) -> str:
    try:
        out = subprocess.run([exe, "--version"], capture_output=True,
                             text=True, timeout=30)
        return (out.stdout or out.stderr).strip().splitlines()[0] if (out.stdout or out.stderr) else ""
    except Exception:
        return ""


def run_arcashla(bam: str | Path,
                 outdir: str | Path,
                 threads: int = 4,
                 single_end: bool = False,
                 arcashla_exe: str = "arcasHLA") -> Path:
    """Run arcasHLA extract + genotype on a STAR BAM. Returns the genotype JSON.

    Raises FileNotFoundError if arcasHLA is not on PATH — callers should catch
    this and fall back to OptiType or synthetic validation, per the module
    docstring. This function does not attempt to download the IMGT reference;
    ``arcasHLA reference`` must have been run once beforehand.
    """
    if shutil.which(arcashla_exe) is None:
        raise FileNotFoundError(
            f"{arcashla_exe} not found on PATH. Install RabadanLab/arcasHLA "
            "(needs kallisto + git-lfs IMGT/HLA reference) or use OptiType."
        )
    bam = Path(bam)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    sample = bam.stem

    # 1) extract chr6 + unmapped HLA reads to FASTQ
    subprocess.run(
        [arcashla_exe, "extract", str(bam), "-o", str(outdir),
         "-t", str(threads), "-v"] + ([] if not single_end else ["--single"]),
        check=True,
    )
    # 2) genotype from the extracted FASTQ(s)
    if single_end:
        fq = [str(outdir / f"{sample}.extracted.fq.gz")]
    else:
        fq = [str(outdir / f"{sample}.extracted.1.fq.gz"),
              str(outdir / f"{sample}.extracted.2.fq.gz")]
    subprocess.run(
        [arcashla_exe, "genotype", *fq, "-g", "A,B,C",
         "-o", str(outdir), "-t", str(threads), "-v"],
        check=True,
    )
    gj = outdir / f"{sample}.genotype.json"
    if not gj.exists():
        raise RuntimeError(f"arcasHLA did not produce {gj}")
    return gj


def type_sample(bam: str | Path,
                run_accession: str,
                cohort: str,
                outdir: str | Path,
                threads: int = 4,
                single_end: bool = False) -> dict:
    """End-to-end: run arcasHLA on one BAM -> tidy output row.

    Convenience wrapper for the Nextflow module / pilot driver. On any failure
    (tool missing, reference missing) the exception propagates so the caller
    decides the fallback — we never silently emit fabricated calls.
    """
    gj = run_arcashla(bam, outdir, threads=threads, single_end=single_end)
    genotype = parse_arcashla_genotype(gj)
    return summarize_genotype(
        genotype, run_accession, cohort,
        tool="arcasHLA", tool_version=_tool_version("arcasHLA"),
    )


def build_hla_table(rows: Iterable[dict]) -> pd.DataFrame:
    """Assemble per-sample rows into the contract-format tidy table."""
    df = pd.DataFrame(list(rows))
    for c in OUTPUT_COLS:
        if c not in df.columns:
            df[c] = pd.NA
    return df[OUTPUT_COLS].reset_index(drop=True)
