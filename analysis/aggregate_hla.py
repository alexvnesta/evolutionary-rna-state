#!/usr/bin/env python3
"""Aggregate arcasHLA Class-I genotype JSONs into a per-sample allele table.

Presentation-layer step 1 output. For each accession, collapses the arcasHLA
genotype.json (A/B/C, 2 alleles each) into:
  - a 2-field allele set (e.g. A*24:02) suitable for NetMHCpan/MHCflurry
  - HLA-I heterozygosity count (Chowell 2018 covariate: n distinct 2-field alleles / 6)
Validates every allele against the mhcflurry supported-allele list and flags
any that need nearest-supported fallback.
"""
import json, os, glob, argparse

def two_field(a):
    # 'A*24:02:124' -> 'A*24:02' ; 'A*66:57' -> 'A*66:57'
    parts = a.split(":")
    return ":".join(parts[:2]) if len(parts) >= 2 else a

def mhcflurry_name(a):
    # 'A*24:02' -> 'HLA-A*24:02'
    return "HLA-" + a

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hla_dir", default="results/hla")
    ap.add_argument("--out", default="results/predictor/hla_alleles.parquet")
    ap.add_argument("--supported", default=None,
                    help="optional path to a text file of mhcflurry supported alleles (one per line)")
    args = ap.parse_args()

    supported = None
    if args.supported and os.path.exists(args.supported):
        supported = set(l.strip() for l in open(args.supported) if l.strip())

    rows = []
    for jpath in sorted(glob.glob(os.path.join(args.hla_dir, "*", "*.genotype.json"))):
        acc = os.path.basename(jpath).replace(".genotype.json", "")
        g = json.load(open(jpath))
        rec = {"run_accession": acc}
        alleles_2f = []
        for locus in ("A", "B", "C"):
            calls = g.get(locus, [])
            tf = [two_field(x) for x in calls]
            # pad to 2 (homozygous reported once by arcasHLA sometimes)
            if len(tf) == 1:
                tf = tf * 2
            rec[f"{locus}1"], rec[f"{locus}2"] = (tf + ["", ""])[:2]
            alleles_2f += tf
        distinct = sorted(set(a for a in alleles_2f if a))
        rec["hla_alleles_2field"] = ";".join(distinct)
        rec["n_distinct_classI"] = len(distinct)
        rec["hla_het_fraction"] = len(distinct) / 6.0  # Chowell 2018 heterozygosity proxy
        rec["fully_heterozygous"] = int(len(distinct) == 6)
        if supported is not None:
            mf = [mhcflurry_name(a) for a in distinct]
            rec["n_supported"] = sum(1 for m in mf if m in supported)
            rec["unsupported"] = ";".join(m for m in mf if m not in supported)
        rows.append(rec)

    import pandas as pd
    df = pd.DataFrame(rows).sort_values("run_accession").reset_index(drop=True)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_parquet(args.out, index=False)
    df.to_csv(args.out.replace(".parquet", ".csv"), index=False)
    print(f"wrote {len(df)} samples -> {args.out}")
    print(df.to_string(index=False))
    return df

if __name__ == "__main__":
    main()