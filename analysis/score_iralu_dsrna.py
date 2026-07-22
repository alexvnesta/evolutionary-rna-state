#!/usr/bin/env python
"""IR-Alu dsRNA feature: score hybridization ΔG of proximal opposite-strand Alu pairs
inside expressed gene bodies. The genomic IR-Alu map is FIXED across samples; only the
per-sample expression weighting varies -> by construction not a bulk immune-abundance readout.

Grounded in the IRAlu / viral-mimicry model: opposite-strand Alu copies in the same
transcript are reverse-complementary and form dsRNA (RIG-I/MDA5/TLR3 ligands).
Output: per-gene ΔG summaries -> aggregated per-sample downstream.
"""
import sys, re, time
import numpy as np, pandas as pd
from multiprocessing import Pool
import RNA
from pyfaidx import Fasta

BASE = "/Users/alex/OrchestratedBiosciences/evolutionary-rna-state"
GAP = 1000
FA_PATH = f"{BASE}/reference/GRCh38/GRCh38.primary_assembly.genome.fa"

def build_pairs():
    alu = pd.read_csv(f"{BASE}/reference/GRCh38/repeats/alu.hg38.bed6", sep="\t",
                      names=["chrom","start","end","name","score","strand"])
    alu = alu[alu.chrom.str.match(r'^chr([0-9]+|X|Y)$')].copy()
    rows = []
    for chrom, g in alu.sort_values(["chrom","start"]).groupby("chrom", sort=False):
        s = g.reset_index(drop=True)
        st, en, strand = s.start.values, s.end.values, s.strand.values
        n = len(s)
        for i in range(n):
            j = i+1
            while j < n and st[j] - en[i] <= GAP:
                if st[j] - en[i] >= 0 and strand[j] != strand[i]:
                    rows.append((chrom, st[i], en[i], st[j], en[j], st[j]-en[i]))
                j += 1
    return pd.DataFrame(rows, columns=["chrom","s1","e1","s2","e2","gap"])

def parse_genes():
    genes = []
    with open(f"{BASE}/reference/GRCh38/gencode.v46.primary_assembly.annotation.gtf") as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.split("\t")
            if f[2] != "gene": continue
            m = re.search(r'gene_id "([^".]+)', f[8])
            if m: genes.append((m.group(1), f[0], int(f[3]), int(f[4])))
    return pd.DataFrame(genes, columns=["ensg","chrom","start","end"])

def assign(pdf, gdf_e):
    assigned = np.full(len(pdf), None, dtype=object)
    posmap = {pi: k for k, pi in enumerate(pdf.index)}
    for chrom, gg in gdf_e.groupby("chrom", sort=False):
        idx = pdf.index[pdf.chrom == chrom]
        if len(idx) == 0: continue
        gs, ge, ens = gg.start.values, gg.end.values, gg.ensg.values
        order = np.argsort(gs); gs, ge, ens = gs[order], ge[order], ens[order]
        for pi in idx:
            a, b = pdf.at[pi,"s1"], pdf.at[pi,"e2"]
            k = np.searchsorted(gs, a, side="right")-1
            j = k
            while j >= 0 and gs[j] >= a-3_000_000:
                if gs[j] <= a and ge[j] >= b:
                    assigned[posmap[pi]] = ens[j]; break
                j -= 1
    pdf = pdf.copy(); pdf["ensg"] = assigned
    return pdf[pdf.ensg.notna()].reset_index(drop=True)

_FA = None
def _init():
    global _FA
    _FA = Fasta(FA_PATH, sequence_always_upper=True)

def _score(args):
    chrom, s1, e1, s2, e2 = args
    a = str(_FA[chrom][int(s1):int(e1)]).replace("T","U")
    b = str(_FA[chrom][int(s2):int(e2)]).replace("T","U")
    if len(a) < 20 or len(b) < 20: return np.nan
    try:
        return RNA.duplexfold(a, b).energy
    except Exception:
        return np.nan

def main():
    t0 = time.time()
    pdf = build_pairs()
    gdf = parse_genes()
    tpm = pd.read_parquet(f"{BASE}/results/features/quant_gene_tpm.parquet")
    tpm_ensg = set(c for c in tpm.columns if c.startswith("ENSG"))
    gdf_e = gdf[gdf.ensg.isin(tpm_ensg)].copy()
    ig = assign(pdf, gdf_e)
    print(f"[{time.time()-t0:.0f}s] scoring {len(ig)} in-gene inverted-Alu pairs", flush=True)
    args = list(zip(ig.chrom, ig.s1, ig.e1, ig.s2, ig.e2))
    with Pool(16, initializer=_init) as pool:
        ig["dG"] = pool.map(_score, args, chunksize=256)
    ig = ig.dropna(subset=["dG"])
    ig.to_parquet(f"{BASE}/results/features/iralu_pairs_scored.parquet")
    # per-gene summary
    gsum = ig.groupby("ensg").agg(
        n_ir_pairs=("dG","size"),
        dG_sum=("dG","sum"),
        dG_mean=("dG","mean"),
        dG_min=("dG","min"),
    ).reset_index()
    gsum.to_parquet(f"{BASE}/results/features/iralu_pergene.parquet")
    print(f"[{time.time()-t0:.0f}s] DONE: {len(ig)} pairs, {gsum.ensg.nunique()} genes", flush=True)

if __name__ == "__main__":
    main()
