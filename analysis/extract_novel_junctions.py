"""Extract NOVEL splice junctions per sample (sequence-visible aberrancy, protocol §1).
Known introns from GENCODE v46; a regtools junction is NOVEL if its (chrom, donor, acceptor)
intron boundary is not in the known set. Emits per-sample novel-junction tables + a deduped
unique-junction table (sequence to score once with Evo2)."""
import os, glob, gzip, json, collections
os.chdir("/Users/alex/OrchestratedBiosciences/evolutionary-rna-state")
GTF="reference/GRCh38/gencode.v46.primary_assembly.annotation.gtf"
MIN_READS=3          # min regtools score (read support) for a junction to count
MIN_INTRON=50        # ignore tiny gaps
MAX_INTRON=500000    # ignore ultra-long (likely artifact)

# 1) known introns: per-transcript, adjacent exon pairs -> intron (donor,acceptor)
# regtools BED: chromStart=last base of upstream exon region span; the thickStart/blockSizes encode
# the actual intron. Standard: intron = (exon_i.end, exon_{i+1}.start). Build known set of (chrom,istart,iend).
known=set()
tx_exons=collections.defaultdict(list)
n=0
with open(GTF) as fh:
    for line in fh:
        if line[0]=="#": continue
        f=line.split("\t")
        if f[2]!="exon": continue
        chrom=f[0]; start=int(f[3]); end=int(f[4])
        # transcript_id
        attr=f[8]; ti=attr.find('transcript_id "')
        tid=attr[ti+15:attr.find('"',ti+15)] if ti>=0 else None
        if tid: tx_exons[(chrom,tid)].append((start,end))
        n+=1
for (chrom,tid),ex in tx_exons.items():
    ex.sort()
    for i in range(len(ex)-1):
        istart=ex[i][1]+1     # first base of intron (1-based)
        iend=ex[i+1][0]-1     # last base of intron
        if iend>istart: known.add((chrom,istart,iend))
print(f"parsed {n} exons, {len(tx_exons)} transcripts, {len(known)} known introns", flush=True)
json.dump({"n_known_introns":len(known)}, open("results/eval/_known_intron_count.json","w"))

# 2) per-sample novel junctions
jdirs=sorted(glob.glob("results/nonref_run/out/*/junctions.bed"))
os.makedirs("results/junctions", exist_ok=True)
uniq=collections.defaultdict(lambda:{"n_samples":0,"total_reads":0})  # (chrom,istart,iend,strand)->agg
per_sample={}
for jd in jdirs:
    acc=jd.split("/")[-2]
    novel=[]
    with open(jd) as fh:
        for line in fh:
            c=line.rstrip("\n").split("\t")
            if len(c)<12: continue
            chrom=c[0]; chromStart=int(c[1]); score=int(c[4]); strand=c[5]
            # regtools: intron = chromStart + blockSizes[0] .. chromEnd - blockSizes[1]
            bsizes=[int(x) for x in c[10].rstrip(",").split(",")]
            if len(bsizes)<2: continue
            istart=chromStart+bsizes[0]+1        # 1-based first intron base
            iend=int(c[2])-bsizes[1]             # last intron base
            ilen=iend-istart+1
            if score<MIN_READS or ilen<MIN_INTRON or ilen>MAX_INTRON: continue
            if (chrom,istart,iend) in known: continue
            novel.append((chrom,istart,iend,strand,score))
            k=(chrom,istart,iend,strand)
            uniq[k]["n_samples"]+=1; uniq[k]["total_reads"]+=score
    per_sample[acc]={"n_novel":len(novel),"total_novel_reads":sum(x[4] for x in novel)}
    # write per-sample novel junctions
    with open(f"results/junctions/{acc}.novel.tsv","w") as o:
        o.write("chrom\tistart\tiend\tstrand\treads\n")
        for x in novel: o.write(f"{x[0]}\t{x[1]}\t{x[2]}\t{x[3]}\t{x[4]}\n")

# 3) unique novel junction table (to score once)
with open("results/junctions/unique_novel.tsv","w") as o:
    o.write("chrom\tistart\tiend\tstrand\tn_samples\ttotal_reads\n")
    for (chrom,istart,iend,strand),v in uniq.items():
        o.write(f"{chrom}\t{istart}\t{iend}\t{strand}\t{v['n_samples']}\t{v['total_reads']}\n")

json.dump(per_sample, open("results/junctions/per_sample_summary.json","w"), indent=2)
print("unique novel junctions:", len(uniq), flush=True)
import numpy as np
nn=[v["n_novel"] for v in per_sample.values()]
print(f"per-sample novel count: median={int(np.median(nn))} min={min(nn)} max={max(nn)}", flush=True)
