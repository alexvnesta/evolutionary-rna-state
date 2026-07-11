"""For each novel junction, extract two sequences for Evo2 delta-likelihood:
  - spliced:    [donor-side flank] + [acceptor-side flank]  (the aberrant transcript across the junction)
  - contiguous: genomic sequence spanning donor..donor+2*FLANK (reference, no splice)
FLANK bp each side. Reverse-complement for '-' strand. Writes FASTA + index."""
import json, os
os.chdir("/Users/alex/OrchestratedBiosciences/evolutionary-rna-state")
import pysam
FLANK=512  # Evo2 handles long context; 512+512 spliced window
fa=pysam.FastaFile("reference/GRCh38/GRCh38.primary_assembly.genome.fa")
uniq=json.load(open("results/junctions/gide_top200_unique.json"))
def rc(s):
    return s.translate(str.maketrans("ACGTNacgtn","TGCANtgcan"))[::-1]
def fetch(chrom,a,b):
    try: return fa.fetch(chrom,max(0,a),b).upper()
    except Exception: return ""
recs=[]
for i,(chrom,istart,iend,strand) in enumerate(uniq):
    istart=int(istart); iend=int(iend)
    # donor side = just before intron start; acceptor side = just after intron end
    donor = fetch(chrom, istart-1-FLANK, istart-1)     # exonic bases up to donor
    acceptor = fetch(chrom, iend, iend+FLANK)           # exonic bases after acceptor
    spliced = donor + acceptor
    contiguous = fetch(chrom, istart-1-FLANK, istart-1+FLANK)  # reference contiguous, same length
    if strand=="-":
        spliced=rc(spliced); contiguous=rc(contiguous)
    if len(spliced)<200 or len(contiguous)<200: continue
    jid=f"J{i}"
    recs.append({"jid":jid,"chrom":chrom,"istart":istart,"iend":iend,"strand":strand,
                 "spliced":spliced,"contiguous":contiguous})
json.dump(recs, open("results/junctions/junction_seqs.json","w"))
# also a compact FASTA for scoring (spliced + contiguous)
with open("results/junctions/junction_seqs.fasta","w") as o:
    for r in recs:
        o.write(f">{r['jid']}_spliced\n{r['spliced']}\n>{r['jid']}_contig\n{r['contiguous']}\n")
print("extracted", len(recs), "junction sequence pairs")
import numpy as np
print("spliced len:", int(np.median([len(r['spliced']) for r in recs])),
      "contig len:", int(np.median([len(r['contiguous']) for r in recs])))
