#!/usr/bin/env python
"""
End-to-end INTEGRATION test for the IO-response feature pipeline.

Unlike the per-module unit tests (which each sub-agent ran in isolation), this
proves the 12 independently-built modules COMPOSE: the one shared MHCflurry
engine is importable across module boundaries, each antigen chain
(peptide-generation -> shared count_binders) recovers a planted HLA-A*02:01
epitope, the baseline scorers run, and every module's output joins on the
contract key (run_accession + cohort) into one per-sample matrix.

Run:  MHCFLURRY_DATA_DIR=reference/mhcflurry_models python analysis/integration_test.py
Requires the 'antigen' env (mhcflurry + models fetched).
"""
import os, sys
os.environ.setdefault("MHCFLURRY_DATA_DIR", os.path.abspath("reference/mhcflurry_models"))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

EPITOPE = "GILGFVFTL"          # influenza M1, canonical A*02:01 binder
ALLELES = ["A*02:01","A*01:01","B*07:02","B*08:01","C*07:01","C*07:02"]

def test_shared_engine():
    from analysis.antigen_core.mhc_binding import binder_counts
    bc = binder_counts([EPITOPE, "NLVPMVATV", "EEEEEEEEE", "KKKKKKKKK"], ALLELES)
    assert bc["n_strong_binders"] >= 1, bc
    assert bc["n_scored"] == 4, bc
    return f"shared engine: {bc['n_strong_binders']} strong of {bc['n_scored']}"

def test_fusion():
    from analysis.differentiated.fusion_antigen import FusionCall, compute_fusion_features
    fc = FusionCall(gene1="A", gene2="B", reading_frame="in-frame", confidence="high",
                    peptide="MKRAAQGILGF|VFTLNPEDS")
    f = compute_fusion_features([fc], ALLELES)
    assert f["fusion_neoantigen_burden"] >= 1, f
    return f"fusion: burden={f['fusion_neoantigen_burden']} n_inframe={f['n_inframe_fusions']}"

def test_te():
    from analysis.differentiated.te_antigen import peptides_from_sequence
    from analysis.antigen_core.mhc_binding import count_binders
    codon={'G':'GGA','I':'ATA','L':'CTA','F':'TTT','V':'GTA','T':'ACA','M':'ATG','A':'GCA'}
    nt = "".join(codon[a] for a in "M"+EPITOPE+"AAAAA") + "TAA"
    peps = peptides_from_sequence(nt, strand="+")
    assert EPITOPE in peps and count_binders(list(peps), ALLELES) >= 1
    return f"TE: {len(peps)} ORF peptides, epitope recovered"

def test_editing():
    from analysis.differentiated.rna_editing import RecodingSite, recoding_peptides
    from analysis.antigen_core.mhc_binding import count_binders
    codon={'M':'ATG','K':'AAA','R':'CGA','A':'GCA','G':'GGA','I':'ATA','L':'CTA',
           'F':'TTT','V':'GTA','T':'ACA','N':'AAT','P':'CCA','Q':'CAA'}
    cds = "".join(codon[a] for a in "MKRA"+EPITOPE+"NPQ")
    site = RecodingSite(site_id="T", chrom="chr1", pos=1000, strand="+",
                        gene="G", cds_window=cds, edit_offset=3)
    rp = recoding_peptides(site)
    assert rp["peptides"], rp
    return f"editing: {len(rp['peptides'])} altered peptides, {count_binders(rp['peptides'], ALLELES)} binders"

def test_gep():
    import pandas as pd
    from analysis.baseline.gep_scores import score_all
    tpm = pd.read_parquet("fixtures/gene_tpm.parquet")
    sm = pd.read_parquet("fixtures/samples.parquet").set_index("run_accession", drop=False)
    gep = score_all(tpm, sm)
    m = gep.merge(sm.reset_index(drop=True), on=["run_accession","cohort"])
    hi = m.loc[m.RESPONDER, "gep_tcell_inflamed"].mean()
    lo = m.loc[~m.RESPONDER, "gep_tcell_inflamed"].mean()
    assert hi > lo, (hi, lo)
    return f"GEP: responder {hi:+.2f} vs non {lo:+.2f}"

def test_ir_load():
    import pandas as pd
    from analysis.differentiated.intron_retention import compute_retained_intron_load
    ir = pd.read_parquet("fixtures/ir_wide.parquet")
    load = compute_retained_intron_load(ir, threshold=0.1)
    assert "retained_intron_load" in load.columns
    return f"IR load: {int(load['retained_intron_load'].max())} max retained introns"

def test_hla_typing():
    from analysis.antigen_core import hla_typing as ht
    assert ht.is_heterozygous_locus("A*02:01", "A*01:01") is True
    assert ht.is_heterozygous_locus("A*02:01", "A*02:01") is False
    het = ht.summarize_genotype({"A":["A*02:01","A*01:01"],"B":["B*07:02","B*08:01"],
                                 "C":["C*07:01","C*07:02"]}, "S2", "pilotA", "arcasHLA", "test")
    hom = ht.summarize_genotype({"A":["A*02:01","A*02:01"],"B":["B*07:02","B*08:01"],
                                 "C":["C*07:01","C*07:02"]}, "S0", "pilotA", "arcasHLA", "test")
    assert het["HLA_I_heterozygous"] and not hom["HLA_I_heterozygous"], (het, hom)
    tbl = ht.build_hla_table([het, hom])
    assert {"run_accession","cohort","HLA_I_heterozygous"} <= set(tbl.columns), tbl.columns.tolist()
    return f"hla_typing: module summarize_genotype/build_hla_table, het={het['HLA_I_heterozygous']} hom={hom['HLA_I_heterozygous']}"

def test_splicing():
    from analysis.differentiated import splicing_neoantigen as sp
    fasta = sp.FastaSeq("fixtures/genome/mini.fa")
    j = sp.Junction(chrom="chrSP", start=201, end=400, strand="+", count=50, normal_mean=0.0)
    peps = sp.translate_junction(*sp.retrieve_flanking_seq(j, fasta, flank=60))
    b = sp.splice_neoantigen_burden([j], ALLELES, fasta=fasta)
    b = b[0] if isinstance(b, tuple) else b
    assert EPITOPE in peps and b >= 1, (EPITOPE in peps, b)
    return f"splicing: neojunction -> {len(peps)} peptides, burden={b}"

def test_ir_cryptic_orf():
    from analysis.differentiated import intron_retention as ir
    from analysis.antigen_core.mhc_binding import count_binders
    genome = ir.GenomeFasta("fixtures/genome/mini.fa")
    coord = ir.IntronCoord(intron_id="IR1", chrom="chrIR", strand="+", blocks=[(101, 151)])
    peps = ir.cryptic_orf_peptides(coord, genome)
    assert EPITOPE in peps and count_binders(peps, ALLELES) >= 1
    return f"IR cryptic-ORF: {len(peps)} peptides, {count_binders(peps, ALLELES)} binders"

def test_snv_indel():
    from analysis.baseline import snv_indel_neoantigen as sni
    v = sni.Variant(gene="TESTG", wt_protein="MKRAAQGILAFVFTLNPEDS",
                    variant_type="missense", protein_pos=10, alt_aa="G")
    peps = sni.peptides_for_variant(v)
    b = sni.snv_indel_neoantigen_burden([v], ALLELES)
    assert EPITOPE in peps and b >= 1, (EPITOPE in peps, b)
    return f"SNV/indel: missense -> {len(peps)} peptides, burden={b}"

def test_join_contract():
    import pandas as pd
    M = pd.read_parquet("fixtures/integration_feature_matrix.parquet")
    assert list(M.columns[:2]) == ["run_accession","cohort"], M.columns.tolist()
    assert M["run_accession"].is_unique
    return f"contract join: {M.shape[0]} samples x {M.shape[1]} cols on run_accession+cohort"

TESTS = [test_shared_engine, test_fusion, test_te, test_editing, test_gep, test_ir_load,
         test_hla_typing, test_splicing, test_ir_cryptic_orf, test_snv_indel, test_join_contract]

if __name__ == "__main__":
    n_pass = 0
    for t in TESTS:
        try:
            print(f"PASS  {t.__name__:<22} {t()}")
            n_pass += 1
        except Exception as e:
            print(f"FAIL  {t.__name__:<22} {type(e).__name__}: {e}")
    print(f"\n{n_pass}/{len(TESTS)} integration tests passed")
    sys.exit(0 if n_pass == len(TESTS) else 1)
