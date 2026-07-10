# Encoder-phase protocol — the learned-representation half (scoped future work)

**Status: NOT a deadline deliverable.** Scoped here so the hybrid half is a designed experiment, not an
afterthought. The interpretable non-reference half was tested this session (`HACKATHON_BRIEF.md`, negative);
this document specifies what the *learned-representation* half must do to be a real test rather than a repeat
of the failures the forensic audit already caught.

## Why it is deliberately deferred
The forensic audit (`ENCODER_EVALUATION_FORENSICS.md`) established two facts that make a rushed encoder run
worse than none:
1. **EVA collapsed to expression.** The per-patient EVA feature was a fixed linear map of gene expression:
   `EVA_SAMP = Wn @ E`, R² = 1.000000 reconstructing 1024 EVA dims from 39 expression PCs, PC1 r=0.957.
   A "learned representation" that is an invertible linear function of expression tests nothing new.
2. **Orthrus never ran at patient scale.** The only Orthrus execution was a 5-gene CPU benchmark
   (B2M/GAPDH/ACTB/CD8A/TP53). The cited "~0.49 Orthrus LOCO" was a misattributed PCA/scVI baseline. Zero
   encoders were ever run on patient-specific or aberrant sequence.

The single lesson: **a learned representation is only worth testing if it ingests information expression
does not contain** — i.e. actual (patient-specific, non-reference) *sequence*, not a re-encoding of the
reference expression vector.

## The one hypothesis this phase tests
> A sequence model that reads a tumor's *own* RNA (including non-reference content — retained introns,
> edited sites, TE-derived transcripts, novel junctions) produces a latent representation that predicts ICB
> response **beyond** both the immune floor AND the interpretable non-reference block tested this session.

If the representation cannot beat the interpretable block, the "learned" part adds nothing; if it cannot beat
the floor, the whole hypothesis is unsupported at the achievable n.

## Design (specified, not run)
1. **Input must be sequence, not expression.** Per patient, extract the actual transcript sequences that
   carry non-reference content — the reads/contigs over retained introns, edited loci, TE-transcribed
   regions, and chimeric junctions already localized by this session's callers. The encoder ingests
   *nucleotide sequence*, never a gene×sample expression matrix.
2. **Encoder candidates (long-context RNA/DNA sequence models):** an SSM-class long-context model
   (e.g. a Mamba/HyenaDNA-style genomic LM) or a specialized RNA foundation model. Selection criterion:
   context length long enough to span a retained-intron+exon unit, and a tokenization that does not silently
   normalize away edits/novel junctions. Verify on a held sequence that the embedding CHANGES when a known
   edited/retained site is toggled — the linearity trap (fact 1) must be actively falsified before any
   downstream fit.
3. **Anti-collapse gate (mandatory pre-registration item):** before any response modeling, fit
   `R²(encoder_embedding ~ expression_PCs)`. **If R² > 0.95 the representation has collapsed to expression
   and the run is aborted** — exactly the EVA failure. Report this number first, always.
4. **Compute:** patient-scale sequence encoding needs GPU. Path = Modal (`byoc:modal`), which needs a
   first-time env image build (documented; not yet built). Budget for a GPU image + a multi-hour encode.
5. **Evaluation:** the SAME pre-registered two-block harness as this session (`analysis/two_block_eval.py`),
   extended to a THIRD block (the encoder embedding), with the identical LOCO / permutation / leakage /
   composition discipline. Success = encoder block beats BOTH floor and interpretable-nonref null, LOCO.

## Honest expectation
Given that the interpretable non-reference features carry no signal at n≈21 (this session) and antigen
*quantity* carried none at n=106, the prior for the encoder half is not favorable at currently-achievable n.
Its value is as a **correctly-designed** test — one that, unlike EVA, cannot pass by re-encoding expression —
so that a negative is informative and a positive is real. It should be run when (a) a GPU image exists and
(b) n is large enough that the positive control (immune floor) is comfortably significant under LOCO.

## Concrete next steps (when resumed)
1. Build the Modal GPU image (`compute_provider` kernel, first-time build).
2. Extract per-patient non-reference transcript sequences from the caller outputs (reuse the localized
   sites in `results/nonref_run/out/<acc>/`).
3. Run the anti-collapse gate on a candidate encoder BEFORE any modeling.
4. Only if the gate passes: encode the cohort, add as block 3, run the pre-registered harness.
