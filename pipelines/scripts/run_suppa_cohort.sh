#!/usr/bin/env bash
# SUPPA2 event-level alternative-splicing analysis for the 40-sample cohort.
#
# WHY SUPPA2 (not DEXSeq/rMATS) for this cohort:
#   - 40 tumor samples, biological replicates, 21 responder vs 19 non-responder.
#   - We already have Salmon transcript quantifications for all 40 (rnaseq spine),
#     which SUPPA2 consumes directly — no re-alignment, no single-threaded htseq
#     counting bottleneck (DEXSeq's dexseq_count is ~25 min/sample serial).
#   - rMATS's paired model (r-pairadise) is conda-unsatisfiable on osx-arm64.
#   - SUPPA2 gives interpretable per-event dPSI for all 7 AS event classes
#     (SE, A5, A3, MX, RI, AF, AL) — the splicing phenotype features we want.
#
# Outputs (results/suppa_cohort/):
#   events/*.ioe            AS event definitions from the GENCODE GTF
#   psi/cohort.psi          PSI per event x 40 samples
#   diff/*.dpsi             dPSI + p-value per event, responder vs non-responder
set -euo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
SENV=/Users/alex/.claude-science/conda/envs/suppa
PY="$SENV/bin/python"; SUPPA="$SENV/bin/suppa.py"
export PYTHONPATH="$SENV/bin"                 # SUPPA ships its modules in bin/
run_suppa(){ "$PY" "$SUPPA" "$@"; }

GTF="$REPO/reference/GRCh38/gencode.v46.primary_assembly.annotation.gtf"
SALMON="$REPO/results/rnaseq_cohort/salmon"
OUT="$REPO/results/suppa_cohort"
MANIFEST="$REPO/data/manifests/selection_manifest.csv"
mkdir -p "$OUT/events" "$OUT/psi" "$OUT/diff"

# 1) Generate AS events from the GTF (once). -f ioe = local events for psiPerEvent.
if [ ! -f "$OUT/events/cohort.events.ioe" ]; then
  echo "[suppa] generateEvents (SE,SS,MX,RI,FL) $(date '+%H:%M:%S')"
  run_suppa generateEvents -i "$GTF" -o "$OUT/events/cohort.events" \
    -f ioe -e SE SS MX RI FL --pool-genes
  # SUPPA writes one .ioe per class (e.g. cohort.events_SE_strict.ioe); concat to one
  # file, keeping a single header line.
  awk 'FNR==1 && NR!=1 {next} {print}' "$OUT"/events/cohort.events_*_strict.ioe \
    > "$OUT/events/cohort.events.ioe"
fi
echo "[suppa] events: $(wc -l < "$OUT/events/cohort.events.ioe") lines"

# 2) Build the combined TPM matrix (rows=transcript, cols=40 samples) from salmon quant.sf.
#    SUPPA needs a TPM table keyed by transcript_id; use its own multipleFieldSelection
#    helper via a small python join (robust to sample ordering).
"$PY" - "$SALMON" "$OUT/psi/cohort.tpm" <<'PYEOF'
import sys, os, glob, pandas as pd
salmon, out = sys.argv[1], sys.argv[2]
qfs = sorted(glob.glob(os.path.join(salmon, "*", "quant.sf")))
mat = None
for q in qfs:
    sample = os.path.basename(os.path.dirname(q))
    df = pd.read_csv(q, sep="\t", usecols=["Name","TPM"]).set_index("Name")["TPM"].rename(sample)
    mat = df.to_frame() if mat is None else mat.join(df, how="outer")
mat = mat.fillna(0.0)
# SUPPA expects the header to list ONLY the sample names (no index/label field),
# while each data row is transcript_id + N values. Write with no leading tab.
with open(out, "w") as fh:
    fh.write("\t".join(mat.columns) + "\n")
    mat.to_csv(fh, sep="\t", header=False)
print(f"[tpm] {mat.shape[0]} transcripts x {mat.shape[1]} samples -> {out}")
PYEOF

# 3) PSI per event across all 40 samples.
echo "[suppa] psiPerEvent $(date '+%H:%M:%S')"
run_suppa psiPerEvent --ioe-file "$OUT/events/cohort.events.ioe" \
  --expression-file "$OUT/psi/cohort.tpm" -o "$OUT/psi/cohort" || true
echo "[suppa] psi: $(wc -l < "$OUT/psi/cohort.psi" 2>/dev/null || echo NA) events"

# 4) Split TPM + PSI into responder / non-responder columns (from manifest resp_NR).
"$PY" - "$MANIFEST" "$OUT" <<'PYEOF'
import sys, pandas as pd, os
manifest, out = sys.argv[1], sys.argv[2]
m = pd.read_csv(manifest)
# map sample_title -> resp_NR (R / N)
lab = dict(zip(m["sample_title"].astype(str), m["resp_NR"].astype(str)))
tpm = pd.read_csv(os.path.join(out,"psi","cohort.tpm"), sep="\t", index_col=0)
psi = pd.read_csv(os.path.join(out,"psi","cohort.psi"), sep="\t", index_col=0)
R = [c for c in psi.columns if lab.get(c,"")=="R"]
N = [c for c in psi.columns if lab.get(c,"")=="N"]
print(f"[split] R={len(R)} N={len(N)} (unlabeled={len(psi.columns)-len(R)-len(N)})")
def write_suppa(df, path):
    # header = column names only (no index label), data = id + values
    with open(path, "w") as fh:
        fh.write("\t".join(df.columns) + "\n")
        df.to_csv(fh, sep="\t", header=False)
for tag, cols in (("R",R),("N",N)):
    write_suppa(tpm[cols], os.path.join(out,"psi",f"{tag}.tpm"))
    write_suppa(psi[cols], os.path.join(out,"psi",f"{tag}.psi"))
PYEOF

# 5) Differential splicing R vs N (empirical, classic SUPPA method).
echo "[suppa] diffSplice $(date '+%H:%M:%S')"
run_suppa diffSplice -m empirical -gc \
  -i "$OUT/events/cohort.events.ioe" \
  -p "$OUT/psi/R.psi" "$OUT/psi/N.psi" \
  -e "$OUT/psi/R.tpm" "$OUT/psi/N.tpm" \
  -o "$OUT/diff/resp_vs_nonresp" || true
echo "[suppa] done $(date '+%H:%M:%S'); dpsi: $(ls "$OUT/diff/"*.dpsi 2>/dev/null)"
