#!/usr/bin/env bash
# Disk-bounded cohort processing: fetch a BATCH of samples -> run rnaseq spine ->
# keep BAMs+counts -> delete that batch's FASTQs -> next batch.
# Keeps peak FASTQ footprint ~ (batch_size * ~6 GB) instead of the whole cohort.
#
# The curated subset is 40 samples / ~219 GB (fits 697 GB free), so batching is
# a SAFETY MARGIN, not a hard requirement — set BATCH large to process all at once.
# It becomes essential if you later add Gide-full (~824 GB) + Riaz-full (~486 GB).
#
# Usage: run_cohort_batched.sh <manifest.csv|json> <outdir> [batch_size=8]
set -euo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
source "$REPO/pipelines/env.sh" >/dev/null
MAN="${1:?manifest}"; OUT="${2:?outdir}"; BATCH="${3:-8}"
FASTQ_DIR="$REPO/data/raw/fastq"; mkdir -p "$FASTQ_DIR" "$OUT"
# reference FASTA for optional archival CRAM (step 5, ARCHIVE_CRAM=1)
REF_FASTA="${REF_FASTA:-$REPO/reference/GRCh38/GRCh38.primary_assembly.genome.fa}"

# Split manifest into batches of run accessions
python3 - "$MAN" "$BATCH" <<'PY'
import sys,json,csv,os
m,bs=sys.argv[1],int(sys.argv[2])
rows=json.load(open(m)) if m.endswith(".json") else list(csv.DictReader(open(m)))
os.makedirs("/tmp/batches",exist_ok=True)
for i in range(0,len(rows),bs):
    chunk=rows[i:i+bs]
    with open(f"/tmp/batches/batch_{i//bs:03d}.csv","w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(chunk)
print(f"{(len(rows)+bs-1)//bs} batches of up to {bs} from {len(rows)} samples")
PY

for bfile in /tmp/batches/batch_*.csv; do
  bname="$(basename "$bfile" .csv)"
  echo "===== $bname ====="
  # 1. fetch this batch's FASTQs (md5-verified)
  bash "$REPO/pipelines/scripts/fetch_fastq.sh" "$bfile" "$FASTQ_DIR"
  # 2. build a samplesheet for just this batch
  python3 "$REPO/pipelines/scripts/_ss_from_manifest.py" "$bfile" "$FASTQ_DIR" "/tmp/ss_${bname}.csv"
  # 3. run the rnaseq spine on the batch (resumable; index reused)
  bash "$REPO/pipelines/scripts/run_rnaseq.sh" "/tmp/ss_${bname}.csv" "$OUT"
  # 4. delete this batch's FASTQs (BAMs+counts already published under $OUT)
  awk -F, 'NR>1{print $1}' /tmp/batches/${bname}.csv >/dev/null 2>&1 || true
  python3 - "$bfile" "$FASTQ_DIR" <<'PY'
import sys,json,csv,os,glob
m,fq=sys.argv[1],sys.argv[2]
rows=list(csv.DictReader(open(m)))
for r in rows:
    for f in glob.glob(f"{fq}/{r['run_accession']}_*.fastq.gz"):
        os.remove(f); print("[del]",os.path.basename(f))
PY
  # 5. OPTIONAL archival CRAM: if ARCHIVE_CRAM=1, losslessly convert this batch's
  #    spine BAMs to CRAM and delete the BAMs (reclaims ~50%). Leave UNSET when the
  #    custom subworkflows (te_erv/intron_retention/rna_editing) still need to read
  #    the BAMs — run them first, then re-invoke with ARCHIVE_CRAM=1, or archive by
  #    hand via pipelines/scripts/archive_bam_to_cram.sh. See RUNBOOK §8.
  if [ "${ARCHIVE_CRAM:-0}" = "1" ]; then
    echo "  [cram] archiving batch BAMs -> CRAM (ARCHIVE_CRAM=1)"
    bash "$REPO/pipelines/scripts/archive_bam_to_cram.sh" \
      "$REF_FASTA" "$OUT/hisat2" -j "${CRAM_THREADS:-6}" || echo "  [cram] WARN: archival reported failures (BAMs kept)"
  fi
  echo "  free after $bname: $(df -h "$REPO" | awk 'NR==2{print $4}')"
done
echo "[done] cohort processed; FASTQs cleared, BAMs+counts under $OUT"
[ "${ARCHIVE_CRAM:-0}" = "1" ] && echo "[done] spine BAMs archived to CRAM (ARCHIVE_CRAM=1)"
