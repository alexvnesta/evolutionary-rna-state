#!/usr/bin/env bash
# pipelines/bcr_repertoire/run_trust4_pilot.sh
# ----------------------------------------------------------------------------
# BCR/TCR repertoire reconstruction arm of the raw-read pipeline. Mirrors
# analysis/pilot/run_salmon_pilot.sh's stream-align-delete pattern: for each run
# in the manifest, stream N read-pairs from ENA, run TRUST4, keep only the small
# report files (_cdr3.out/_report.tsv/_airr.tsv), delete FASTQ. Peak disk ~1
# sample. The downstream feature extraction lives in
# analysis/differentiated/bcr_shm.py (build_bcr_features on this OUT dir).
#
# TRUST4 (Song et al., Nat Methods 2021, doi:10.1038/s41592-021-01142-2) is
# built from source (tools/TRUST4-1.1.5, arm64-clean) — the osx-arm64 conda
# build path that broke STAR/SNAF is avoided.
#
# DEPTH NOTE: BCR reads are a tiny fraction of a bulk library, so this arm
# streams a LARGER subsample (default 15M pairs) than the 3M the gene-TPM pilot
# uses. Contig yield/SHM stability still scale with depth and B-cell content;
# bcr_n_clonotypes/bcr_n_reads are emitted per sample so reliability is visible.
#
# Usage: run_trust4_pilot.sh <manifest.csv> <out_dir>
#   manifest columns (superset OK): run_accession, cohort, fastq_ftp (mate URLs
#   ';'-separated), and any others (ignored). Env: NREADS, THREADS, T4DIR.
set -uo pipefail

MANIFEST="${1:?manifest csv (needs run_accession,cohort,fastq_ftp columns)}"
OUT="${2:?output dir}"
T4DIR="${T4DIR:-tools/TRUST4-1.1.5}"
TRUST4="$T4DIR/run-trust4"
REFCOORD="$T4DIR/hg38_bcrtcr.fa"
REFIMGT="$T4DIR/human_IMGT+C.fa"
NREADS="${NREADS:-15000000}"
THREADS="${THREADS:-8}"

for f in "$TRUST4" "$REFCOORD" "$REFIMGT"; do
  [ -e "$f" ] || { echo "FATAL: missing $f (build TRUST4: extract tools/trust4_built_arm64.tar.gz)"; exit 1; }
done
mkdir -p "$OUT" "$OUT/_fastq_tmp"
LOG="$OUT/trust4_run.log"; : > "$LOG"
IDX="$OUT/trust4_index.csv"
echo "run_accession,cohort,status,n_igh_contigs,trust4_seconds" > "$IDX"

# Resolve the header -> column positions (robust to column order / extra cols).
header=$(head -1 "$MANIFEST")
col() { echo "$header" | tr ',' '\n' | grep -n -x "$1" | head -1 | cut -d: -f1; }
C_RUN=$(col run_accession); C_COH=$(col cohort); C_FTP=$(col fastq_ftp)
[ -n "$C_RUN" ] && [ -n "$C_FTP" ] || { echo "FATAL: manifest needs run_accession + fastq_ftp columns"; exit 1; }

# Stream one mate: curl | zcat | head N*4 | gzip. Judge success by OUTPUT
# integrity (valid gzip, enough reads), not the pipe exit (head SIGPIPEs curl).
stream_subsample() {
  local url="$1" dest="$2" nreads="$3"
  case "$url" in http*) : ;; *) url="https://$url" ;; esac
  ( set +o pipefail
    curl -sS --retry 3 --retry-delay 2 "$url" | zcat 2>/dev/null | head -n $((nreads*4)) | gzip > "$dest" )
  gzip -t "$dest" 2>/dev/null || return 1
  local n; n=$(zcat < "$dest" 2>/dev/null | wc -l); n=$((n/4))
  [ "$n" -ge $((nreads/2)) ] || return 2   # accept if >= half (sample may be shallower than N)
  return 0
}

# Manifest body on FD 3, NOT stdin, so curl inside the loop can't eat lines.
tail -n +2 "$MANIFEST" > "$OUT/_manifest_body.csv"
while IFS=, read -r -a F <&3; do
  run="${F[$((C_RUN-1))]}"; run="${run//\"/}"
  [ -z "$run" ] && continue
  cohort="${F[$((C_COH-1))]:-NA}"; cohort="${cohort//\"/}"
  ftp="${F[$((C_FTP-1))]}"; ftp="${ftp//\"/}"
  odir="$OUT/$run"
  if [ -f "$odir/${run}_cdr3.out" ]; then echo "[$run] done — skip" | tee -a "$LOG"; continue; fi
  echo "==== [$run] $cohort ====" | tee -a "$LOG"
  IFS=';' read -r u1 u2 <<< "$ftp"
  r1="$OUT/_fastq_tmp/${run}_1.fq.gz"; r2="$OUT/_fastq_tmp/${run}_2.fq.gz"
  ok=1
  for pair in "$u1|$r1" "$u2|$r2"; do
    IFS='|' read -r url dest <<< "$pair"
    [ -z "$url" ] && continue
    stream_subsample "$url" "$dest" "$NREADS" || { echo "  STREAM FAIL rc=$? $url" | tee -a "$LOG"; ok=0; break; }
  done
  if [ "$ok" -ne 1 ]; then rm -f "$r1" "$r2"; echo "$run,$cohort,stream_fail,,"  >> "$IDX"; continue; fi
  mkdir -p "$odir"; t0=$(date +%s)
  if [ -n "${u2:-}" ]; then
    "$TRUST4" -1 "$r1" -2 "$r2" -f "$REFCOORD" --ref "$REFIMGT" -t "$THREADS" -o "$run" --od "$odir" >> "$LOG" 2>&1
  else
    "$TRUST4" -u "$r1" -f "$REFCOORD" --ref "$REFIMGT" -t "$THREADS" -o "$run" --od "$odir" >> "$LOG" 2>&1
  fi
  rc4=$?; t1=$(date +%s)
  nigh=$(awk -F'\t' '$3 ~ /^IGHV/ || $6 ~ /^IGH[MDGAE]/ {c++} END{print c+0}' "$odir/${run}_cdr3.out" 2>/dev/null)
  echo "$run,$cohort,$([ $rc4 -eq 0 ] && echo ok || echo trust4_fail),$nigh,$((t1-t0))" >> "$IDX"
  rm -f "$r1" "$r2"   # stream-align-delete
  echo "  [$run] exit=$rc4 ${nigh} IGH contigs in $((t1-t0))s" | tee -a "$LOG"
done 3< "$OUT/_manifest_body.csv"
rm -f "$OUT/_manifest_body.csv"
echo "TRUST4 PILOT COMPLETE — index at $IDX" | tee -a "$LOG"
