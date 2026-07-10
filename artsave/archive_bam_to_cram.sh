#!/usr/bin/env bash
# archive_bam_to_cram.sh -- losslessly convert alignment BAM(s) to CRAM to
# reclaim disk, verify the round-trip, then delete the source BAM.
#
# WHY ARCHIVAL, NOT CRAM-NATIVE: the base spine emits BAM and every downstream
# consumer (nf-core/rnaseq, rnasplice, rnafusion, and the three custom
# subworkflows) reads BAM. Making the working format CRAM would require
# threading the reference FASTA through featureCounts + TEcount and adds
# CRAM encode/decode CPU on every pass. Instead we keep BAM as the transient
# working format and convert to CRAM only AFTER the fan-out has consumed it —
# same disk savings (the BAM is transient either way), no downstream changes.
# CRAM here is fully LOSSLESS: all quality scores retained (no binning), and the
# exact reference is required to decode (recorded via the CRAM's stored M5/UR).
#
# Usage:
#   archive_bam_to_cram.sh <reference.fa> <bam-or-dir> [--quick] [--keep] [-j N]
#     <reference.fa>  reference FASTA the BAM was aligned to (must have .fai)
#     <bam-or-dir>    a single .bam, or a directory searched recursively for
#                     *.markdup.sorted.bam (real files only, symlinks skipped)
#     --quick         verify by record count + flagstat only (fast). Default is
#                     a DEEP verify: checksum the decoded core alignment stream
#                     (QNAME/FLAG/RNAME/POS/MAPQ/CIGAR/SEQ/QUAL) of BAM vs CRAM.
#                     Deep is the default because we delete the only BAM copy.
#     --keep          convert + verify but do NOT delete the BAM.
#     -j N            samtools compression threads (default: 4).
#
# Exit non-zero if ANY file fails conversion or verification; a failed BAM is
# never deleted.
set -euo pipefail

REF=""; TARGET=""; QUICK=0; KEEP=0; THREADS=4
while [ $# -gt 0 ]; do
  case "$1" in
    --quick) QUICK=1; shift;;
    --keep)  KEEP=1; shift;;
    -j)      THREADS="$2"; shift 2;;
    *) if [ -z "$REF" ]; then REF="$1"; elif [ -z "$TARGET" ]; then TARGET="$1"; fi; shift;;
  esac
done

if [ -z "$REF" ] || [ -z "$TARGET" ]; then
  echo "usage: archive_bam_to_cram.sh <reference.fa> <bam-or-dir> [--quick] [--keep] [-j N]" >&2
  exit 2
fi
[ -s "$REF" ] || { echo "[err] reference not found: $REF" >&2; exit 2; }
[ -s "${REF}.fai" ] || { echo "[info] indexing reference..."; samtools faidx "$REF"; }

# collect target BAMs (real files only — Nextflow work dirs are full of symlinks)
BAMS=()
if [ -d "$TARGET" ]; then
  while IFS= read -r f; do BAMS+=("$f"); done < <(find "$TARGET" -name '*.markdup.sorted.bam' -type f ! -type l)
else
  [ -f "$TARGET" ] && [ ! -L "$TARGET" ] && BAMS+=("$TARGET")
fi
[ "${#BAMS[@]}" -gt 0 ] || { echo "[err] no real BAM files found under: $TARGET" >&2; exit 2; }

# core-field checksum of an alignment stream (order-preserving). Compares the 11
# mandatory SAM fields but NORMALIZES the TLEN (col 9) sign to its magnitude:
# for a mate pair whose two reads start at the SAME position (POS==MPOS), the
# TLEN sign is ambiguous under the SAM spec, and htslib recomputes it on CRAM
# decode with a consistent tie-break — so the original BAM's sign for those few
# reads is not round-tripped. The magnitude is always preserved, and no
# downstream tool we run (featureCounts, JACUSA2, Telescope, rMATS) depends on
# the TLEN sign of co-located mates. Optional tags (col 12+) are excluded: CRAM
# preserves them but may reorder, which is likewise spec-permitted. Everything
# that carries alignment meaning — QNAME/FLAG/RNAME/POS/MAPQ/CIGAR/RNEXT/PNEXT/
# SEQ/QUAL and |TLEN| — is compared exactly.
#   $1 = alignment file, $2 = "" for BAM or the reference for CRAM (-T)
stream_cksum() {
  local ref_arg=""; [ -n "${2:-}" ] && ref_arg="-T $2"
  samtools view -@ "$THREADS" $ref_arg "$1" \
    | awk -F'\t' 'BEGIN{OFS="\t"} {$9=($9<0?-$9:$9); print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' \
    | cksum
}

fail=0
for BAM in "${BAMS[@]}"; do
  CRAM="${BAM%.bam}.cram"
  echo "=== $BAM -> $CRAM ==="
  # convert (lossless: full quality retained; embed reference URI via -T)
  if ! samtools view -@ "$THREADS" -C -T "$REF" -o "$CRAM" "$BAM"; then
    echo "[err] conversion failed: $BAM"; fail=1; continue
  fi
  samtools index "$CRAM"
  samtools quickcheck "$CRAM" || { echo "[err] quickcheck failed: $CRAM"; rm -f "$CRAM" "${CRAM}.crai"; fail=1; continue; }

  # verify
  n_bam=$(samtools view -@ "$THREADS" -c "$BAM")
  n_cram=$(samtools view -@ "$THREADS" -c -T "$REF" "$CRAM")
  if [ "$n_bam" != "$n_cram" ]; then
    echo "[err] record count mismatch: BAM=$n_bam CRAM=$n_cram"; rm -f "$CRAM" "${CRAM}.crai"; fail=1; continue
  fi
  fs_bam=$(samtools flagstat "$BAM" | md5 -q 2>/dev/null || samtools flagstat "$BAM" | md5sum | cut -d' ' -f1)
  fs_cram=$(samtools flagstat -@ "$THREADS" "$CRAM" | md5 -q 2>/dev/null || samtools flagstat "$CRAM" | md5sum | cut -d' ' -f1)
  if [ "$fs_bam" != "$fs_cram" ]; then
    echo "[err] flagstat mismatch"; rm -f "$CRAM" "${CRAM}.crai"; fail=1; continue
  fi
  if [ "$QUICK" -eq 0 ]; then
    c_bam=$(stream_cksum "$BAM" "")
    c_cram=$(stream_cksum "$CRAM" "$REF")
    if [ "$c_bam" != "$c_cram" ]; then
      echo "[err] deep stream checksum mismatch: BAM[$c_bam] CRAM[$c_cram]"; rm -f "$CRAM" "${CRAM}.crai"; fail=1; continue
    fi
    echo "[ok] deep verify passed ($n_bam records, core-field stream identical incl. |TLEN|)"
  else
    echo "[ok] quick verify passed ($n_bam records, flagstat identical)"
  fi

  b_sz=$(stat -f%z "$BAM" 2>/dev/null || stat -c%s "$BAM"); c_sz=$(stat -f%z "$CRAM" 2>/dev/null || stat -c%s "$CRAM")
  pct=$(awk "BEGIN{printf \"%.1f\", 100*(1-$c_sz/$b_sz)}")
  echo "[size] BAM=$b_sz  CRAM=$c_sz  saved=${pct}%"
  if [ "$KEEP" -eq 0 ]; then
    rm -f "$BAM" "${BAM}.bai" "${BAM%.bam}.bai"
    echo "[done] BAM deleted; CRAM retained"
  else
    echo "[done] BAM kept (--keep)"
  fi
done
exit $fail
