#!/usr/bin/env bash
# Fetch TE/ERV annotations for the subworkflow.
#
#   LOCUS-LEVEL (Telescope): retro.hg38.v1 (HERV_rmsk + L1Base), from the
#   mlbendall telescope_annotation_db repo (Git LFS -> media.githubusercontent).
#   VERIFIED reachable: HTTP 200, 18.9 MB, 28,513 distinct loci
#   (72,169 exon-feature lines), UCSC 'chr' naming.
#
#   FAMILY-LEVEL (TEtranscripts): the curated GRCh38 GENCODE rmsk TE GTF is
#   distributed from the Hammell lab site (mghlab.org / labshare.cshl.edu).
#   That download path was NOT reachable/confirmed at authoring time, so this
#   script leaves --te_gtf_family as a parameter and prints the source. A
#   reproducible fallback from UCSC RepeatMasker is provided below (opt-in).
set -euo pipefail
REPO="$(cd "$(dirname "$0")/../../.." && pwd)"
OUT="$REPO/reference/te"; mkdir -p "$OUT"

# ---- locus-level (confirmed) ----
LOCUS_URL="https://github.com/mlbendall/telescope_annotation_db/raw/master/builds/retro.hg38.v1/transcripts.gtf"
LOCUS_OUT="$OUT/retro.hg38.v1.transcripts.gtf"
if [ ! -s "$LOCUS_OUT" ]; then
  echo "[fetch] Telescope retro.hg38.v1 -> $LOCUS_OUT"
  curl -sL -o "$LOCUS_OUT" "$LOCUS_URL"
fi
echo "[ok] locus GTF: $(grep -vc '^#' "$LOCUS_OUT") features"

# ---- family-level (fallback build from UCSC rmsk, opt-in via BUILD_FAMILY=1) ----
# The canonical file is GRCh38_GENCODE_rmsk_TE.gtf from the Hammell lab; if you
# have it, drop it at $OUT/GRCh38_GENCODE_rmsk_TE.gtf and skip this.
FAMILY_OUT="$OUT/GRCh38_rmsk_TE.gtf"
if [ "${BUILD_FAMILY:-0}" = "1" ] && [ ! -s "$FAMILY_OUT" ]; then
  echo "[fetch] UCSC hg38 RepeatMasker table (fallback family GTF)"
  curl -sL -o "$OUT/rmsk.txt.gz" "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"
  # rmsk.txt cols: bin swScore ... genoName genoStart genoEnd ... strand repName repClass repFamily
  zcat "$OUT/rmsk.txt.gz" | awk 'BEGIN{OFS="\t"}
    $12!="Simple_repeat" && $12!="Low_complexity" && $12!="Satellite" && $12!="rRNA" && $12!="tRNA" && $12!="srpRNA" && $12!="snRNA" && $12!="scRNA" {
      cls=$12; fam=$13; nm=$11;
      print $6, "rmsk", "exon", $7+1, $8, ".", $10, ".",
        "gene_id \""nm"\"; transcript_id \""nm"\"; family_id \""fam"\"; class_id \""cls"\";"
    }' > "$FAMILY_OUT"
  echo "[ok] fallback family GTF: $(grep -vc '^#' "$FAMILY_OUT") features -> $FAMILY_OUT"
fi

cat <<EOF

Locus-level (Telescope) : $LOCUS_OUT   [--te_gtf_locus]
Family-level (TEtranscripts):
  Preferred : GRCh38_GENCODE_rmsk_TE.gtf from the Hammell lab
              (https://www.mghlab.org/software/tetranscripts ->
               labshare.cshl.edu/shares/mhammelllab/... ; path not confirmed here)
  Fallback  : re-run with BUILD_FAMILY=1 to build $FAMILY_OUT from UCSC rmsk
  Pass whichever you use as  --te_gtf_family <path>
EOF
