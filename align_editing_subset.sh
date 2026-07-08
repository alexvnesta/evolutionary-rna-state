#!/usr/bin/env bash
# Download + STAR-align a balanced subset for A-to-I editing analysis.
# Serial (one STAR load at a time) to bound RAM; deletes FASTQ after each align to bound disk.
set -uo pipefail
REPO=/Users/alex/OrchestratedBiosciences/evolutionary-rna-state
IDX=$REPO/reference/GRCh38/star_index
OUT=$REPO/results/editing_bams
FQ=$REPO/data/raw/fastq
MAN=$REPO/editing_subset_manifest.csv
mkdir -p "$OUT" "$FQ"
LOG=$OUT/align.log
echo "=== align run started $(date) ===" > "$LOG"

# parse manifest: sample, ... , fastq_ftp (last col). fastq_ftp is 'url1;url2'
tail -n +2 "$MAN" | while IFS=, read -r line; do
  samp=$(echo "$line" | cut -d, -f1)
  ftp=$(echo "$line" | awk -F, '{print $NF}')
  bam="$OUT/${samp}.Aligned.sortedByCoord.out.bam"
  if [[ -f "$bam" ]]; then echo "[$samp] BAM exists, skip" >>"$LOG"; continue; fi
  u1=$(echo "$ftp" | cut -d';' -f1); u2=$(echo "$ftp" | cut -d';' -f2)
  f1="$FQ/${samp}_1.fastq.gz"; f2="$FQ/${samp}_2.fastq.gz"
  echo "[$samp] downloading $(date)" >>"$LOG"
  [[ -f "$f1" ]] || curl -sS -L --retry 3 "https://${u1}" -o "$f1" 2>>"$LOG"
  [[ -f "$f2" ]] || curl -sS -L --retry 3 "https://${u2}" -o "$f2" 2>>"$LOG"
  if [[ ! -s "$f1" || ! -s "$f2" ]]; then echo "[$samp] DOWNLOAD FAILED" >>"$LOG"; continue; fi
  echo "[$samp] STAR $(date)" >>"$LOG"
  STAR --runThreadN 14 --genomeDir "$IDX" \
       --readFilesIn "$f1" "$f2" --readFilesCommand gunzip -c \
       --outSAMtype BAM SortedByCoordinate \
       --outFileNamePrefix "$OUT/${samp}." \
       --outSAMattributes NH HI AS nM MD \
       --outFilterMultimapNmax 1 \
       --outSAMprimaryFlag AllBestScore \
       --limitBAMsortRAM 20000000000 >>"$LOG" 2>&1
  if [[ -f "$bam" ]]; then
     samtools index "$bam"
     echo "[$samp] DONE $(date) $(samtools view -c "$bam" 2>/dev/null) reads" >>"$LOG"
     rm -f "$f1" "$f2"
     rm -rf "$OUT/${samp}._STARtmp" "$OUT/${samp}.Aligned.out.sam"
  else
     echo "[$samp] STAR FAILED" >>"$LOG"
  fi
done
echo "=== align run finished $(date) ===" >>"$LOG"
