#!/usr/bin/env bash
# Full-depth salmon quant with robust resumable download + MD5 verify.
# Same quant recipe as run_salmon_pilot.sh (salmon -l A --gcBias --seqBias -g tx2gene),
# but downloads the FULL FASTQ pair (no subsample) with:
#   - resume on truncation (curl -C -) up to N attempts
#   - exact Content-Length size check before accepting
#   - MD5 verification against the manifest
# Stream-align-delete: peak disk ~= one sample's FASTQ pair.
# Usage: MANIFEST INDEX OUTDIR   (env: SALMON THREADS GENEMAP)
set -uo pipefail
MANIFEST="$1"; INDEX="$2"; OUT="$3"
SALMON="${SALMON:-salmon}"; THREADS="${THREADS:-8}"; GENEMAP="${GENEMAP:-}"
TMP="$OUT/_fastq_tmp"; mkdir -p "$TMP"
LOG="$OUT/fulldepth.log"; : > "$LOG"
[ -f "$OUT/pilot_index.csv" ] || echo "run_accession,status,mapping_rate,n_processed,quant_path" > "$OUT/pilot_index.csv"
md5of() { if command -v md5sum >/dev/null; then md5sum "$1" | awk '{print $1}'; else md5 -q "$1"; fi; }

# resumable download to exact size; returns 0 on full+correct size
robust_get() {
  local url="$1" dest="$2" want="$3" tries=0
  while [ $tries -lt 8 ]; do
    tries=$((tries+1))
    curl -sS -C - --retry 5 --retry-delay 3 --retry-all-errors \
         --max-time 3600 --speed-time 120 --speed-limit 2000 \
         -o "$dest" "$url" >/dev/null 2>>"$LOG" || true
    local have; have=$(stat -f%z "$dest" 2>/dev/null || stat -c%s "$dest" 2>/dev/null || echo 0)
    if [ -n "$want" ] && [ "$have" = "$want" ]; then return 0; fi
    if [ -z "$want" ] && [ "$have" -gt 0 ]; then
      # no known size: accept if a second attempt adds nothing
      curl -sS -C - --retry 3 --retry-all-errors --max-time 3600 -o "$dest" "$url" >/dev/null 2>>"$LOG" || true
      local have2; have2=$(stat -f%z "$dest" 2>/dev/null || stat -c%s "$dest" 2>/dev/null || echo 0)
      [ "$have2" = "$have" ] && return 0
    fi
    echo "  [retry $tries] $dest have=$have want=${want:-?}" | tee -a "$LOG"
  done
  return 1
}
content_len() { curl -sSI --max-time 120 "$1" 2>/dev/null | awk 'BEGIN{IGNORECASE=1}/^content-length:/{v=$2} END{gsub(/\r/,"",v); print v}'; }

# disk guard: free MB on the results filesystem
free_mb() { df -m "$OUT" | tail -1 | awk '{print $4}'; }
MIN_FREE_MB="${MIN_FREE_MB:-30000}"   # yield to sibling jobs below this

while IFS=, read -r cohort run sample title patient arm tp recist resp gb ftp md5 rest <&3; do
  sdir="$OUT/$run"
  if [ -s "$sdir/quant.sf" ]; then echo "  [$run] exists, skip" | tee -a "$LOG"; continue; fi
  # wait for disk headroom before starting a sample (peak ~10-11GB/pair)
  waited=0
  while [ "$(free_mb)" -lt "$MIN_FREE_MB" ]; do
    echo "  [disk-guard] free $(free_mb)MB < ${MIN_FREE_MB}MB — waiting 120s (waited ${waited}s)" | tee -a "$LOG"
    sleep 120; waited=$((waited+120))
    if [ "$waited" -ge 7200 ]; then echo "  [disk-guard] gave up after 2h; stopping" | tee -a "$LOG"; break 2; fi
  done
  echo "==== [$run] ($resp) full-depth ====" | tee -a "$LOG"
  IFS=';' read -r u1 u2 <<< "$(echo "$ftp" | tr -d '\"')"
  IFS=';' read -r e1 e2 <<< "$(echo "$md5" | tr -d '\"')"
  u1="https://${u1#ftp://}"; u1="${u1/https:\/\/https:\/\//https://}"
  u2="https://${u2#ftp://}"; u2="${u2/https:\/\/https:\/\//https://}"
  r1="$TMP/${run}_1.fastq.gz"; r2="$TMP/${run}_2.fastq.gz"
  ok=1
  for pair in "$u1|$r1|$e1" "$u2|$r2|$e2"; do
    IFS='|' read -r url dest expmd5 <<< "$pair"
    echo "  downloading (full) $url" | tee -a "$LOG"
    want=$(content_len "$url")
    if ! robust_get "$url" "$dest" "$want"; then echo "  DL FAIL $dest" | tee -a "$LOG"; ok=0; break; fi
    got=$(md5of "$dest")
    if [ -n "$expmd5" ] && [ "$got" != "$expmd5" ]; then
      echo "  MD5 MISMATCH $dest exp=$expmd5 got=$got" | tee -a "$LOG"; ok=0; break
    fi
    echo "  ok $dest ($(stat -f%z "$dest" 2>/dev/null || stat -c%s "$dest") bytes, md5 verified)" | tee -a "$LOG"
  done
  if [ "$ok" -ne 1 ]; then echo "$run,dl_or_md5_fail,,," >> "$OUT/pilot_index.csv"; rm -f "$r1" "$r2"; continue; fi
  GENEMAP_ARG=""; [ -n "$GENEMAP" ] && [ -f "$GENEMAP" ] && GENEMAP_ARG="-g $GENEMAP"
  "$SALMON" quant -i "$INDEX" -l A -1 "$r1" -2 "$r2" -p "$THREADS" --gcBias --seqBias $GENEMAP_ARG \
      -o "$sdir" </dev/null >> "$LOG" 2>&1 || { echo "$run,quant_fail,,," >> "$OUT/pilot_index.csv"; rm -f "$r1" "$r2"; continue; }
  meta="$sdir/aux_info/meta_info.json"
  mr=$(grep -o '"percent_mapped"[^,]*' "$meta" 2>/dev/null | grep -o '[0-9.]*' | head -1 || echo "")
  np=$(grep -o '"num_processed"[^,]*' "$meta" 2>/dev/null | grep -o '[0-9]*' | head -1 || echo "")
  echo "$run,ok,${mr:-NA},${np:-NA},$sdir/quant.sf" >> "$OUT/pilot_index.csv"
  echo "  [$run] done mapping_rate=${mr:-NA}% n=${np:-NA}" | tee -a "$LOG"
  rm -f "$r1" "$r2"
done 3< <(tail -n +2 "$MANIFEST")
rmdir "$TMP" 2>/dev/null || true
echo "ALL DONE" | tee -a "$LOG"
