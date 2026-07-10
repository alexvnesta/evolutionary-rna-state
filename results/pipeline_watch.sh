#!/usr/bin/env bash
# Autonomous watcher: canonical Nextflow cohort build + legacy nonref_run drain.
# Logs one status line per cycle; flags REAL errors (not transient SSL retries).
# Written by audit session 64079601, 2026-07-10. Safe/read-only: never kills or deletes.
set -uo pipefail
REPO="/Users/alex/OrchestratedBiosciences/evolutionary-rna-state"
cd "$REPO"
RUN="results/nonref_run"
WLOG="results/pipeline_watch.log"
CYCLES="${CYCLES:-240}"      # 240 * 60s = 4h default
SLEEP="${SLEEP:-60}"
ts(){ date '+%Y-%m-%d %H:%M:%S'; }
echo "[$(ts)] watcher start (canonical=Nextflow cohort, draining=nonref_run)" >> "$WLOG"
for i in $(seq 1 "$CYCLES"); do
  # --- canonical Nextflow cohort ---
  spine=$(ls results/rnaseq_cohort/hisat2/*.markdup.sorted.bam 2>/dev/null | wc -l | tr -d ' ')
  nfrun=$(lsof 2>/dev/null | grep -ciE 'nf_work|\.nextflow/')
  fetching=$(ls data/raw/fastq/*.part 2>/dev/null | wc -l | tr -d ' ')
  part=$(ls -la data/raw/fastq/*.part 2>/dev/null | awk '{s+=$5} END{printf "%.1fG", s/1e9}')
  # real (non-transient) errors in the cohort resume log: exclude the retry-recover lines
  cohorterr=$(grep -aiE "error|fail|exit status [1-9]|Command error" results/rnaseq_cohort_resume.log 2>/dev/null \
              | grep -viE "SSL_read|unexpected eof|Problem \(retrying|retry|Warning|retrying in|--retry" | tail -1)
  # --- legacy drain ---
  inflight=$(ls "$RUN"/work/ 2>/dev/null | tr '\n' ' ')
  legacyproc=$(lsof +D "$RUN" 2>/dev/null | awk 'NR>1{print $2}' | sort -u | wc -l | tr -d ' ')
  orchn=$(tail -c 300 "$RUN"/orchestrate_run.log 2>/dev/null | tr '\r' '\n' | grep -oE '[0-9]+/105' | tail -1)
  load=$(uptime | sed -E 's/.*load averages?: //' | awk '{print $1}')
  free=$(df -h "$REPO" | awk 'NR==2{print $4}')
  drained="no"; [ -z "$inflight" ] && [ "$legacyproc" -le 1 ] && drained="YES"
  echo "[$(ts)] canon: spine=$spine/40 nf_active=$nfrun fetch=${fetching}(${part}) | legacy: inflight=[$inflight] procs=$legacyproc orch=$orchn drained=$drained | load=$load free=$free" >> "$WLOG"
  [ -n "$cohorterr" ] && echo "[$(ts)] ⚠ CANON-ERR: $cohorterr" >> "$WLOG"
  # stop early once legacy is fully drained AND canonical clearly advancing
  if [ "$drained" = "YES" ]; then
    echo "[$(ts)] ✓ legacy nonref_run fully drained; canonical has the machine. watcher continuing to monitor canonical." >> "$WLOG"
  fi
  sleep "$SLEEP"
done
echo "[$(ts)] watcher exit after $CYCLES cycles" >> "$WLOG"