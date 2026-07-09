#!/usr/bin/env bash
# Fetch paired FASTQs for a set of runs from ENA, verifying md5 against the manifest.
# Usage: fetch_fastq.sh <manifest.csv|pilot_manifest.json> [dest_dir]
# Reads ENA FTP URLs + md5 from the manifest (fastq_ftp / fastq_md5 fields).
set -euo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
MAN="${1:?manifest required}"
DEST="${2:-$REPO/data/raw/fastq}"
mkdir -p "$DEST"; cd "$DEST"
python3 - "$MAN" <<'PY'
import sys,json,csv,os
m=sys.argv[1]
rows=json.load(open(m)) if m.endswith(".json") else list(csv.DictReader(open(m)))
out=[]
for r in rows:
    ftp=r["fastq_ftp"].split(";"); md5=r["fastq_md5"].split(";")
    for u,h in zip(ftp,md5): out.append((u.strip(),h.strip()))
json.dump(out, open("/tmp/_fastq_urls.json","w"))
print(f"{len(out)} files across {len(rows)} runs")
PY
python3 - <<'PY'
import json,os,subprocess,hashlib,time
urls=json.load(open("/tmp/_fastq_urls.json"))
def md5(p):
    h=hashlib.md5()
    with open(p,'rb') as fh:
        for c in iter(lambda: fh.read(1<<20), b''): h.update(c)
    return h.hexdigest()

# Resilient per-file fetch: ENA/curl can drop mid-stream ("SSL_read: unexpected
# eof"), which must NOT abort a multi-hour cohort fetch. Resume the .part file
# with `curl -C -`, retry the download several times with backoff, and only
# fail (raise) if md5 is still wrong after all attempts.
MAX_TRIES=6
def fetch_one(url, part):
    for attempt in range(1, MAX_TRIES+1):
        # -C - resumes a partial .part; --retry/--retry-all-errors handle
        # transient HTTP/TLS errors within a single curl invocation too.
        r=subprocess.run(
            ["curl","-fSL","-C","-","--retry","5","--retry-delay","5",
             "--retry-all-errors","--connect-timeout","60","-o",part,url])
        if r.returncode==0:
            return True
        wait=min(60, 5*attempt)
        print(f"        curl exit {r.returncode} (attempt {attempt}/{MAX_TRIES}); retrying in {wait}s")
        time.sleep(wait)
    return False

failed=[]
for u,h in urls:
    fn=os.path.basename(u); url="https://"+u if not u.startswith("http") else u
    part=fn+".part"
    if os.path.exists(fn) and md5(fn)==h:
        print(f"[ok  ] {fn} (verified)"); continue
    print(f"[get ] {fn}")
    ok=False
    for md5_try in range(1,3):  # re-download from scratch once if md5 mismatches
        if not fetch_one(url, part):
            break
        got=md5(part)
        if got==h:
            os.rename(part,fn); print(f"[ok  ] {fn} (md5 {got[:8]}…)"); ok=True; break
        print(f"        MD5 mismatch {fn}: {got[:8]} != {h[:8]} (try {md5_try}); re-fetching clean")
        if os.path.exists(part): os.remove(part)  # force clean re-download
    if not ok:
        failed.append(fn); print(f"[FAIL] {fn} after retries")
if failed:
    raise SystemExit(f"[error] {len(failed)} file(s) failed to fetch: {failed}")
print("[done] all fastq verified")
PY
