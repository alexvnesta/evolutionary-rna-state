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
import json,os,subprocess,hashlib
urls=json.load(open("/tmp/_fastq_urls.json"))
def md5(p):
    h=hashlib.md5()
    with open(p,'rb') as fh:
        for c in iter(lambda: fh.read(1<<20), b''): h.update(c)
    return h.hexdigest()
for u,h in urls:
    fn=os.path.basename(u); url="https://"+u if not u.startswith("http") else u
    if os.path.exists(fn) and md5(fn)==h:
        print(f"[ok  ] {fn} (verified)"); continue
    print(f"[get ] {fn}")
    subprocess.run(["curl","-fsSL","--retry","3","-o",fn+".part",url],check=True)
    got=md5(fn+".part")
    if got!=h:
        os.remove(fn+".part"); raise SystemExit(f"MD5 MISMATCH {fn}: {got} != {h}")
    os.rename(fn+".part",fn); print(f"[ok  ] {fn} (md5 {got[:8]}…)")
print("[done] all fastq verified")
PY
