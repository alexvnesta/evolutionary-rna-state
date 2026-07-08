#!/usr/bin/env bash
# Fetch GRCh38 primary assembly + GENCODE v46 annotation to reference/ (git-ignored).
# Records sizes + md5 to reference/GRCh38/checksums.json. Idempotent (skips existing).
set -euo pipefail
REL=46
GEN="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${REL}"
DEST="$(cd "$(dirname "$0")/../../reference/GRCh38" && pwd)"
FA="GRCh38.primary_assembly.genome.fa.gz"
GTF="gencode.v${REL}.primary_assembly.annotation.gtf.gz"
cd "$DEST"
for pair in "$FA:$GEN/$FA" "$GTF:$GEN/$GTF"; do
  f="${pair%%:*}"; url="${pair#*:}"
  if [ -s "$f" ]; then echo "[skip] $f exists ($(du -h "$f"|cut -f1))"; else
    echo "[get ] $f"; curl -fsSL --retry 3 -o "$f.part" "$url" && mv "$f.part" "$f"; fi
done
echo "[gunzip] keeping .gz + producing uncompressed for STAR"
[ -s "${FA%.gz}" ]  || gunzip -kf "$FA"
[ -s "${GTF%.gz}" ] || gunzip -kf "$GTF"
python3 - <<'PY'
import hashlib,json,os,glob
def md5(p):
    h=hashlib.md5()
    with open(p,'rb') as fh:
        for c in iter(lambda: fh.read(1<<20), b''): h.update(c)
    return h.hexdigest()
rec={}
for p in sorted(glob.glob("*.fa*")+glob.glob("*.gtf*")):
    rec[p]={"bytes":os.path.getsize(p),"md5":md5(p)}
json.dump(rec, open("checksums.json","w"), indent=2)
print(json.dumps(rec, indent=2))
PY
echo "[done] reference staged at $DEST"
