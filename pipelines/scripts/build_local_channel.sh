#!/usr/bin/env bash
# Build a local conda channel of repackaged, sandbox-safe builds.
#
# Problem: the `referencing` conda package (a MultiQC dependency) ships an
# `etc/conda/test-files/referencing/1/suite/.git` directory. The sandbox's
# git-protection blocks creating any `.git` path, so micromamba's extraction
# aborts and corrupts the package cache ("index.json ... empty input").
#
# Fix: download the package, strip the `.git` tree from BOTH the payload tar and
# the info/paths.json manifest, repackage as a valid .conda, and serve it from a
# file:// channel prepended (highest priority) to conda.channels.
set -euo pipefail
_REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "$_REPO_ROOT/pipelines/env.sh"
ZSTD="$(command -v zstd || echo /Users/alex/.claude-science/conda/envs/nextflow/bin/zstd)"
LC="$_REPO_ROOT/.nextflow_home/local_channel"
PKG="referencing-0.37.0-pyhcf101f3_0"
URL="https://conda.anaconda.org/conda-forge/noarch/${PKG}.conda"

mkdir -p "$LC/noarch" "$LC/osx-arm64"
W="$(mktemp -d)"; trap 'rm -rf "$W"' EXIT
cd "$W"
curl -sL --max-time 120 -o orig.conda "$URL"
unzip -o orig.conda >/dev/null
"$ZSTD" -d -c info-*.tar.zst > info.tar && mkdir infoex && tar -xf info.tar -C infoex
"$ZSTD" -d -c pkg-*.tar.zst  > pkg.tar  && mkdir pkgex  && tar -xf pkg.tar  -C pkgex \
    --exclude='*/suite/.git' --exclude='*/suite/.git/*' 2>/dev/null || true

python3 - "$W" <<'PY'
import json,os,sys
W=sys.argv[1]
p=os.path.join(W,"infoex/info/paths.json")
d=json.load(open(p))
d["paths"]=[x for x in d["paths"] if os.path.exists(os.path.join(W,"pkgex",x["_path"]))]
json.dump(d,open(p,"w"))
PY

( cd infoex && tar -cf ../info.new.tar info )
( cd pkgex  && tar -cf ../pkg.new.tar  . )
"$ZSTD" -q -f -19 info.new.tar -o info-${PKG}.tar.zst
"$ZSTD" -q -f -19 pkg.new.tar  -o pkg-${PKG}.tar.zst
rm -f "${PKG}.conda"
zip -0 -q "${PKG}.conda" metadata.json info-${PKG}.tar.zst pkg-${PKG}.tar.zst
cp "${PKG}.conda" "$LC/noarch/"

python3 - "$LC" "$W" "$PKG" <<'PY'
import json,hashlib,os,sys
LC,W,PKG=sys.argv[1],sys.argv[2],sys.argv[3]
b=open(os.path.join(LC,"noarch",PKG+".conda"),"rb").read()
idx=json.load(open(os.path.join(W,"infoex/info/index.json")))
rec={"build":idx["build"],"build_number":idx.get("build_number",0),"depends":idx.get("depends",[]),
     "license":idx.get("license",""),"md5":hashlib.md5(b).hexdigest(),"name":"referencing","noarch":"python",
     "sha256":hashlib.sha256(b).hexdigest(),"size":len(b),"subdir":"noarch","timestamp":idx.get("timestamp",0),"version":"0.37.0"}
json.dump({"info":{"subdir":"noarch"},"packages":{},"packages.conda":{PKG+".conda":rec},"repodata_version":1},
          open(os.path.join(LC,"noarch","repodata.json"),"w"),indent=1)
json.dump({"info":{"subdir":"osx-arm64"},"packages":{},"packages.conda":{},"repodata_version":1},
          open(os.path.join(LC,"osx-arm64","repodata.json"),"w"))
print("local channel built:",LC)
PY
