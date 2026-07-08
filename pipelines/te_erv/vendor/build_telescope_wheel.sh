#!/usr/bin/env bash
# build_telescope_wheel.sh -- build an arm64 wheel of Telescope v1.0.3 that the
# te_erv environment.yml installs (see the pip: section there).
#
# WHY A PRE-BUILT WHEEL, not a pip source/git install:
#   1. git+https  -> sandbox forbids .git dir creation ("Operation not permitted")
#   2. remote URL -> conda-env pip can't verify the sandbox proxy TLS cert
#                    (SSLCertVerificationError OSStatus -26276)
#   3. source build -> Telescope v1.0.3's calignment.pyx does
#        `from calignment cimport AlignedPair`
#      i.e. the module cimports ITSELF; modern Cython (0.29.36) rejects that
#      (the sibling calignment.pxd is auto-applied, so the line is redundant).
#      We strip that one line, then the single Cython extension compiles cleanly
#      against the conda-provided numpy/cython (needs --no-build-isolation).
#
# Usage:  build_telescope_wheel.sh <python-executable-of-target-env>
#   e.g.  build_telescope_wheel.sh $E/bin/python   (E = the te_erv conda env)
# Produces: pipelines/te_erv/vendor/telescope_ngs-1.0.3-*-arm64.whl
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
PY="${1:?usage: build_telescope_wheel.sh <python-of-target-env>}"
SRC_TARBALL="$HERE/telescope-1.0.3.tar.gz"
BUILD="$HERE/build"

[ -s "$SRC_TARBALL" ] || {
  echo "[info] fetching Telescope v1.0.3 source tarball (system curl)..."
  curl -fsSL -o "$SRC_TARBALL" \
    "https://github.com/mlbendall/telescope/archive/refs/tags/v1.0.3.tar.gz"
}

rm -rf "$BUILD"; mkdir -p "$BUILD"
tar xzf "$SRC_TARBALL" -C "$BUILD"
cd "$BUILD/telescope-1.0.3"

# (a) add the source dir to Cython's include_path + pin language_level
python3 - <<'PYEOF'
s = open("setup.py").read()
old = "    extensions = cythonize(extensions)"
new = ('    extensions = cythonize(extensions,\n'
       '        include_path=["telescope/utils"],\n'
       '        compiler_directives={"language_level": "2"})')
if old in s:
    open("setup.py","w").write(s.replace(old, new))
    print("[patch] cythonize include_path + language_level=2")
PYEOF

# (b) strip the self-cimport line that modern Cython rejects
python3 - <<'PYEOF'
p = "telescope/utils/calignment.pyx"
lines = open(p).read().splitlines(keepends=True)
out = [l for l in lines if l.strip() != "from calignment cimport AlignedPair"]
open(p,"w").write("".join(out))
print(f"[patch] removed self-cimport ({len(lines)} -> {len(out)} lines)")
PYEOF

# (c) replace deprecated numpy aliases removed in NumPy >=1.24 (np.int / np.float).
#     Telescope v1.0.3 predates the removal; our conda env has numpy 1.26.
#     Precision-preserving substitutions (np.int -> np.int64, np.float -> np.float64).
python3 - <<'PYEOF'
import re, glob
n = 0
for f in glob.glob("telescope/**/*.py", recursive=True):
    s = open(f).read()
    s2 = re.sub(r"np\.int\b(?!\d)", "np.int64", s)
    s2 = re.sub(r"np\.float\b(?!\d)", "np.float64", s2)
    s2 = re.sub(r"np\.bool\b(?!_|\d)", "np.bool_", s2)
    if s2 != s:
        open(f,"w").write(s2); n += 1
print(f"[patch] numpy-alias fix applied to {n} file(s)")
PYEOF

rm -rf build telescope/utils/calignment.c
"$PY" -m pip wheel --no-build-isolation --no-deps -w "$HERE" .
echo "[done] wheel(s) in $HERE:"
ls -la "$HERE"/telescope_ngs-*.whl
