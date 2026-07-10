# Native arm64 STAR build — root cause of the "0 reads ingested" bug and the fix

**Date:** 2026-07-10
**Machine:** Apple Silicon (arm64), macOS, Apple clang 21.0.0, GNU Make 3.81
**STAR version built:** 2.7.11b (source tarball from github.com/alexdobin/STAR tag `2.7.11b`)

## Verdict (headline)

A native arm64 recompile **alone does NOT fix** the 0-reads bug — a fresh, unmodified
2.7.11b built with Apple clang *still* reports `Number of input reads | 0`, identical to
the bioconda osx-arm64 binary and the osx-64 build under Rosetta.

The bug is **not** a conda packaging problem, an arm64 SIMD/codegen problem, or a
char-signedness problem. It is a **libc++ vs libstdc++ standard-library difference** that
affects *every* STAR binary built against libc++ (which is the only practical C++ stdlib
on macOS). A **one-line source patch** fixes it completely: after the patch, STAR ingests
all **1,665,315** test reads, emits a full BAM, `SJ.out.tab`, and a non-empty
`Chimeric.out.junction` (the input STAR-Fusion consumes).

## Root cause

STAR distributes reads to worker threads through a fixed char buffer wired to an
`istringstream` via `pubsetbuf`:

```cpp
// ReadAlignChunk.cpp (constructor)
readInStream[ii]->rdbuf()->pubsetbuf(chunkIn[ii], P.chunkInSizeBytesArray);
```

`std::basic_stringbuf::setbuf` (what `pubsetbuf` calls) is a **documented no-op in
libc++** (Apple/LLVM's C++ standard library), whereas **libstdc++ (GNU/Linux) honors it**.
Under libc++ the `istringstream` therefore stays permanently empty: the `ifstream` reads
the FASTQ into `chunkIn` correctly, but that data never reaches the read parser, so the
first `peek()` returns EOF (-1) and STAR reports 0 reads while exiting 0 "successfully".

This was confirmed with a 6-line standalone program: `istringstream::pubsetbuf` returns
`peek() == -1` under libc++, vs the real first byte under libstdc++. HISAT2/bowtie2 don't
use this buffer trick, which is why they ingest the same reads fine.

`Log.out` fingerprint of the bug: `Thread #N end of input stream, nextChar=-1` before a
single read is processed.

## The fix (one line, in `source/ReadAlignChunk_mapChunk.cpp`)

In `ReadAlignChunk::mapChunk()`, before the streams are cleared/rewound, load the chunk
buffer content into the istringstream with `.str()` (which libc++ **does** honor) instead
of relying on the no-op `pubsetbuf`:

```cpp
for (uint ii=0;ii<P.readNends;ii++) {//clear eof and rewind the input streams
    // ARM64/macOS libc++ fix: basic_stringbuf::setbuf (pubsetbuf) is a no-op in libc++,
    // so the pubsetbuf() wiring in the constructor leaves the istringstream empty
    // (=> "0 input reads"). Load the chunk buffer content directly via str() instead.
    readInStream[ii]->str(string(chunkIn[ii], chunkInSizeBytesTotal[ii]));   // <-- added
    RA->readInStream[ii]->clear();
    RA->readInStream[ii]->seekg(0,ios::beg);
};
```

## Build recipe (arm64 macOS, Apple clang)

Prereqs: Xcode CLT (clang++), GNU make, Homebrew `libomp` (Apple clang ships no OpenMP
runtime). Build in /tmp.

```bash
# 1. source (no git clone — sandbox blocks .git; use release tarball)
curl -L -o star.tar.gz https://github.com/alexdobin/STAR/archive/refs/tags/2.7.11b.tar.gz
tar --exclude='.git*' -xzf star.tar.gz
cd STAR-2.7.11b/source

# 2. Makefile patch: Apple clang has no builtin -fopenmp; use -Xpreprocessor form + libomp headers
sed -i.bak 's#-std=c++11 -fopenmp #-std=c++11 -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include #' Makefile

# 3. apply the one-line libc++ read-ingest fix to ReadAlignChunk_mapChunk.cpp (see above)

# 4. build htslib (bundled)
make -C htslib lib-static CC=clang

# 5. build STAR — key flags:
LIBOMP=/opt/homebrew/opt/libomp
make STAR \
  CXX=clang++ \
  CXXFLAGS_SIMD="" \
  CXXFLAGSextra="-DCOMPILE_FOR_MAC" \
  "LDFLAGS_shared=-pthread -L$LIBOMP/lib -lomp -lz htslib/libhts.a" \
  -j8
```

### Why each flag matters
- **`CXXFLAGS_SIMD=""`** — the default is `-mavx2` (x86-only; Apple clang rejects it). STAR's
  `opal` aligner ships the amalgamated **SIMDe** header (`simde_avx2.h`), which transparently
  translates the x86 AVX2/SSE intrinsics to ARM NEON, so no `-m` SIMD flag is needed.
- **`-Xpreprocessor -fopenmp -I<libomp>/include`** (Makefile patch) + **`-L<libomp>/lib -lomp`**
  — Apple clang has no builtin `-fopenmp` and bundles no OpenMP runtime; use Homebrew libomp.
- **`CXXFLAGSextra="-DCOMPILE_FOR_MAC"`** — guards Linux-only `SHM_NORESERVE` (shared-memory
  flag) so it compiles on macOS.
- Do **not** override `GIT_BRANCH_COMMIT_DIFF` or `COMPTIMEPLACE` on the command line — the
  Makefile turns them into `-D` defines internally; overriding them to empty breaks the build.
- `-fsigned-char` was tested and is **not** needed (bug is libc++, not char signedness).

Resulting binary: `file STAR` → `Mach-O 64-bit executable arm64`.

## Portability

The binary links Homebrew's `libomp.dylib` by absolute path. For a self-contained artifact,
bundle it and rewrite the load path to an rpath, then re-sign (arm64 refuses a binary whose
signature `install_name_tool` invalidated — it gets `Killed: 9`):

```bash
install_name_tool -change /opt/homebrew/opt/libomp/lib/libomp.dylib '@rpath/libomp.dylib' STAR
install_name_tool -add_rpath '@loader_path/../lib' STAR
install_name_tool -id '@rpath/libomp.dylib' libomp.dylib
codesign --force -s - libomp.dylib
codesign --force -s - STAR
```
Layout: `bin/STAR-arm64-native` + `lib/libomp.dylib`. Runs with no `DYLD_*` env.

## STAR-Fusion setup

- STAR-Fusion **v1.15.1** FULL tarball (bundled deps). The `STAR-Fusion` executable is a Perl
  script; it runs on system Perl 5.34 and reports `STAR-Fusion version: 1.15.0`, and finds the
  native STAR on `PATH`.
- Needed non-core Perl modules: `DB_File`, `URI::Escape`, `JSON::XS`, `PerlIO::gzip` (all
  present in system Perl) + **`Set::IntervalTree`** (used by core `STAR-Fusion.map_chimeric_reads_to_genes`
  and `.handle_multimapping_reads`).
- `Set::IntervalTree` is bundled in `plugins/Set-IntervalTree-0.01.tar.gz` but its XS won't
  compile as-is on modern macOS. Two patches to `IntervalTree.xs`:
  1. `#include <tr1/memory>` → `#include <memory>`, and `std::tr1::shared_ptr` → `std::shared_ptr`
     (tr1 was removed from the C++ stdlib years ago).
  2. Add `#undef do_open` / `#undef do_close` (+ `seed bind Copy Move Zero New`) after the Perl
     headers and before the C++ headers — Perl's function-like macros collide with libc++
     `<locale>`/`std::messages` methods.
  Then `perl Makefile.PL && make`, and place `Set/` + `auto/Set/` under STAR-Fusion's `PerlLib/`.

## What a full nf-core/rnafusion (STAR-Fusion) run still needs

- **CTAT genome resource library** (`GRCh38_gencode_v??_CTAT_lib_*.plug-n-play.tar.gz`, ~30–40 GB
  download, expands larger). This bundles the STAR index, ref annotation, fusion-annotation
  databases (FusionAnnotator), PFAM/coding-effect data, etc. It is the large piece deliberately
  **not** built here.
- A full-genome STAR index (part of the CTAT lib) — needs ~30 GB RAM to build; use the prebuilt
  CTAT plug-n-play lib rather than building it.
- For nf-core/rnafusion specifically: point the pipeline at this STAR binary (arm64-native,
  patched) and a CTAT lib matching the GENCODE annotation.
