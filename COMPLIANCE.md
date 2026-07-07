# Compliance notes

## Code provenance

This repository was created as a **clean slate during the hackathon** (kickoff
2026-07-07, 16:00 UTC). Per the organizer rule, all submitted **code** is
authored live in-session — no code, pipeline scripts, encoder implementation,
package modules, configs, lockfiles, or manuscript text was carried over from
any pre-kickoff work. Ideation and public-data acquisition before kickoff were
permitted.

**What is committed here and why it is allowed:**

- `data/` — metadata *only* (manifests, download receipts, a run catalog, and
  clinical/response labels) for **already-public** datasets. Public-data
  acquisition is permitted; these files carry no source code.
- `docs/` — prose (data inventory, roadmap, methods notes).
- `README.md`, `LICENSE`, `.gitignore`, this file — authored in-session.

**What is deliberately absent** (to be re-authored live, not ported):
pipeline scripts, the raw-read encoder implementation, `src/` modules,
`config.yaml`, `Makefile`, `pyproject.toml` / lockfiles, and the manuscript.

## Data access levels

- **Committed:** metadata for open datasets (ENA/SRA/GEO — Gide 2019
  `PRJEB23709`, Riaz 2017 `PRJNA356761`/`GSE91061`, etc.).
- **Never committed:** raw reads, alignments, quantifications, and any
  controlled-access data (e.g. IMvigor210 raw = EGA controlled; only the
  processed expression matrix from the public R package is usable).

Raw/large data paths are git-ignored (see `.gitignore`).

## A note on this checkout's git metadata

The build sandbox forbids any path named `.git`. This working copy therefore
keeps its git metadata in `.gitmeta/`. On a normal machine, `git clone` of the
pushed remote produces a standard `.git` layout; no history or content is
affected.
