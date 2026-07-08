#!/usr/bin/env python
"""
registry_update.py — safe concurrent updates to feature_registry.json.

The 12 feature modules were built by parallel sub-agents that each wrote
feature_registry.json by overwrite, so the last writer clobbered the others'
entries (the "concurrent-writer collision" that forked the registry). This
module replaces overwrite-semantics with **atomic read-merge-write under an
advisory file lock**, so any number of writers can register their features
without losing each other's:

- merge unions feature entries by key (later registration wins for the SAME key,
  never drops OTHER keys),
- `_meta` is merged shallowly with a running `reconciled_from` provenance list,
- the write is atomic (temp file + os.replace) so a reader never sees a partial
  file, and serialized by an flock so two writers cannot interleave.

Usage (each module, instead of json.dump over the file):
    from analysis.registry_update import register_features
    register_features("results/features/feature_registry.json",
                      {"te_antigen_burden": {...}, ...},
                      meta_note="te module", source_id="<artifact vid>")
"""
import os, json, tempfile, time

try:
    import fcntl  # POSIX advisory locks
    _HAVE_FCNTL = True
except ImportError:  # pragma: no cover - non-POSIX
    _HAVE_FCNTL = False


def _empty_registry():
    return {"_meta": {"schema": "feature_registry/v2", "reconciled_from": [],
                      "buckets": {}, "note": ""}, "features": {}}


def _load(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return _empty_registry()
    try:
        with open(path) as fh:
            d = json.load(fh)
    except (json.JSONDecodeError, ValueError):
        # corrupt/partial (a legacy non-atomic write) — start clean, don't crash
        return _empty_registry()
    d.setdefault("_meta", {}).setdefault("reconciled_from", [])
    d.setdefault("features", {})
    return d


def _merge(base, new_features, meta_note=None, source_id=None, bucket_of=None):
    for k, v in new_features.items():
        base["features"][k] = v            # union by key; same key -> latest wins
        if bucket_of:
            b = bucket_of.get(k)
            if b:
                base["_meta"].setdefault("buckets", {}).setdefault(b, [])
                if k not in base["_meta"]["buckets"][b]:
                    base["_meta"]["buckets"][b].append(k)
    prov = {"at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "n_features": len(new_features), "keys": sorted(new_features)}
    if meta_note:
        prov["note"] = meta_note
    if source_id:
        prov["source_id"] = source_id
    base["_meta"]["reconciled_from"].append(prov)
    return base


def _atomic_write(path, data):
    d = os.path.dirname(os.path.abspath(path))
    os.makedirs(d, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=d, prefix=".reg_", suffix=".json")
    try:
        with os.fdopen(fd, "w") as fh:
            json.dump(data, fh, indent=1, sort_keys=True)
        os.replace(tmp, path)              # atomic on POSIX
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def register_features(path, new_features, meta_note=None, source_id=None,
                      bucket_of=None, timeout=30.0):
    """Atomically merge ``new_features`` into the registry at ``path`` under an
    advisory lock. Returns the merged registry dict. Never drops entries written
    by a concurrent writer."""
    lock_path = path + ".lock"
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    deadline = time.time() + timeout
    lockf = open(lock_path, "w")
    try:
        if _HAVE_FCNTL:
            while True:
                try:
                    fcntl.flock(lockf, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except BlockingIOError:
                    if time.time() > deadline:
                        raise TimeoutError(f"registry lock busy > {timeout}s: {path}")
                    time.sleep(0.05)
        reg = _load(path)
        reg = _merge(reg, new_features, meta_note=meta_note, source_id=source_id,
                     bucket_of=bucket_of)
        _atomic_write(path, reg)
        return reg
    finally:
        if _HAVE_FCNTL:
            try:
                fcntl.flock(lockf, fcntl.LOCK_UN)
            except OSError:
                pass
        lockf.close()
