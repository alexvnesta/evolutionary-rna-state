#!/usr/bin/env python
"""Tests for registry_update — the concurrent-writer-safe feature registry."""
import os, json, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis.registry_update import register_features


def test_sequential_merge_unions_keys():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "reg.json")
        register_features(p, {"a": {"bucket": "baseline"}}, bucket_of={"a": "baseline"})
        register_features(p, {"b": {"bucket": "differentiated"}}, bucket_of={"b": "differentiated"})
        register_features(p, {"c": {"bucket": "differentiated"}}, bucket_of={"c": "differentiated"})
        reg = json.load(open(p))
        assert set(reg["features"]) == {"a", "b", "c"}, reg["features"]
        assert reg["_meta"]["buckets"]["differentiated"] == ["b", "c"]
        assert len(reg["_meta"]["reconciled_from"]) == 3
        print("PASS test_sequential_merge_unions_keys")


def test_same_key_latest_wins():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "reg.json")
        register_features(p, {"a": {"v": 1}})
        register_features(p, {"a": {"v": 2}})
        reg = json.load(open(p))
        assert reg["features"]["a"]["v"] == 2
        assert set(reg["features"]) == {"a"}
        print("PASS test_same_key_latest_wins")


def test_corrupt_file_recovered():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "reg.json")
        open(p, "w").write("{partial garbage")  # simulate a torn legacy write
        register_features(p, {"a": {"v": 1}})
        reg = json.load(open(p))
        assert reg["features"] == {"a": {"v": 1}}
        print("PASS test_corrupt_file_recovered")


def test_concurrent_writers_no_clobber():
    """8 parallel processes each register a distinct feature; all must survive."""
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "reg.json")
        procs = []
        for i in range(8):
            code = (f"import sys; sys.path.insert(0, {os.getcwd()!r});"
                    f"from analysis.registry_update import register_features;"
                    f"register_features({p!r}, {{'feat_{i}': {{'bucket': 'differentiated'}}}},"
                    f" bucket_of={{'feat_{i}': 'differentiated'}})")
            procs.append(subprocess.Popen([sys.executable, "-c", code]))
        for pr in procs:
            pr.wait()
        reg = json.load(open(p))
        assert set(reg["features"]) == {f"feat_{i}" for i in range(8)}, reg["features"]
        assert len(reg["_meta"]["reconciled_from"]) == 8
        print("PASS test_concurrent_writers_no_clobber")


if __name__ == "__main__":
    test_sequential_merge_unions_keys()
    test_same_key_latest_wins()
    test_corrupt_file_recovered()
    test_concurrent_writers_no_clobber()
    print("\nALL REGISTRY-UPDATE TESTS PASSED")
