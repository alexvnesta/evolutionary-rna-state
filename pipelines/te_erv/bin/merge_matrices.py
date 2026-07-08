#!/usr/bin/env python3
"""Merge per-sample TE/ERV count files into a single count matrix.

Usage:
    merge_matrices.py telescope      sampleA-telescope_report.tsv sampleB-... > locus_matrix.tsv
    merge_matrices.py tetranscripts  sampleA.cntTable sampleB.cntTable ...      > family_matrix.tsv

Telescope report format: commented header line, then a header row with columns
    transcript  transcript_length  final_count  final_conf  final_prop  ...
We take the 'final_count' column keyed by 'transcript' (locus).

TEcount .cntTable format: two columns  'gene/TE\\t<sample>'  (feature<TAB>count).
Sample name is taken from the file's own header column, falling back to the
filename stem.
"""
import sys, os, csv


def sample_stem(path, suffix):
    b = os.path.basename(path)
    for s in suffix:
        if b.endswith(s):
            return b[: -len(s)]
    return os.path.splitext(b)[0]


def parse_telescope(path):
    """Return (sample, {locus: count})."""
    sample = sample_stem(path, ["-telescope_report.tsv", ".tsv"])
    counts = {}
    with open(path) as fh:
        rows = [l.rstrip("\n") for l in fh if not l.startswith("##")]
    if not rows:
        return sample, counts
    header = rows[0].split("\t")
    try:
        i_id = header.index("transcript")
    except ValueError:
        i_id = 0
    # prefer final_count, else final_conf, else 2nd col
    i_cnt = header.index("final_count") if "final_count" in header else 2
    for line in rows[1:]:
        f = line.split("\t")
        if len(f) <= max(i_id, i_cnt):
            continue
        locus = f[i_id]
        if locus == "__no_feature":
            continue
        try:
            counts[locus] = int(round(float(f[i_cnt])))
        except ValueError:
            pass
    return sample, counts


def parse_tecount(path):
    """Return (sample, {feature: count}) from a TEcount .cntTable."""
    with open(path) as fh:
        reader = list(csv.reader(fh, delimiter="\t"))
    if not reader:
        return sample_stem(path, [".cntTable"]), {}
    header = reader[0]
    sample = header[1] if len(header) > 1 and header[1] else sample_stem(path, [".cntTable"])
    counts = {}
    for row in reader[1:]:
        if len(row) < 2:
            continue
        try:
            counts[row[0]] = int(round(float(row[1])))
        except ValueError:
            pass
    return sample, counts


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    kind, files = sys.argv[1], sys.argv[2:]
    parse = parse_telescope if kind == "telescope" else parse_tecount
    feature_label = "locus" if kind == "telescope" else "feature"

    per_sample = {}
    all_features = []
    seen = set()
    for f in files:
        s, c = parse(f)
        per_sample[s] = c
        for k in c:
            if k not in seen:
                seen.add(k)
                all_features.append(k)

    samples = list(per_sample.keys())
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow([feature_label] + samples)
    for feat in all_features:
        w.writerow([feat] + [per_sample[s].get(feat, 0) for s in samples])


if __name__ == "__main__":
    main()
