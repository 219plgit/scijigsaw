#!/usr/bin/env python3
"""VAMP2-centred real-structure micro-benchmark (Level 1: chain-pair interface detection).

This evaluates the structure-extraction tier against curated chain-pair labels over a
cutoff sweep, and reports precision / recall / F1 per structure and in aggregate.

It does NOT fabricate results: it runs the real extractor on real deposited structures
that you download, and compares its calls to reference_labels.csv.

USAGE
    1. Download the mmCIF/PDB files into ./structures/ (see README.md), one per PDB id,
       named <PDB>.cif (or .pdf). e.g. structures/1SFC.cif
    2. python run_benchmark.py --structures ./structures --labels reference_labels.csv \
           --cutoffs 4 5 6 8 --min-residues 3 5 10 --out results.csv
    3. Inspect results.csv and benchmark_figure.pdf.

Level 2 (residue-level agreement vs PDBePISA) and the site-overlap experiment are
described in README.md; they require the PISA/RCSB reference residue sets, which are not
bundled here.
"""
import argparse
import csv
import itertools
import os
import subprocess
import sys
from collections import defaultdict


def load_labels(path):
    ref = defaultdict(dict)  # pdb -> {frozenset({a,b}): "positive"/"negative"}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            ref[row["pdb"]][frozenset((row["chain_a"], row["chain_b"]))] = row["label"]
    return ref


def run_extractor(struct_path, cutoff, min_res, workdir):
    """Call scijigsaw-extract and return the set of detected chain pairs (frozensets)."""
    out = os.path.join(workdir, "interactions.csv")
    cmd = [
        "scijigsaw-extract", struct_path,
        "--contact-cutoff", str(cutoff),
        "--min-residues", str(min_res),
        "--out", out,
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        sys.stderr.write(f"extractor failed for {struct_path} @ {cutoff}A/{min_res}: {e}\n")
        return None
    pairs = set()
    if os.path.exists(out):
        with open(out) as fh:
            for row in csv.DictReader(fh):
                pairs.add(frozenset((row["protein_a"], row["protein_b"])))
    return pairs


def prf(tp, fp, fn):
    p = tp / (tp + fp) if tp + fp else float("nan")
    r = tp / (tp + fn) if tp + fn else float("nan")
    f = 2 * p * r / (p + r) if p and r and (p + r) else float("nan")
    return p, r, f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--structures", default="./structures")
    ap.add_argument("--labels", default="reference_labels.csv")
    ap.add_argument("--cutoffs", nargs="+", type=float, default=[4, 5, 6, 8])
    ap.add_argument("--min-residues", nargs="+", type=int, default=[3, 5, 10])
    ap.add_argument("--out", default="results.csv")
    args = ap.parse_args()

    ref = load_labels(args.labels)
    workdir = os.path.join(os.path.dirname(args.out) or ".", "_tmp")
    os.makedirs(workdir, exist_ok=True)

    rows = []
    for cutoff, min_res in itertools.product(args.cutoffs, args.min_residues):
        agg_tp = agg_fp = agg_fn = 0
        for pdb, labels in ref.items():
            spath = None
            for ext in (".cif", ".mmcif", ".pdf", ".ent"):
                cand = os.path.join(args.structures, pdb + ext)
                if os.path.exists(cand):
                    spath = cand
                    break
            if spath is None:
                sys.stderr.write(f"missing structure for {pdb}; skipping\n")
                continue
            detected = run_extractor(spath, cutoff, min_res, workdir)
            if detected is None:
                continue
            positives = {pair for pair, lab in labels.items() if lab == "positive"}
            negatives = {pair for pair, lab in labels.items() if lab == "negative"}
            tp = len(positives & detected)
            fn = len(positives - detected)
            fp = len(negatives & detected)  # false positive = detected a labelled negative
            agg_tp += tp; agg_fp += fp; agg_fn += fn
            p, r, f = prf(tp, fp, fn)
            rows.append(dict(cutoff=cutoff, min_residues=min_res, pdb=pdb,
                             tp=tp, fp=fp, fn=fn, precision=p, recall=r, f1=f))
        p, r, f = prf(agg_tp, agg_fp, agg_fn)
        rows.append(dict(cutoff=cutoff, min_residues=min_res, pdb="AGGREGATE",
                         tp=agg_tp, fp=agg_fp, fn=agg_fn, precision=p, recall=r, f1=f))

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["cutoff", "min_residues", "pdb",
                                           "tp", "fp", "fn", "precision", "recall", "f1"])
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {args.out} ({len(rows)} rows)")
    print("Aggregate rows:")
    for row in rows:
        if row["pdb"] == "AGGREGATE":
            print("  {cutoff}A/{min_residues}res  P={precision:.3f} R={recall:.3f} "
                  "F1={f1:.3f}".format(**row))


if __name__ == "__main__":
    main()
