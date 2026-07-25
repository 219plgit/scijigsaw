#!/usr/bin/env python
"""Evaluate the structure-extraction tier against curated chain-pair labels.

This does NOT fabricate results. It runs the same interface-detection code used
by the extractor (scijigsaw.extract.interface_residues) on real deposited
structures that you download, and compares its calls to reference_labels.csv
over a sweep of contact cutoffs and minimum interface sizes.

The labels name PROTEINS (VAMP2, Syntaxin-1A, ...), while structures contain
CHAINS (A, B, ...). The mapping between them is therefore required. Run

    python run_benchmark.py --list-chains

first: it prints every chain in every structure with the molecule description
taken from the file header. Save a mapping as CSV with columns

    pdb,chain,protein

and pass it with --chain-map. A best-effort automatic match is attempted when no
map is supplied, and every match it makes is reported, but headers frequently use
names that differ from the labels (for example SYNAPTOBREVIN-2 for VAMP2), so the
explicit map is the reliable path.

Usage:
    python run_benchmark.py --list-chains
    python run_benchmark.py --chain-map chain_map.csv
"""
from __future__ import annotations

import argparse
import csv
import itertools
import os
import re
import sys
from collections import defaultdict

STRUCT_EXT = (".pdb", ".ent", ".cif", ".mmcif")


def load_labels(path):
    ref = defaultdict(dict)          # pdb -> {frozenset({a,b}): "positive"/"negative"}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            ref[row["pdb"]][frozenset((row["chain_a"], row["chain_b"]))] = row["label"]
    return ref


def find_structure(directory, pdb):
    for ext in STRUCT_EXT:
        for name in (pdb, pdb.lower(), pdb.upper()):
            cand = os.path.join(directory, name + ext)
            if os.path.exists(cand):
                return cand
    return None


def load_model(path):
    from Bio.PDB import PDBParser, MMCIFParser
    parser = (MMCIFParser(QUIET=True) if path.endswith((".cif", ".mmcif"))
              else PDBParser(QUIET=True))
    structure = parser.get_structure("x", path)
    return next(iter(structure)), structure


def chain_descriptions(path, structure):
    """chain id -> molecule description from the file header (may be empty)."""
    desc = {}
    if path.endswith((".cif", ".mmcif")):
        try:
            from Bio.PDB.MMCIF2Dict import MMCIF2Dict
            d = MMCIF2Dict(path)
            names = d.get("_entity.pdbx_description", [])
            strands = d.get("_entity_poly.pdbx_strand_id", [])
            if isinstance(names, str):
                names = [names]
            if isinstance(strands, str):
                strands = [strands]
            for name, strand in zip(names, strands):
                for ch in str(strand).split(","):
                    desc[ch.strip()] = str(name)
        except Exception:
            pass
    else:
        for mol in (structure.header.get("compound") or {}).values():
            name = mol.get("molecule", "")
            for ch in str(mol.get("chain", "")).split(","):
                if ch.strip():
                    desc[ch.strip().upper()] = name
    return desc


def _norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def auto_map(desc, proteins):
    """Best-effort protein -> [chains], by substring match on normalised names."""
    out = defaultdict(list)
    for ch, d in desc.items():
        nd = _norm(d)
        for p in proteins:
            np_ = _norm(p)
            if np_ and nd and (np_ in nd or nd in np_):
                out[p].append(ch)
    return out


def load_chain_map(path):
    m = defaultdict(lambda: defaultdict(list))   # pdb -> protein -> [chains]
    with open(path) as fh:
        for row in csv.DictReader(fh):
            m[row["pdb"]][row["protein"]].append(row["chain"].strip())
    return m


def prf(tp, fp, fn):
    p = tp / (tp + fp) if tp + fp else float("nan")
    r = tp / (tp + fn) if tp + fn else float("nan")
    f = 2 * p * r / (p + r) if p and r and (p + r) else float("nan")
    return p, r, f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--structures", default="./structures")
    ap.add_argument("--labels", default="reference_labels.csv")
    ap.add_argument("--chain-map", default=None,
                    help="CSV with columns pdb,chain,protein")
    ap.add_argument("--list-chains", action="store_true",
                    help="print chains and header descriptions, then exit")
    ap.add_argument("--cutoffs", nargs="+", type=float, default=[4, 5, 6, 8])
    ap.add_argument("--min-residues", nargs="+", type=int, default=[3, 5, 10])
    ap.add_argument("--out", default="results.csv")
    a = ap.parse_args()

    from scijigsaw.extract import interface_residues

    ref = load_labels(a.labels)

    # ---- inventory / mapping -------------------------------------------------
    if a.list_chains:
        for pdb in ref:
            path = find_structure(a.structures, pdb)
            if path is None:
                print(f"{pdb}: NO STRUCTURE FILE in {a.structures}")
                continue
            model, structure = load_model(path)
            desc = chain_descriptions(path, structure)
            print(f"\n{pdb}  ({os.path.basename(path)})")
            for ch in model:
                n = sum(1 for _ in ch.get_residues())
                print(f"   chain {ch.id:<3} {n:>5} residues   {desc.get(ch.id.upper(), '')}")
            print(f"   labels mention: "
                  f"{sorted({p for pair in ref[pdb] for p in pair})}")
        return

    cmap = load_chain_map(a.chain_map) if a.chain_map else None

    rows = []
    for cutoff, min_res in itertools.product(a.cutoffs, a.min_residues):
        agg_tp = agg_fp = agg_fn = 0
        for pdb, labels in ref.items():
            path = find_structure(a.structures, pdb)
            if path is None:
                sys.stderr.write(f"missing structure for {pdb}; skipping\n")
                continue
            model, structure = load_model(path)
            proteins = sorted({p for pair in labels for p in pair})

            if cmap and pdb in cmap:
                mapping = {p: list(cmap[pdb].get(p, [])) for p in proteins}
            else:
                desc = chain_descriptions(path, structure)
                mapping = auto_map(desc, proteins)
                matched = {p: v for p, v in mapping.items() if v}
                sys.stderr.write(f"{pdb}: auto-mapped {matched or 'NOTHING'}\n")

            unmapped = [p for p in proteins if not mapping.get(p)]
            if unmapped:
                sys.stderr.write(
                    f"{pdb}: no chains for {unmapped}; skipping. "
                    f"Run --list-chains and supply --chain-map.\n")
                continue

            detected = set()
            for pa, pb in itertools.combinations(proteins, 2):
                best = 0
                for ca, cb in itertools.product(mapping[pa], mapping[pb]):
                    if ca not in model or cb not in model:
                        continue
                    res, _ = interface_residues(model[ca], model[cb], cutoff)
                    best = max(best, len(res))
                if best >= min_res:
                    detected.add(frozenset((pa, pb)))

            positives = {pair for pair, lab in labels.items() if lab == "positive"}
            negatives = {pair for pair, lab in labels.items() if lab == "negative"}
            tp = len(positives & detected)
            fn = len(positives - detected)
            fp = len(negatives & detected)
            agg_tp += tp; agg_fp += fp; agg_fn += fn
            rows.append(dict(cutoff=cutoff, min_residues=min_res, pdb=pdb,
                             tp=tp, fp=fp, fn=fn))

        p, r, f = prf(agg_tp, agg_fp, agg_fn)
        rows.append(dict(cutoff=cutoff, min_residues=min_res, pdb="AGGREGATE",
                         tp=agg_tp, fp=agg_fp, fn=agg_fn,
                         precision=round(p, 3) if p == p else "",
                         recall=round(r, 3) if r == r else "",
                         f1=round(f, 3) if f == f else ""))

    fields = ["cutoff", "min_residues", "pdb", "tp", "fp", "fn",
              "precision", "recall", "f1"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fields})
    print(f"wrote {a.out} ({len(rows)} rows)")
    print("Aggregate rows:")
    for row in rows:
        if row["pdb"] == "AGGREGATE":
            print(f"  {row['cutoff']:g}A/{row['min_residues']}res  "
                  f"TP={row['tp']} FP={row['fp']} FN={row['fn']}  "
                  f"P={row['precision']} R={row['recall']} F1={row['f1']}")


if __name__ == "__main__":
    main()
