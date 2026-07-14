#!/usr/bin/env python3
"""Regenerate every NUMBER reported in the paper.  ~25 s.

    python scripts/reproduce_numbers.py

Exits non-zero if the code and the manuscript disagree.
"""
import argparse, json, os, sys, time
import numpy as np
from scijigsaw import INFLAMMASOME, VAMP2, count_30S
from scijigsaw.benchmark import analyse, run

CHECKS = []


def check(label, got, expect, tol=0.0):
    ok = abs(got - expect) <= tol if isinstance(expect, (int, float)) else got == expect
    CHECKS.append((label, got, expect, ok))
    print(f"  [{'ok ' if ok else 'FAIL'}] {label:<42} {got}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="outputs")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    print("RESULTS 1 -- three encoded assemblies")
    v, i, r = VAMP2.summary(), INFLAMMASOME.summary(), count_30S()
    check("VAMP2 permitted / 5040", v["permitted"], 252)
    check("VAMP2 depth", v["depth"], 3)
    check("inflammasome permitted / 3,628,800", i["permitted"], 2)
    check("inflammasome depth", i["depth"], 9)
    check("30S units", r["n"], 20)
    check("30S reduction (~6x)", round(r["reduction"], 1), 6.0, 0.6)

    print("\nRESULTS 2 -- random-poset benchmark (900 posets, seed 0)")
    t0 = time.perf_counter()
    A = run(n_posets=900, seed=0)
    B = analyse(A)
    print(f"  ({time.perf_counter()-t0:.0f} s)")
    check("Pearson r, subunit count", round(B["pearson"]["n"], 2), 0.49, 0.02)
    check("Pearson r, depth", round(B["pearson"]["depth"], 2), 0.93, 0.02)
    check("R2 (n, depth, width, density)", round(B["r2_full"], 2), 0.97, 0.02)
    check("depth adds to R2", round(B["depth_adds"], 3), 0.016, 0.01)
    check("VIF depth", round(B["vif"]["depth"], 1), 6.9, 1.0)

    n_, d_, w_, dens_, br_, y = A.T
    print("\n  depth at fixed n (mean log10 reduction)")
    for n0 in (10, 12, 14):
        m = n_ == n0
        cells = []
        for lo, hi in [(1, 3), (4, 6), (7, 20)]:
            k = m & (d_ >= lo) & (d_ <= hi)
            cells.append(f"{y[k].mean():5.2f}" if k.sum() >= 5 else "  -  ")
        print(f"    n={n0}: 1-3 {cells[0]}  4-6 {cells[1]}  7+ {cells[2]}")

    json.dump({"vamp2": v, "inflammasome": i, "ribosome_30S": r, "benchmark": B},
              open(os.path.join(a.out, "results.json"), "w"), indent=2, default=float)

    bad = [c for c in CHECKS if not c[3]]
    print()
    if bad:
        print(f"{len(bad)} CHECK(S) FAILED - manuscript disagrees with the code")
        for lab, got, exp, _ in bad:
            print(f"   {lab}: got {got}, manuscript says {exp}")
        sys.exit(1)
    print(f"All {len(CHECKS)} reported values reproduced -> {a.out}/results.json")


if __name__ == "__main__":
    main()
