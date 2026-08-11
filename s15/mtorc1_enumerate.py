#!/usr/bin/env python3
"""
Step 3 of the mTORC1 / rapamycin example: the actual Scientific Jigsaw run.

Builds two DECLARED model variants over the mTORC1 board and enumerates each
across the (cell line, drug) contexts produced by tahoe_availability.py.

  M0      untreated. Contacts extracted from apo mTORC1 coordinates.

  M_Rapa  rapamycin-treated. Two declared changes, each cited separately:
            (a) CONTACT   FKBP12-mTOR(FRB), ligand-conditioned on rapamycin.
                          Evidence: 1FAP / 5FLC.
            (b) EXCLUSION mTOR-RAPTOR. Evidence: the biochemical and cryo-EM
                          literature reporting that FKBP12-rapamycin disrupts
                          the mTOR-RAPTOR association and the mTORC1 dimer.
                          This is NOT extracted from a structure and must not
                          be cited to one.

The perturbation therefore indexes a different admissible space rather than
reweighting a fixed one. Both changes are declared and evidence-linked; the
Tahoe condition label selects which variant applies, and supplies nothing else.

Usage
-----
    python mtorc1_enumerate.py --contexts contexts.csv --out results.csv
"""
import argparse
import csv
from itertools import combinations

SUBUNITS = ["mTOR", "RAPTOR", "mLST8", "PRAS40", "DEPTOR"]

# PLACEHOLDER contact graph -- REPLACE with the output of mtorc1_contacts.py.
# Retained only so the pipeline is runnable end to end before extraction.
M0_CONTACTS = [
    ("mTOR", "RAPTOR"),
    ("mTOR", "mLST8"),
    ("mTOR", "DEPTOR"),
    ("RAPTOR", "PRAS40"),
]
RAPA_CONTACT = ("FKBP12", "mTOR")
RAPA_EXCLUSION = ("mTOR", "RAPTOR")

OBSERVED = [("mTOR", "RAPTOR", "mLST8")]   # the mTORC1 heterotrimer


def _conn(A, B, edges):
    return any(frozenset((u, v)) in edges for u in A for v in B)


def _subsets(S):
    S = sorted(S)
    for r in range(1, len(S)):
        for c in combinations(S, r):
            yield frozenset(c)


def merger_trees(S, edges, _memo=None):
    """Exact count of binary merger trees over S under the contact graph."""
    if _memo is None:
        _memo = {}
    key = (S, frozenset(edges))
    if key in _memo:
        return _memo[key]
    if len(S) <= 1:
        return 1
    total, seen = 0, set()
    for A in _subsets(S):
        B = S - A
        k = frozenset((A, B))
        if k in seen:
            continue
        seen.add(k)
        if not _conn(A, B, edges):
            continue
        total += merger_trees(A, edges, _memo) * merger_trees(B, edges, _memo)
    _memo[key] = total
    return total


def fragments(S, edges):
    out, rem = [], set(S)
    while rem:
        seed = rem.pop()
        comp, stack = {seed}, [seed]
        while stack:
            u = stack.pop()
            for v in list(rem):
                if frozenset((u, v)) in edges:
                    rem.discard(v)
                    comp.add(v)
                    stack.append(v)
        out.append(frozenset(comp))
    return out


def build(model, present):
    comps = set(present) & (set(SUBUNITS) | {"FKBP12"})
    edges = {frozenset(e) for e in M0_CONTACTS}
    if model == "M_Rapa":
        edges.discard(frozenset(RAPA_EXCLUSION))
        if "FKBP12" in comps:
            edges.add(frozenset(RAPA_CONTACT))
    else:
        comps.discard("FKBP12")
    comps = frozenset(comps)
    edges = {e for e in edges if e <= comps}
    return comps, edges


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contexts", required=True)
    ap.add_argument("--out", default="results.csv")
    args = ap.parse_args()

    rows = []
    with open(args.contexts) as fh:
        for ctx in csv.DictReader(fh):
            present = ctx["components_present"].split("|")
            comps, edges = build(ctx["model"], present)
            trees = merger_trees(comps, edges)
            frags = fragments(comps, edges)
            for obs in OBSERVED:
                o = frozenset(obs)
                verdict = ("impossible" if not (o <= comps)
                           or not any(o <= f for f in frags) else "formable")
                rows.append({
                    "cell_line": ctx["cell_line"],
                    "drug": ctx["drug"],
                    "threshold_cpm": ctx["threshold_cpm"],
                    "model": ctx["model"],
                    "holo_merger_trees": trees,
                    "n_fragments": len(frags),
                    "fragments": " ; ".join(
                        "+".join(sorted(f)) for f in
                        sorted(frags, key=lambda x: -len(x))),
                    "observed": "+".join(obs),
                    "verdict": verdict,
                })

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    imp = sum(1 for r in rows if r["verdict"] == "impossible")
    print(f"wrote {len(rows)} rows -> {args.out}")
    print(f"impossibility calls: {imp}/{len(rows)}")
    print("\nReport in the manuscript:\n"
          "  (1) |Omega_0| vs |Omega_Rapa| for the reference context;\n"
          "  (2) the fragmentation induced by the declared exclusion;\n"
          "  (3) contexts where a subunit is unavailable and the heterotrimer\n"
          "      is IMPOSSIBLE -- the robust direction of the evidence.")


if __name__ == "__main__":
    main()
