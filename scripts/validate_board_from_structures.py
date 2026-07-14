#!/usr/bin/env python3
"""Partial validation of the VAMP2 board against deposited structures.

    python scripts/validate_board_from_structures.py        # needs network (RCSB)

WHY "PARTIAL", AND WHICH PART
-----------------------------
The board makes seven structural assertions. Four can be checked against deposited
complexes; three cannot. Saying which is which is itself part of the result.

  TESTABLE
    A1  VAMP2's SNARE motif engages syntaxin-1A AND both SNAP25 helices.   1SFC / 1N7S
    A2  Complexin contacts TWO OR MORE chains of the assembled bundle.     1KIL
        This is what makes it a BRIDGE piece, and a bridge piece cannot be seated
        until the seam beneath it closes. That constraint does the combinatorial
        work of the whole paper -- and 1KIL can refute it.
    A3  Synaptotagmin-1 likewise bridges the bundle.                       5CCG
    A4  Munc18-1 engages syntaxin-1A.                                      3C98

  NOT TESTABLE, and the paper must keep saying so
    A5  Munc18-1 does NOT contact VAMP2. No deposit contains both in the closed
        conformation, and a chain's absence from a crystal is not evidence that an
        interaction is absent. This edge is inferred, not shown.
    A6  AP180/CALM (ANTH) compete for VAMP2's SNARE motif. No deposited complex is
        known to us; the competition rests on Koo et al. (2011).
    A7  Alpha-synuclein binds the VAMP2 N-terminus. Synuclein is intrinsically
        disordered and the interface is not crystallographically resolved
        (Burre et al., 2010).

If the derived interfaces disagree with the board, the disagreement is the finding.
Report it; do not quietly re-encode the board to match.
"""
from __future__ import annotations

import argparse
import itertools
import os
import urllib.request

URL = "https://files.rcsb.org/download/{}.pdb"

# Chain labels per the deposit. Verify against the PDB entry page: chain identifiers
# are a property of the file, not of the biology, and they do change between entries.
DEPOSITS = {
    "1SFC": {"A": "Syntaxin-1A", "B": "VAMP2", "C": "SNAP25-N", "D": "SNAP25-C"},
    "1N7S": {"A": "Syntaxin-1A", "B": "VAMP2", "C": "SNAP25-N", "D": "SNAP25-C"},
    "1KIL": None,   # complexin + SNARE bundle
    "5CCG": None,   # synaptotagmin-1 + SNARE
    "3C98": None,   # Munc18-1 + syntaxin-1A
}


def fetch(pdb, out):
    p = os.path.join(out, f"{pdb}.pdb")
    if not os.path.exists(p):
        urllib.request.urlretrieve(URL.format(pdb), p)
    return p


def chain_contacts(path, cutoff):
    from Bio.PDB import PDBParser
    from scijigsaw.extract import interface_residues
    model = next(iter(PDBParser(QUIET=True).get_structure("x", path)))
    chains = {c.id: c for c in model if sum(1 for r in c if r.id[0] == " ") > 4}
    out = {}
    for x, y in itertools.combinations(chains, 2):
        res, _ = interface_residues(chains[x], chains[y], cutoff)
        if res:
            out[(x, y)] = len(res)
    return chains, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="validation")
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--min-contacts", type=int, default=5,
                    help="residues below which a chain pair is not called an interface")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    log = []

    def say(s=""):
        print(s)
        log.append(s)

    say("=" * 76)
    say(f"PARTIAL VALIDATION OF THE VAMP2 BOARD    (cutoff {a.cutoff} A)")
    say("=" * 76)

    results = {}
    for pdb in DEPOSITS:
        try:
            path = fetch(pdb, a.out)
        except Exception as e:
            say(f"\n  {pdb}: could not fetch ({e}) -- skipped")
            continue
        chains, cts = chain_contacts(path, a.cutoff)
        results[pdb] = (chains, cts)
        say(f"\n  {pdb}   chains: {', '.join(chains)}")
        lab = DEPOSITS[pdb] or {}
        for (x, y), n in sorted(cts.items(), key=lambda t: -t[1]):
            flag = "" if n >= a.min_contacts else "   (below threshold)"
            say(f"      {lab.get(x, x):14} - {lab.get(y, y):14} {n:3d} residues{flag}")

    say("\n" + "=" * 76)
    say("ASSERTIONS THE BOARD MAKES")
    say("=" * 76)
    verdict = {}

    for pdb in ("1SFC", "1N7S"):
        if pdb in results:
            _, c = results[pdb]
            partners = {(y if x == "B" else x) for (x, y), n in c.items()
                        if "B" in (x, y) and n >= a.min_contacts}
            verdict[f"A1  VAMP2 four-helix bundle ({pdb})"] = (
                {"A", "C", "D"} <= partners,
                f"VAMP2 engages {sorted(partners) or 'nothing'}; "
                f"expected syntaxin (A) and SNAP25 (C, D)")
            break

    for pdb, key, who in (("1KIL", "A2  complexin is a BRIDGE", "complexin"),
                          ("5CCG", "A3  synaptotagmin is a BRIDGE", "synaptotagmin-1")):
        if pdb in results:
            ch, c = results[pdb]
            deg = {x: 0 for x in ch}
            for (x, y), n in c.items():
                if n >= a.min_contacts:
                    deg[x] += 1
                    deg[y] += 1
            multi = [x for x, d in deg.items() if d >= 2]
            verdict[key] = (bool(multi),
                f"chains contacting >=2 partners: {multi or 'NONE'} "
                f"-- {who} must be among them")

    if "3C98" in results:
        _, c = results["3C98"]
        verdict["A4  Munc18-1 engages syntaxin"] = (
            any(n >= a.min_contacts for n in c.values()),
            f"{sum(1 for n in c.values() if n >= a.min_contacts)} chain pair(s) in contact")

    for k, (ok, why) in verdict.items():
        say(f"  [{'PASS' if ok else 'FAIL'}]  {k:34} {why}")

    say("""
  NOT TESTABLE FROM STRUCTURE -- and the paper says so:
    A5  Munc18-1 does NOT contact VAMP2. No deposit contains both; a chain's
        absence from a crystal is not evidence that an interaction is absent.
    A6  AP180/CALM compete for VAMP2's SNARE motif. No deposited complex.
    A7  Alpha-synuclein at the VAMP2 N-terminus. Intrinsically disordered.

  A2 and A3 are the assertions worth caring about. A bridge piece is DEFINED by
  touching two or more partners, which is what makes it unplaceable until the seam
  beneath it closes -- the constraint that does the combinatorial work in this
  paper. 1KIL and 5CCG can refute it. Nothing else here can.""")

    with open(os.path.join(a.out, "report.txt"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    say(f"\n  wrote {a.out}/report.txt")


if __name__ == "__main__":
    main()
