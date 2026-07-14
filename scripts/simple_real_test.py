#!/usr/bin/env python3
"""The simplest real test that could refute the paper.

    python scripts/simple_real_test.py          # ~30 seconds, needs network

ONE STRUCTURE. ONE QUESTION.

The paper's combinatorial result rests on BRIDGE PIECES. Complexin is drawn with two
sockets, spanning the seam between VAMP2 and SNAP25, and therefore cannot be seated
until that seam is closed. That single constraint is a large part of why the VAMP2
board eliminates 93% of assembly orders.

The claim is structural and it is falsifiable: if complexin is a bridge, then in a
deposited complexin-SNARE complex it must contact TWO OR MORE chains of the bundle.
If it contacts only one, it is an ordinary single-socket partner, the bridge constraint
is wrong, and the pruning figure is wrong with it.

1KIL is the complexin/SNARE complex (Chen et al., Neuron 2002). We ask the structure.

WHAT A PASS MEANS
    The bridge constraint is supported by structure for this piece. Not that the whole
    board is validated -- see validate_board_from_structures.py for the fuller set, and
    note that three of the board's seven assertions cannot be tested from structure at
    all (Munc18-not-VAMP2, AP180/CALM competition, alpha-synuclein at the N-terminus).

WHAT A FAIL MEANS
    Better you find it than a referee. Report it.
"""
import itertools
import os
import sys
import urllib.request

PDB = "1KIL"
CUTOFF = 5.0            # angstroms; the standard heavy-atom contact distance
MIN_RES = 5             # below this, a chain pair is noise, not an interface


def main():
    path = f"{PDB}.pdb"
    if not os.path.exists(path):
        print(f"fetching {PDB} from RCSB ...")
        try:
            urllib.request.urlretrieve(
                f"https://files.rcsb.org/download/{PDB}.pdb", path)
        except Exception as e:
            sys.exit(f"could not fetch {PDB}: {e}\n"
                     f"Download it by hand from https://www.rcsb.org/structure/{PDB} "
                     f"and put {PDB}.pdb here.")

    from Bio.PDB import PDBParser
    from scijigsaw.extract import interface_residues

    # which chain is which? read it from the file rather than assuming
    names = {}
    cur = None
    for line in open(path):
        if line.startswith("COMPND"):
            t = line[10:].strip()
            if t.startswith("MOLECULE:"):
                cur = t.split(":", 1)[1].strip().rstrip(";")
            elif t.startswith("CHAIN:") and cur:
                for c in t.split(":", 1)[1].strip().rstrip(";").split(","):
                    names[c.strip()] = cur

    model = next(iter(PDBParser(QUIET=True).get_structure("x", path)))
    chains = {c.id: c for c in model if sum(1 for r in c if r.id[0] == " ") > 4}

    print(f"\n{PDB}: {len(chains)} chains")
    for c in chains:
        print(f"    {c}  {names.get(c, '(unnamed)')}")

    print(f"\ninterfaces (heavy-atom contacts within {CUTOFF} A):")
    degree = {c: 0 for c in chains}
    contacts_of = {c: set() for c in chains}
    for x, y in itertools.combinations(chains, 2):
        res, _ = interface_residues(chains[x], chains[y], CUTOFF)
        if len(res) >= MIN_RES:
            degree[x] += 1
            degree[y] += 1
            contacts_of[x].add(y)
            contacts_of[y].add(x)
            print(f"    {x} - {y}   {len(res):3d} residues")

    print("\npartners contacted by each chain:")
    for c, d in degree.items():
        tag = "  <-- BRIDGE (>=2 partners)" if d >= 2 else ""
        print(f"    {c}  {names.get(c, '?')[:34]:34} {d}{tag}")

    # the verdict
    cpx = [c for c, n in names.items()
           if "COMPLEXIN" in n.upper() and c in chains]
    print("\n" + "=" * 66)
    if not cpx:
        print("  Could not identify the complexin chain from COMPND. Read the chain")
        print("  list above and see which chains complexin contacts.")
    else:
        c = cpx[0]
        partners = sorted(p for p in contacts_of[c])
        pn = [names.get(p, p) for p in partners]
        print(f"  complexin is chain {c}, and contacts {len(partners)} chain(s):")
        for p in partners:
            print(f"      {p}  {names.get(p, '?')}")
        print()
        up = " ".join(pn).upper()
        has_vamp = "SYNAPTOBREVIN" in up or "VAMP" in up
        has_stx = "SYNTAXIN" in up
        has_snap = "SNAP-25" in up or "SNAP25" in up
        if len(partners) < 2:
            print("  [FAIL] complexin is NOT a bridge. The precedence encoding is wrong,")
            print("         the permitted-order count is wrong, and 93% must be recomputed.")
        elif has_vamp and has_stx and not has_snap:
            print("  [!] The structure supports  Complexin <- {VAMP2, Syntaxin-1A},")
            print("      NOT {VAMP2, SNAP25} as encoded in cases.py.")
            print("      -> permitted orders become 252, not 336; 95.0% eliminated, not 93.3%.")
            print("      Fix cases.py, rerun scripts/reproduce_numbers.py, and correct the")
            print("      abstract, Results, Figure 3 and the README.")
        elif has_vamp and has_snap:
            print("  [PASS] The structure supports the encoding as written:")
            print("         Complexin <- {VAMP2, SNAP25}. 336 / 5040 stands.")
        else:
            print("  [?] complexin bridges, but not the pair encoded. Read the contacts")
            print("      above and re-encode cases.py to match. Then recompute.")
    print("=" * 66)
    print("\n  This tests ONE assertion on ONE structure. It does not validate the board.")
    print("  Three of the board's seven assertions cannot be tested from structure at")
    print("  all. See scripts/validate_board_from_structures.py.")


if __name__ == "__main__":
    main()
