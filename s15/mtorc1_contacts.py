#!/usr/bin/env python3
"""
Step 1 of the mTORC1 / rapamycin example.

Extracts the mTORC1 contact graph from deposited coordinates using the SAME
criterion as the rest of the manuscript (heavy-atom separation within 5 A,
minimum three unique interface residues per partner), and writes the Scientific
Jigsaw protein and interaction tables.

Structures
----------
  6BCX  apo mTORC1, 3.0 A (Yang et al. 2017)        -> M0 contacts
  6BCU  RHEB-mTORC1, 3.4 A (Yang et al. 2017)       -> sanity check on M0
  5FLC  mTORC1 with FKBP12-rapamycin, 5.9 A
          (Aylett et al. 2016)                       -> FKBP12-FRB contact
  1FAP  FKBP12-rapamycin-FRB ternary, 2.7 A          -> FKBP12-FRB at high res

VERIFY BEFORE USE: confirm each PDB ID resolves to the entry described above and
that the chain->subunit mapping below matches the deposited entity records. Do
NOT assume chain identities from deposition order -- derive them from entities,
as done for the proteasome case in Section S14.2.

Usage
-----
    python mtorc1_contacts.py --outdir boards/mtorc1
"""
import argparse
import csv
import os
import sys
import urllib.request

CIF_URL = "https://files.rcsb.org/download/{}.cif"

# ---------------------------------------------------------------------------
# Chain -> subunit mapping. PLACEHOLDER: derive from entity records and verify
# against residue counts before trusting. Left explicit so the mapping is
# auditable rather than assumed.
# ---------------------------------------------------------------------------
CHAIN_MAP = {
    "6BCX": {},   # e.g. {"A": "mTOR", "B": "RAPTOR", "C": "mLST8", ...}
    "6BCU": {},
    "5FLC": {},
    "1FAP": {},   # FKBP12 / FRB / rapamycin ligand
}

CUTOFF = 5.0        # angstrom, heavy atoms
MIN_RESIDUES = 3    # unique interface residues required on each partner


def fetch(pdb_id, cache="cif"):
    os.makedirs(cache, exist_ok=True)
    path = os.path.join(cache, f"{pdb_id}.cif")
    if not os.path.exists(path):
        urllib.request.urlretrieve(CIF_URL.format(pdb_id), path)
    return path


def extract(path, chain_map):
    """Return {(subunitA, subunitB): (nA, nB)} for pairs passing the criterion.

    Uses the project extractor so that this example is subject to exactly the
    same rules -- and the same PDBePISA-validated behaviour -- as every other
    structural claim in the manuscript.
    """
    try:
        from scijigsaw.extract import interfaces_from_structure
    except ImportError:
        sys.exit(
            "scijigsaw not importable. Install the package (pip install -e .) so "
            "that this example uses the validated extractor rather than a "
            "reimplementation."
        )
    return interfaces_from_structure(
        path,
        chain_map=chain_map,
        contact_cutoff=CUTOFF,
        min_interface_residues=MIN_RESIDUES,
        heavy_atoms_only=True,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="boards/mtorc1")
    ap.add_argument("--structures", nargs="+",
                    default=["6BCX", "6BCU", "5FLC", "1FAP"])
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    rows = []
    for pdb_id in args.structures:
        cmap = CHAIN_MAP.get(pdb_id) or {}
        if not cmap:
            print(f"[skip] {pdb_id}: chain map empty -- fill CHAIN_MAP first",
                  file=sys.stderr)
            continue
        path = fetch(pdb_id)
        for (a, b), (na, nb) in sorted(extract(path, cmap).items()):
            rows.append({
                "component_a": a,
                "component_b": b,
                "relation": "contact",
                "residues_a": na,
                "residues_b": nb,
                "evidence": f"PDB {pdb_id}",
                "provenance": f"extracted, {CUTOFF} A, >= {MIN_RESIDUES} residues",
            })
            print(f"{pdb_id}  {a:<8} {b:<8} {na}+{nb}")

    out = os.path.join(args.outdir, "interactions_extracted.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else
                           ["component_a", "component_b", "relation",
                            "residues_a", "residues_b", "evidence", "provenance"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {len(rows)} extracted contacts -> {out}")
    print("\nNEXT: hand-check every row against the deposited entity records, then "
          "add the two declared relations that are NOT extracted from coordinates:\n"
          "  - FKBP12 contact is ligand-conditioned on rapamycin (declared)\n"
          "  - mTOR-RAPTOR exclusion under rapamycin (declared, cited to the\n"
          "    disruption literature, NOT to a structure)")


if __name__ == "__main__":
    main()
