"""
run_cryoem_real.py -- full pipeline on the deposited Mark et al. 2026 structures.

Downloads six yeast proteasome precursor mmCIF files, extracts chain-pair
interfaces at the declared threshold (5 A heavy-atom, >=3 qualifying residues),
derives each structure's subunit composition, and runs the exact assembly
analysis: which observed intermediates are admissible, and how many distinct
mechanisms are required to explain them all.

No dependencies beyond the standard library and sj_prob / sj_navigate.

Usage:
    python3 run_cryoem_real.py            # downloads on first run, then caches
    python3 run_cryoem_real.py --cutoff 4.5 --minres 5
"""

from __future__ import annotations

import argparse
import gzip
import itertools
import os
import sys
import urllib.request
from collections import defaultdict
from typing import Dict, FrozenSet, List, Set, Tuple

# Mark E, Ramos PC, Nunes MM, Matias AC, Dohmen RJ, Wendler P.
# Structural transitions in the stepwise assembly of proteasome core particles.
# Nat Commun 17 (2026). doi:10.1038/s41467-026-70525-w
ENTRIES = {
    "9rl3": ("13S-PC",            "EMD-54029", 3.31),
    "9rla": ("13S+b1-PC",         "EMD-54032", 3.16),
    "9rm0": ("13S+b5+b6-PC",      "EMD-54047", 3.29),
    "9rm1": ("13S+b1+b5-PC",      "EMD-54048", 4.11),
    "9rlz": ("15S-PC",            "EMD-54046", 3.12),
    "9rlt": ("13S-13S+b5 dimer",  "EMD-54045", 3.09),
}
CACHE = "cryoem_cache"
URLS = ["https://files.rcsb.org/download/{pdb}.cif",
        "https://files.wwpdb.org/pub/pdb/data/structures/divided/mmCIF/"
        "{mid}/{pdb}.cif.gz",
        "https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb}.cif"]

# Yeast 20S beta subunits of interest, by common gene/systematic naming.
# The extractor matches these case-insensitively against entity descriptions.
SUBUNIT_ALIASES = {
    "b1": ["pre3", "beta1", "beta-1", "psb6", "subunit beta type-1"],
    "b2": ["pup1", "beta2", "beta-2", "subunit beta type-2"],
    "b3": ["pup3", "beta3", "beta-3", "subunit beta type-3"],
    "b4": ["pre1", "beta4", "beta-4", "subunit beta type-4"],
    "b5": ["pre2", "doa3", "beta5", "beta-5", "subunit beta type-5"],
    "b6": ["pre7", "beta6", "beta-6", "subunit beta type-6"],
    "b7": ["pre4", "beta7", "beta-7", "subunit beta type-7"],
    "Ump1": ["ump1"],
    "Pba1": ["pba1", "poc1"],
    "Pba2": ["pba2", "poc2", "add66"],
}
ALPHA_HINTS = ["alpha", "scl1", "pre8", "pre9", "pre6", "pup2", "pre5", "pre10",
               "subunit alpha type"]


# ---------------------------------------------------------------- download
def fetch(pdb: str) -> str:
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, f"{pdb}.cif")
    if os.path.exists(path) and os.path.getsize(path) > 10000:
        return path
    for tmpl in URLS:
        url = tmpl.format(pdb=pdb, mid=pdb[1:3])
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                raw = r.read()
            if url.endswith(".gz"):
                raw = gzip.decompress(raw)
            if len(raw) < 10000:
                continue
            with open(path, "wb") as fh:
                fh.write(raw)
            print(f"    downloaded {pdb} ({len(raw)//1024} kB) from "
                  f"{url.split('/')[2]}")
            return path
        except Exception as exc:                                  # noqa: BLE001
            last = exc
    raise RuntimeError(f"could not download {pdb}: {last}")


# ------------------------------------------------------- minimal mmCIF parse
def parse_cif(path: str):
    """Return (atoms, chain_desc).

    atoms: list of (chain, resseq, x, y, z) for heavy atoms of polymers.
    chain_desc: chain id -> entity description string (lower-case).
    """
    entity_desc: Dict[str, str] = {}
    chain_entity: Dict[str, str] = {}
    atoms: List[Tuple[str, str, float, float, float]] = []

    with open(path, "r", errors="ignore") as fh:
        lines = fh.readlines()

    # --- entity descriptions (loop or single-item form) ---------------------
    i = 0
    while i < len(lines):
        ln = lines[i].strip()
        if ln == "loop_":
            j, tags = i + 1, []
            while j < len(lines) and lines[j].lstrip().startswith("_"):
                tags.append(lines[j].strip()); j += 1
            if any(t.startswith("_entity.") for t in tags):
                idx_id = tags.index("_entity.id") if "_entity.id" in tags else None
                idx_de = (tags.index("_entity.pdbx_description")
                          if "_entity.pdbx_description" in tags else None)
                while j < len(lines) and not lines[j].startswith(("#", "loop_")):
                    row = split_cif_row(lines[j])
                    if idx_id is not None and idx_de is not None and len(row) > max(idx_id, idx_de):
                        entity_desc[row[idx_id]] = row[idx_de].lower()
                    j += 1
            i = j
            continue
        if ln.startswith("_entity.id"):
            eid = ln.split(None, 1)[1].strip() if len(ln.split()) > 1 else None
            k = i
            while k < len(lines) and not lines[k].startswith("#"):
                if lines[k].strip().startswith("_entity.pdbx_description") and eid:
                    entity_desc[eid] = lines[k].split(None, 1)[1].strip().strip("'\"").lower()
                k += 1
        i += 1

    # --- atom_site loop -----------------------------------------------------
    i = 0
    while i < len(lines):
        if lines[i].strip() == "loop_":
            j, tags = i + 1, []
            while j < len(lines) and lines[j].lstrip().startswith("_"):
                tags.append(lines[j].strip()); j += 1
            if any(t.startswith("_atom_site.") for t in tags):
                col = {t: k for k, t in enumerate(tags)}
                need = ["_atom_site.label_asym_id", "_atom_site.Cartn_x",
                        "_atom_site.Cartn_y", "_atom_site.Cartn_z",
                        "_atom_site.type_symbol", "_atom_site.label_seq_id",
                        "_atom_site.label_entity_id", "_atom_site.group_PDB"]
                if not all(n in col for n in need):
                    i = j; continue
                while j < len(lines) and not lines[j].startswith(("#", "loop_")):
                    row = split_cif_row(lines[j])
                    if len(row) < len(tags):
                        j += 1; continue
                    if row[col["_atom_site.group_PDB"]] != "ATOM":
                        j += 1; continue
                    if row[col["_atom_site.type_symbol"]].upper() == "H":
                        j += 1; continue
                    ch = row[col["_atom_site.label_asym_id"]]
                    chain_entity.setdefault(ch, row[col["_atom_site.label_entity_id"]])
                    try:
                        atoms.append((ch, row[col["_atom_site.label_seq_id"]],
                                      float(row[col["_atom_site.Cartn_x"]]),
                                      float(row[col["_atom_site.Cartn_y"]]),
                                      float(row[col["_atom_site.Cartn_z"]])))
                    except ValueError:
                        pass
                    j += 1
            i = j
            continue
        i += 1

    chain_desc = {ch: entity_desc.get(ent, f"entity_{ent}")
                  for ch, ent in chain_entity.items()}
    return atoms, chain_desc


def split_cif_row(line: str) -> List[str]:
    out, cur, quote = [], "", None
    for ch in line.rstrip("\n"):
        if quote:
            if ch == quote:
                out.append(cur); cur = ""; quote = None
            else:
                cur += ch
        elif ch in "'\"":
            quote = ch
        elif ch.isspace():
            if cur:
                out.append(cur); cur = ""
        else:
            cur += ch
    if cur:
        out.append(cur)
    return out


# --------------------------------------------------------- interface extract
def interfaces(atoms, cutoff: float, minres: int) -> Set[FrozenSet[str]]:
    """Chain pairs with >= minres residues within cutoff (grid-accelerated)."""
    cell = cutoff
    grid: Dict[Tuple[int, int, int], List[int]] = defaultdict(list)
    for idx, (_, _, x, y, z) in enumerate(atoms):
        grid[(int(x // cell), int(y // cell), int(z // cell))].append(idx)

    pair_res: Dict[FrozenSet[str], Set[Tuple[str, str, str]]] = defaultdict(set)
    c2 = cutoff * cutoff
    for key, idxs in grid.items():
        gx, gy, gz = key
        neigh: List[int] = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    neigh.extend(grid.get((gx + dx, gy + dy, gz + dz), ()))
        for a in idxs:
            ca, ra, xa, ya, za = atoms[a]
            for b in neigh:
                if b <= a:
                    continue
                cb, rb, xb, yb, zb = atoms[b]
                if ca == cb:
                    continue
                if (xa - xb) ** 2 + (ya - yb) ** 2 + (za - zb) ** 2 <= c2:
                    key2 = frozenset({ca, cb})
                    pair_res[key2].add((ca, ra, "a"))
                    pair_res[key2].add((cb, rb, "b"))
    return {p for p, rs in pair_res.items() if len(rs) >= 2 * minres}


# ------------------------------------------------------------- name mapping
def label_chain(desc: str) -> str:
    d = desc.lower()
    for name, aliases in SUBUNIT_ALIASES.items():
        if any(a in d for a in aliases):
            return name
    if any(a in d for a in ALPHA_HINTS):
        return "alpha"
    return "other"


# --------------------------------------------------------------------- main
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--minres", type=int, default=3)
    args = ap.parse_args()

    print("=" * 78)
    print("REAL cryo-EM extraction: Mark et al. 2026 proteasome precursors")
    print(f"threshold: {args.cutoff} A heavy-atom, >= {args.minres} qualifying residues")
    print("=" * 78)

    comps: Dict[str, Set[str]] = {}
    contacts: Set[FrozenSet[str]] = set()
    per_entry_contacts: Dict[str, Set[FrozenSet[str]]] = {}

    for pdb, (name, emd, res) in ENTRIES.items():
        print(f"\n  {pdb}  {name}  ({emd}, {res} A)")
        path = fetch(pdb)
        atoms, desc = parse_cif(path)
        chain_lab = {ch: label_chain(d) for ch, d in desc.items()}
        present = {v for v in chain_lab.values() if v not in ("other",)}
        comps[pdb] = present
        print(f"    {len(atoms)} heavy atoms, {len(chain_lab)} chains")
        print(f"    composition: {' '.join(sorted(present))}")
        ifs = interfaces(atoms, args.cutoff, args.minres)
        named = {frozenset({chain_lab[a], chain_lab[b]})
                 for a, b in (tuple(p) for p in ifs)
                 if chain_lab.get(a) and chain_lab.get(b)
                 and chain_lab[a] != chain_lab[b]
                 and "other" not in (chain_lab[a], chain_lab[b])}
        per_entry_contacts[pdb] = named
        contacts |= named
        print(f"    {len(ifs)} chain-pair interfaces -> {len(named)} named contacts")

    print("\n" + "=" * 78)
    print("EXTRACTED CONTACT SET (union over all six structures)")
    print("=" * 78)
    for c in sorted(contacts, key=lambda s: sorted(s)):
        u, v = sorted(c)
        seen = [p for p in ENTRIES if c in per_entry_contacts[p]]
        print(f"  {u:>6} -- {v:<6}   in {len(seen)}/6 structures: {' '.join(seen)}")

    # ---- assembly analysis on the extracted model -------------------------
    from sj_prob import TypedModel
    from sj_navigate import Navigator

    base = comps["9rl3"]                      # 13S core = the seed set
    variable = sorted(set().union(*comps.values()) - base)
    print("\n" + "=" * 78)
    print("ASSEMBLY ANALYSIS")
    print("=" * 78)
    print(f"  seed (13S composition): {' '.join(sorted(base))}")
    print(f"  variable subunits:      {' '.join(variable)}")

    SEED = "13S"
    V = tuple([SEED] + variable)
    C = frozenset(frozenset({SEED, v}) for v in variable) | \
        frozenset(c for c in contacts if c <= set(variable))
    model = TypedModel(V=V, C=C, P=())
    nav = Navigator(model, SEED)
    print(f"  contacts used: {len(C)}   admissible mechanisms: {nav.n_histories()}")

    observed = {p: (comps[p] - base) | {SEED} for p in ENTRIES
                if comps[p] != base and "dimer" not in ENTRIES[p][0]}
    print("\n  observed intermediates as queries:")
    for p, E in observed.items():
        r = nav.ask(nav.forms_exact(E))
        n = round(r["prob"] * nav.n_histories())
        print(f"    {p} {ENTRIES[p][0]:<16} {sorted(E)}  "
              f"{n}/{nav.n_histories()} = {r['prob']:.3f}")

    hist = [h for h, _ in nav.masses]
    cover = {}
    for h in hist:
        got = set()
        for p, E in observed.items():
            cur = {SEED}
            if cur == E:
                got.add(p)
            for x in h:
                cur.add(x)
                if cur == E:
                    got.add(p)
        cover[h] = got
    best = max((len(v) for v in cover.values()), default=0)
    need = None
    for k in range(1, len(hist) + 1):
        for combo in itertools.combinations(hist, k):
            if set().union(*(cover[c] for c in combo)) == set(observed):
                need = (k, combo); break
        if need:
            break
    print(f"\n  best single mechanism explains {best}/{len(observed)} observations")
    if need:
        k, combo = need
        print(f"  MINIMUM distinct mechanisms required: {k}")
        for c in combo:
            print(f"    {SEED}->{'->'.join(c)}   covers "
                  f"{' '.join(sorted(cover[c]))}")


if __name__ == "__main__":
    main()
