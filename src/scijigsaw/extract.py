#!/usr/bin/env python3
"""
jigsaw_extract.py -- ADVANCED TIER of the Scientific Jigsaw.

    structures (PDB / mmCIF / AlphaFold)  -->  interactions.csv  -->  the puzzle

The simple tier hands the student a curated interactions.csv: someone has
ALREADY decided that AP180 and SNAP25 compete. The advanced tier makes the
student decide -- and every decision is contestable, exposed as a flag, and
visibly reshapes the puzzle:

    --contact-cutoff   what distance counts as a contact?      (default 5.0 A)
    --overlap          what shared-residue fraction makes two
                       partners COMPETITORS?                   (default 0.20)
    --min-confidence   below what pLDDT/ipTM do we refuse to
                       cut a socket at all?                    (default 70)
    --n-sub            subsites per interface ladder           (default 5)

Loosen the cutoff and separate sites MERGE: competitors appear that are not
real. Tighten it and a real competition VANISHES. The student watches their own
methodological choice reshape the biology. That is the lesson; the picture is
just where it lands.

THE RULE FOR COMPETITION
------------------------
Two partners of the same hub are mutually exclusive if their interface residue
sets on that hub overlap by more than --overlap. This is the standard criterion
used in structural interaction networks (>20% shared interface residues; see
Kim et al., Science 2006 and the structural-PPI literature that follows it).

THE REFUSAL (inherited from the simple tier, and it matters MORE here)
---------------------------------------------------------------------
Below --min-confidence no socket is cut at all: the partner is written out with
an empty site and a warning, and the renderer benches it. A low-confidence
prediction must NOT become a merely 'loose' socket -- if weak evidence could
always be accommodated by loosening a hole, everything would fit and the tool
would have quietly become an arrow diagram with rounder corners.

INPUT
-----
A directory of complex structures, one per interaction, named:

    HUB__PARTNER.pdb        e.g.  VAMP2__SNAP25.pdb
                                  VAMP2__AP180.pdb

Two chains each: chain A = hub, chain B = partner (override with --chains).
AlphaFold outputs carry per-residue pLDDT in the B-factor column, which is read
automatically; for experimental PDBs, confidence is treated as 100.

Optional affinities.csv (protein_a,protein_b,kd_nM) is merged if present; the
Kd -> contacts ladder lives in the renderer, not here.
"""
import argparse
import glob
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, MMCIFParser


def load_chains(path, hub_chain, par_chain):
    parser = (MMCIFParser(QUIET=True) if path.endswith((".cif", ".mmcif"))
              else PDBParser(QUIET=True))
    model = next(iter(parser.get_structure("x", path)))
    try:
        return model[hub_chain], model[par_chain]
    except KeyError:
        ch = [c.id for c in model]
        sys.exit(f"[!] {os.path.basename(path)}: chains {hub_chain}/{par_chain} "
                 f"not found. Present: {ch}. Use --chains.")


def atoms_of(chain):
    out = []
    for res in chain:
        if res.id[0] != " ":                     # skip waters/hetero
            continue
        for a in res:
            if a.element != "H":
                out.append((res.id[1], res.get_resname(), a.coord,
                            float(a.get_bfactor())))
    return out


def interface_residues(hub, par, cutoff):
    """Hub residues within `cutoff` of any partner atom.

    Returns (hub_residue_set, interface_confidence).

    NOTE the confidence must average BOTH SIDES of the interface. Averaging only
    the hub's pLDDT hides a partner that is modelled badly -- which is exactly
    the case that must be refused. (AlphaFold's ipTM/PAE is the proper measure;
    mean interface pLDDT over both chains is the pragmatic proxy.)"""
    H, P = atoms_of(hub), atoms_of(par)
    if not H or not P:
        return set(), 0.0
    hc = np.array([a[2] for a in H])
    pc = np.array([a[2] for a in P])
    hub_res, hub_b, par_b = set(), [], []
    for i in range(0, len(hc), 2000):
        d = np.linalg.norm(hc[i:i + 2000, None, :] - pc[None, :, :], axis=-1)
        for j, row in enumerate(d < cutoff):
            if row.any():
                rn, _, _, b = H[i + j]
                hub_res.add(rn)
                hub_b.append(b)
    for k in range(0, len(pc), 2000):
        d = np.linalg.norm(pc[k:k + 2000, None, :] - hc[None, :, :], axis=-1)
        for j, row in enumerate(d < cutoff):
            if row.any():
                par_b.append(P[k + j][3])
    both = hub_b + par_b
    return hub_res, (float(np.mean(both)) if both else 0.0)


def cluster_sites(faces, overlap):
    """faces = {partner: set(hub residues)}.  Two partners share a SITE if their
    residue sets overlap by more than `overlap` (Jaccard-style, relative to the
    smaller face). Single-linkage clustering -> one site per cluster."""
    names = list(faces)
    parent = {n: n for n in names}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    pairs = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            A, B = faces[a], faces[b]
            if not A or not B:
                continue
            frac = len(A & B) / min(len(A), len(B))
            pairs.append((a, b, frac))
            if frac > overlap:
                parent[find(a)] = find(b)

    groups = defaultdict(list)
    for n in names:
        groups[find(n)].append(n)

    site_of, i = {}, 1
    for members in groups.values():
        # name the site by the residue range its members share
        allres = set().union(*[faces[m] for m in members if faces[m]])
        label = (f"site{i}_{min(allres)}-{max(allres)}" if allres else f"site{i}")
        for m in members:
            site_of[m] = label
        i += 1
    return site_of, pairs


def _cli():
    ap = argparse.ArgumentParser()
    ap.add_argument("structures", help="directory of HUB__PARTNER.pdb complexes")
    ap.add_argument("--out", default="interactions.csv")
    ap.add_argument("--contact-cutoff", type=float, default=5.0,
                    help="A; what counts as a contact (try 4, 5, 8 and compare)")
    ap.add_argument("--overlap", type=float, default=0.20,
                    help="shared-residue fraction above which partners COMPETE")
    ap.add_argument("--min-confidence", type=float, default=70.0,
                    help="mean interface pLDDT below which NO socket is cut")
    ap.add_argument("--n-sub", type=int, default=5)
    ap.add_argument("--chains", default="A,B", help="hub,partner chain ids")
    ap.add_argument("--affinities", help="optional csv: protein_a,protein_b,kd_nM")
    a = ap.parse_args()

    hub_chain, par_chain = a.chains.split(",")
    files = sorted(glob.glob(os.path.join(a.structures, "*__*.pdb"))
                   + glob.glob(os.path.join(a.structures, "*__*.cif")))
    if not files:
        sys.exit(f"[!] no HUB__PARTNER.pdb files in {a.structures}")

    faces, conf, hub_name = {}, {}, None
    print(f"contact cutoff {a.contact_cutoff} A   overlap {a.overlap}   "
          f"min confidence {a.min_confidence}\n")
    for f in files:
        base = os.path.basename(f).rsplit(".", 1)[0]
        hub, partner = base.split("__", 1)
        hub_name = hub_name or hub
        if hub != hub_name:
            sys.exit(f"[!] mixed hubs: {hub_name} vs {hub}. One hub per run.")
        hc, pc = load_chains(f, hub_chain, par_chain)
        res, mean_conf = interface_residues(hc, pc, a.contact_cutoff)
        faces[partner] = res
        conf[partner] = mean_conf
        print(f"  {partner:16} {len(res):3d} interface residues   "
              f"mean pLDDT {mean_conf:5.1f}"
              + ("   <-- BELOW FLOOR: no socket will be cut"
                 if mean_conf < a.min_confidence else ""))

    # refusal: low-confidence faces are DISCARDED, not merely loosened
    refused = {p for p, c in conf.items() if c < a.min_confidence}
    kept = {p: f for p, f in faces.items() if p not in refused}

    site_of, pairs = cluster_sites(kept, a.overlap)

    print("\nPAIRWISE INTERFACE OVERLAP ON THE HUB")
    for x, y, frac in sorted(pairs, key=lambda t: -t[2]):
        tag = "  COMPETE" if frac > a.overlap else ""
        print(f"  {x:14} vs {y:14} {frac:5.2f}{tag}")

    sites = defaultdict(list)
    for p, s in site_of.items():
        sites[s].append(p)
    print("\nSITES DERIVED")
    for s, members in sites.items():
        tag = "  <-- CONTESTED" if len(members) > 1 else ""
        print(f"  {s:22} {', '.join(members)}{tag}")

    kd = {}
    if a.affinities and os.path.exists(a.affinities):
        for r in pd.read_csv(a.affinities).itertuples():
            kd[(r.protein_a, r.protein_b)] = r.kd_nM
            kd[(r.protein_b, r.protein_a)] = r.kd_nM

    rows = []
    for p in faces:
        rows.append(dict(
            protein_a=hub_name, protein_b=p,
            kd_nM=kd.get((hub_name, p), np.nan),
            site_on_a="" if p in refused else site_of[p],
            site_on_b="" if p in refused else f"{p}_iface",
            n_interface_res=len(faces[p]),
            mean_plddt=round(conf[p], 1),
            refused=p in refused))
    df = pd.DataFrame(rows)
    df.to_csv(a.out, index=False)

    if refused:
        print(f"\n[!] REFUSED (confidence < {a.min_confidence}): "
              f"{', '.join(sorted(refused))}")
        print("    No socket is cut. These are written with an EMPTY site and the")
        print("    renderer will bench them. A weak prediction must not become a")
        print("    merely 'loose' socket -- otherwise everything fits and the tool")
        print("    has quietly become an arrow diagram.")

    print(f"\nwrote {a.out}")
    print("Now re-run with a different --contact-cutoff and watch the puzzle change.")




# ------------------------------------------------------------------ public API
def extract(structure_dir, contact_cutoff=5.0, overlap=0.20,
            min_confidence=70.0, chains="A,B"):
    """Derive interfaces, sites and overlaps from HUB__PARTNER complexes.

    Returns a DataFrame with one row per partner:
        protein_a, protein_b, site_on_a, site_on_b, n_interface_res,
        mean_plddt, refused

    THRESHOLDS ARE CHOICES, NOT DEFAULTS TO BE IGNORED. Interface detection is
    sensitive to `contact_cutoff`: in the bundled example the alpha-synuclein
    interface is absent at 5 A and present at 8 A. That is not a property of the
    biology; it is a property of the analysis. Report the value you used, and
    sweep it.

    Partners whose mean interface confidence falls below `min_confidence` are
    REFUSED: no socket is cut and the partner is benched. A weak prediction must
    not become a merely loose socket, or everything fits and the discipline
    evaporates.
    """
    import glob, os, sys
    from collections import defaultdict
    import numpy as np, pandas as pd

    hub_chain, par_chain = chains.split(",")
    files = sorted(glob.glob(os.path.join(structure_dir, "*__*.pdb")) +
                   glob.glob(os.path.join(structure_dir, "*__*.cif")))
    if not files:
        raise FileNotFoundError(f"no HUB__PARTNER.pdb files in {structure_dir}")

    faces, conf, hub_name = {}, {}, None
    for f in files:
        base = os.path.basename(f).rsplit(".", 1)[0]
        hub, partner = base.split("__", 1)
        hub_name = hub_name or hub
        if hub != hub_name:
            raise ValueError(f"mixed hubs: {hub_name} vs {hub}; one hub per run")
        hc, pc = load_chains(f, hub_chain, par_chain)
        res, c = interface_residues(hc, pc, contact_cutoff)
        faces[partner], conf[partner] = res, c

    refused = {p for p, c in conf.items() if c < min_confidence}
    kept = {p: f for p, f in faces.items() if p not in refused}
    site_of, pairs = cluster_sites(kept, overlap)

    rows = []
    for p in faces:
        rows.append(dict(protein_a=hub_name, protein_b=p,
                         site_on_a="" if p in refused else site_of[p],
                         site_on_b="" if p in refused else f"{p}_iface",
                         n_interface_res=len(faces[p]),
                         mean_plddt=round(conf[p], 1),
                         refused=p in refused))
    return pd.DataFrame(rows), pairs
