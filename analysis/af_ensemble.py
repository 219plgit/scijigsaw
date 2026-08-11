"""Assembly-space analysis over an AlphaFold 3 prediction ensemble.

Turns a set of predicted structures into a distribution over admissible
assembly spaces, separating STRUCTURAL uncertainty (does the interface appear?)
from MECHANISTIC uncertainty (given the interfaces, which histories are
admissible?).

Quantities reported
-------------------
q_uv     interface recurrence: fraction of samples in which u-v qualifies
R(E)     ensemble robustness: fraction of samples in which subcomplex E is
         formable at all
G_E/F    within-model support: fraction of admissible trees containing E,
         computed per sample and summarised across the ensemble
decisive an interface whose presence/absence changes whether E is formable

Usage
-----
    python af_ensemble.py <dir-with-cif-and-json> --chains a=A b=B ...
"""
from __future__ import annotations
import argparse, glob, json, os, re, sys
from collections import defaultdict
import numpy as np

CUTOFF, MIN_RES = 5.0, 3


# ------------------------------------------------------------------ parsing
def parse_cif(path):
    """Return {auth_chain: (coords Nx3, resids)} for heavy protein atoms.

    Reads atom_site directly: the label/auth chain-id distinction has caused
    silent chain loss with library parsers, so column positions are resolved
    from the header rather than assumed.
    """
    cols, in_loop, idx = [], False, {}
    co, ri = defaultdict(list), defaultdict(list)
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                cols.append(line.strip().split(".")[-1]); in_loop = True; continue
            if in_loop and not idx:
                idx = {c: i for i, c in enumerate(cols)}
                need = ["type_symbol", "Cartn_x", "Cartn_y", "Cartn_z",
                        "auth_asym_id", "auth_seq_id"]
                missing = [n for n in need if n not in idx]
                if missing:
                    raise SystemExit(f"{path}: missing columns {missing}")
            if line.startswith(("ATOM", "HETATM")):
                f = line.split()
                if len(f) <= max(idx.values()): continue
                if f[idx["type_symbol"]] == "H": continue
                c = f[idx["auth_asym_id"]]
                co[c].append((float(f[idx["Cartn_x"]]), float(f[idx["Cartn_y"]]),
                              float(f[idx["Cartn_z"]])))
                ri[c].append(f[idx["auth_seq_id"]])
    return {c: (np.array(v), ri[c]) for c, v in co.items()}


def interfaces(chains, cutoff=CUTOFF, min_res=MIN_RES):
    """{(u,v): (n_u, n_v)} for pairs meeting the criterion."""
    out, ids = {}, sorted(chains)
    for i, u in enumerate(ids):
        for v in ids[i + 1:]:
            A, ra = chains[u]; B, rb = chains[v]
            su, sv = set(), set()
            for k in range(len(A)):
                d = np.linalg.norm(B - A[k], axis=1)
                if (d <= cutoff).any():
                    su.add(ra[k])
                    for j in np.where(d <= cutoff)[0]: sv.add(rb[j])
            if len(su) >= min_res or len(sv) >= min_res:
                out[(u, v)] = (len(su), len(sv))
    return out


def read_confidence(path):
    """Pull whatever confidence summary is present; tolerate schema variation."""
    try:
        d = json.load(open(path))
    except Exception:
        return {}
    keep = {}
    for k in ("iptm", "ptm", "ranking_score", "fraction_disordered",
              "has_clash", "chain_iptm", "chain_pair_iptm"):
        if k in d: keep[k] = d[k]
    return keep


# ------------------------------------------------- enumeration over a sample
def trees_and_support(units, contacts, prereqs, queries):
    """Admissible assembly trees for one contact graph, plus support of each
    queried subset. Uses the same recurrence as the manuscript."""
    idx = {u: i for i, u in enumerate(units)}
    n = len(units)
    con = [0] * n
    for (u, v) in contacts:
        if u in idx and v in idx:
            con[idx[u]] |= 1 << idx[v]; con[idx[v]] |= 1 << idx[u]
    pre = [0] * n
    for (u, v) in prereqs:                       # u must precede v
        if u in idx and v in idx: pre[idx[v]] |= 1 << idx[u]

    def connected(S):
        if S == 0: return False
        f = (S & -S).bit_length() - 1; seen = 1 << f; st = [f]
        while st:
            i = st.pop(); nb = con[i] & S & ~seen
            while nb:
                j = (nb & -nb).bit_length() - 1
                seen |= 1 << j; st.append(j); nb &= nb - 1
        return seen == S

    def merge_ok(A, B):
        C = A | B; linked = False; m = A
        while m:
            i = (m & -m).bit_length() - 1
            if con[i] & B: linked = True
            m &= m - 1
        if not linked: return False
        m = C
        while m:
            i = (m & -m).bit_length() - 1
            if pre[i] & ~C: return False
            m &= m - 1
        return True

    memoF = {}
    def F(S):
        if S & (S - 1) == 0: return 1
        v = memoF.get(S)
        if v is not None: return v
        tot = 0; low = S & -S; sub = (S - 1) & S
        while sub:
            if sub & low:
                A, B = sub, S ^ sub
                if B and connected(A) and connected(B) and merge_ok(A, B):
                    tot += F(A) * F(B)
            sub = (sub - 1) & S
        memoF[S] = tot; return tot

    def G(E, S, memoG):
        if S == E: return F(S)
        v = memoG.get(S)
        if v is not None: return v
        if E & ~S or S & (S - 1) == 0:
            memoG[S] = 0; return 0
        tot = 0; low = S & -S; sub = (S - 1) & S
        while sub:
            if sub & low:
                A, B = sub, S ^ sub
                if B and connected(A) and connected(B) and merge_ok(A, B):
                    if not (E & ~A):   tot += G(E, A, memoG) * F(B)
                    elif not (E & ~B): tot += F(A) * G(E, B, memoG)
            sub = (sub - 1) & S
        memoG[S] = tot; return tot

    full = (1 << n) - 1
    total = F(full)
    res = {}
    for name, subset in queries.items():
        E = sum(1 << idx[u] for u in subset if u in idx)
        k = G(E, full, {}) if (total and connected(E)) else 0
        res[name] = (k, total)
    return total, res


# ------------------------------------------------------------------- driver
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("directory")
    ap.add_argument("--chains", nargs="+", required=True,
                    help="label=chainID, e.g. b=A i=B g=C a=D c=E")
    ap.add_argument("--prereq", nargs="*", default=[],
                    help="temporal relations u:v meaning u precedes v")
    ap.add_argument("--query", nargs="*", default=[],
                    help="subcomplexes to test, e.g. g:i b:i b:g:i")
    ap.add_argument("--cutoff", type=float, default=CUTOFF)
    a = ap.parse_args()

    lab2chain = dict(x.split("=") for x in a.chains)
    chain2lab = {v: k for k, v in lab2chain.items()}
    units = sorted(lab2chain)
    prereqs = [tuple(x.split(":")) for x in a.prereq]
    queries = {q: q.split(":") for q in a.query}

    cifs = sorted(glob.glob(os.path.join(a.directory, "*.cif")))
    if not cifs: sys.exit(f"no .cif files in {a.directory}")
    print(f"{len(cifs)} predictions, units {units}\n")

    recur = defaultdict(int); per_sample = []
    for c in cifs:
        chains = parse_cif(c)
        missing = [ch for ch in lab2chain.values() if ch not in chains]
        if missing:
            print(f"  ! {os.path.basename(c)}: chains {missing} absent "
                  f"(found {sorted(chains)}) -- skipped"); continue
        sub = {chain2lab[ch]: chains[ch] for ch in lab2chain.values()}
        ifs = interfaces(sub, a.cutoff)
        edges = set(ifs)
        for e in edges: recur[tuple(sorted(e))] += 1
        tot, res = trees_and_support(units, edges, prereqs, queries)
        j = os.path.splitext(c)[0] + ".json"
        conf = read_confidence(j) if os.path.exists(j) else {}
        per_sample.append((os.path.basename(c), edges, tot, res, conf))
        note = ""
        if tot == 0:
            iso = [u for u in units
                   if not any(u in e for e in edges)]
            note = ("  <- no admissible tree: "
                    + (f"units {iso} have no interface" if iso
                       else "contact graph disconnected"))
        print(f"  {os.path.basename(c):40} {len(edges):2} interfaces  "
              f"{tot:>6} trees"
              + (f"  ipTM={conf.get('iptm'):.2f}" if isinstance(conf.get('iptm'), float) else "")
              + note)

    n = len(per_sample)
    if not n: sys.exit("no usable predictions")

    print(f"\n=== interface recurrence q_uv over {n} samples ===")
    for e, k in sorted(recur.items(), key=lambda t: -t[1]):
        print(f"   {e[0]}-{e[1]:<4} {k}/{n} = {k/n:.2f}")

    print(f"\n=== ensemble robustness R(E) and within-model support ===")
    for q in queries:
        form = [1 for _, _, tot, res, _ in per_sample if res[q][0] > 0]
        sup = [res[q][0] / res[q][1] for _, _, tot, res, _ in per_sample if res[q][1]]
        R = len(form) / n
        print(f"   {q:10} R(E) = {len(form)}/{n} = {R:.2f}"
              + (f"   support {min(sup):.3f}-{max(sup):.3f} "
                 f"(mean {sum(sup)/len(sup):.3f})" if sup else ""))

    print(f"\n=== decisive interfaces ===")
    for q in queries:
        yes = set.intersection(*[e for _, e, _, r, _ in per_sample if r[q][0] > 0]) \
              if any(r[q][0] > 0 for _, _, _, r, _ in per_sample) else set()
        no = set.union(*[e for _, e, _, r, _ in per_sample if r[q][0] == 0]) \
             if any(r[q][0] == 0 for _, _, _, r, _ in per_sample) else set()
        dec = yes - no
        print(f"   {q:10} present in every formable sample and absent from every "
              f"non-formable one: {sorted(dec) if dec else 'none / all samples agree'}")


if __name__ == "__main__":
    main()
