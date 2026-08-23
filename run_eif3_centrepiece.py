"""
run_eif3_centrepiece.py -- the §8 real-data pipeline for eIF3.

Encodings transcribed from Scientific Jigsaw Supplementary S9.10 / S8.7;
AlphaFold 3 recurrence and median interface residues from the Supplementary
AF3 table. Implements the ten-step pipeline with the §9 lambda sweep and the
§8.4 output set.

Structural score (declared, not calibrated):   s_m = a*r_m + b*log(1 + n_m)
Potential:                                     phi_m(lambda) = exp(lambda * s_m)
Tree weight:                                   W(omega) = prod_m phi_m
"""
import math
from itertools import combinations
from sj_prob import TypedModel, merger_F, merger_support

# --------------------------------------------------------------- evidence
# AF3: (recurrence out of 10, median qualifying interface residues)
AF3 = {
    frozenset({"Tif32", "Prt1"}):  (10 / 10, 195),
    frozenset({"Tif32", "Nip1"}):  (10 / 10, 206),
    frozenset({"Tif35", "Tif34"}): (10 / 10, 84),
    frozenset({"Prt1", "Tif34"}):  (10 / 10, 73),
    frozenset({"Prt1", "Tif35"}):  (10 / 10, 70),
    frozenset({"Prt1", "Nip1"}):   (10 / 10, 60),
    frozenset({"Tif32", "Tif35"}): (10 / 10, 40),
    frozenset({"Nip1", "Tif35"}):  (8 / 10, 92),
    frozenset({"Nip1", "Tif34"}):  (1 / 10, 24),
}
V = ("Tif32", "Prt1", "Nip1", "Tif35", "Tif34")
DECLARED_C = frozenset([frozenset({"Prt1", "Tif34"}), frozenset({"Tif32", "Prt1"}),
                        frozenset({"Tif35", "Tif34"}), frozenset({"Tif32", "Nip1"})])
AF3_10 = frozenset(e for e, (r, _) in AF3.items() if r >= 1.0)
AF3_8 = frozenset(e for e, (r, _) in AF3.items() if r >= 0.8)
PREREQ = (("Tif34", "Tif35"),)          # the single declared temporal claim

QUERIES = {"g:i  Tif35-Tif34": frozenset({"Tif35", "Tif34"}),
           "b:i  Prt1-Tif34": frozenset({"Prt1", "Tif34"}),
           "b:g:i Prt1-Tif35-Tif34": frozenset({"Prt1", "Tif35", "Tif34"})}

A_COEF, B_COEF = 1.0, 1.0                # declared, sensitivity-tested


# ------------------------------------------------- explicit weighted trees
def enumerate_trees(model: TypedModel):
    """All admissible merger trees on V, each as (canonical_form, [merge list])."""
    def rec(S):
        S = frozenset(S)
        if len(S) == 1:
            return [(next(iter(S)), [])]
        if not model.valid_intermediate(S):
            return []
        out, anchor = [], min(S)
        for r in range(len(S)):
            for combo in combinations(sorted(S - {anchor}), r):
                A = frozenset({anchor}) | frozenset(combo)
                B = S - A
                if not B or not (model.valid_intermediate(A)
                                 and model.valid_intermediate(B)):
                    continue
                if not any(model.has_contact(u, v) for u in A for v in B):
                    continue
                for ta, ma in rec(A):
                    for tb, mb in rec(B):
                        form = tuple(sorted([repr(ta), repr(tb)]))
                        out.append((form, ma + mb + [(A, B)]))
        return out
    return rec(frozenset(model.V))


def merge_score(model, A, B, additive=False):
    """Structural score for the merge interface (A,B).

    additive=True   s_m = SUM over crossing edges of (a*r + b*log(1+n)).
                    EDGE-ADDITIVE: provably order-invariant (see lemma below).
    additive=False  s_m = a*mean(r) + b*log(1 + SUM n)  -- an INTERFACE-LEVEL
                    score. Non-additive in the crossing edges, so it can
                    discriminate between trees whenever some merge interface is
                    crossed by more than one contact (i.e. on cyclic graphs).
    """
    ks = [frozenset({u, v}) for u in A for v in B
          if frozenset({u, v}) in model.C and frozenset({u, v}) in AF3]
    if not ks:
        return 0.0
    if additive:
        return sum(A_COEF * AF3[e][0] + B_COEF * math.log(1 + AF3[e][1]) for e in ks)
    r = sum(AF3[e][0] for e in ks) / len(ks)
    n = sum(AF3[e][1] for e in ks)
    return A_COEF * r + B_COEF * math.log(1 + n)


def weighted_trees(model, lam, additive=False):
    rows = []
    for form, merges in enumerate_trees(model):
        logw = lam * sum(merge_score(model, A, B, additive) for A, B in merges)
        rows.append((form, merges, logw))
    mx = max(r[2] for r in rows)
    tot = sum(math.exp(r[2] - mx) for r in rows)
    return [(f, m, math.exp(lw - mx) / tot) for f, m, lw in rows]


def describe(merges):
    """Readable tree: the non-trivial subcomplexes it builds, largest last."""
    parts = sorted({frozenset(A) for A, _ in merges} | {frozenset(B) for _, B in merges},
                   key=lambda s: (len(s), sorted(s)))
    return " ; ".join("+".join(sorted(p)) for p in parts if 2 <= len(p) < len(V))


def h2_neff(probs):
    h = -sum(p * math.log2(p) for p in probs if p > 0)
    return h, 2 ** h


def support_weighted(rows, E):
    """Probability that E occurs as a subcomplex, under the weighted measure."""
    tot = 0.0
    for _, merges, p in rows:
        subs = {frozenset(A) for A, _ in merges} | {frozenset(B) for _, B in merges}
        subs.add(frozenset(V))
        if E in subs or len(E) == 1:
            tot += p
    return tot


# ------------------------------------------------------------------ report
def analyse(name, C, lams=(0.0, 0.25, 0.5, 1.0, 2.0, 4.0), additive=False):
    model = TypedModel(V=V, C=C, P=PREREQ)
    n = merger_F(model)(frozenset(V))
    cyc = len(C) - (len(V) - 1)
    print("=" * 76)
    print(f"{name}: {len(C)} contacts, {n} admissible trees, "
          f"cyclomatic number {cyc} ({'TREE topology' if cyc == 0 else 'contains cycles'})")
    print("=" * 76)

    print("  exact (unit-weight) supports:")
    for lbl, E in QUERIES.items():
        g, f = merger_support(model, E)
        print(f"    {lbl:<24} {g}/{f} = {g/f:.3f}")

    print(f"\n  lambda sweep  (N_eff, top-3 and top-5 cumulative mass, weighted supports)")
    hdr = f"    {'lam':>5} {'N_eff':>7} {'top3':>6} {'top5':>6} " + \
          " ".join(f"{lbl.split()[0]:>7}" for lbl in QUERIES)
    print(hdr)
    ranks = {}
    for lam in lams:
        rows = weighted_trees(model, lam, additive)
        ps = sorted((p for _, _, p in rows), reverse=True)
        h, neff = h2_neff([p for _, _, p in rows])
        sup = [support_weighted(rows, E) for E in QUERIES.values()]
        print(f"    {lam:>5} {neff:>7.2f} {sum(ps[:3]):>6.3f} {sum(ps[:5]):>6.3f} " +
              " ".join(f"{s:>7.3f}" for s in sup))
        order = [f for f, _, _ in sorted(rows, key=lambda r: -r[2])]
        for i, f in enumerate(order):
            ranks.setdefault(f, []).append(i + 1)

    rows1 = sorted(weighted_trees(model, 1.0, additive), key=lambda r: -r[2])
    print(f"\n  top 5 weighted trees at lambda=1 (rank range across the sweep):")
    for f, merges, p in rows1[:5]:
        rr = ranks[f]
        stable = "stable" if max(rr) - min(rr) == 0 else f"ranks {min(rr)}-{max(rr)}"
        print(f"    P={p:.4f}  [{stable}]  builds: {describe(merges) or '(no intermediates)'}")
    return model


# ------------------------------------------------------------------- main
print("LEMMA CHECK: edge-additive vs interface-level potentials (AF3 10/10)")
_m = TypedModel(V=V, C=AF3_10, P=PREREQ)
for _add in (True, False):
    _row = [f"{h2_neff([p for _,_,p in weighted_trees(_m, l, _add)])[1]:.2f}" for l in (0,1,4)]
    _sup = [f"{support_weighted(weighted_trees(_m, l, _add), QUERIES['g:i  Tif35-Tif34']):.3f}" for l in (0,1,4)]
    print(f"  {'edge-additive ' if _add else 'interface-level'}: N_eff(lam=0,1,4) = {_row}   P(g:i) = {_sup}")
print()
print("STEP 1-8: weighted merger-tree inference on the real eIF3 encodings\n")
m_dec = analyse("DECLARED literature encoding", DECLARED_C)
print()
m_af3 = analyse("AF3 recurrent contacts (10/10)", AF3_10)
print()
m_af8 = analyse("AF3 contacts (>=8/10)", AF3_8)

print()
print("=" * 76)
print("STEP 9: robustness across structural hypotheses; STEP 10: hypothesis query")
print("=" * 76)
# hypothesis H: the off-scaffold g:i dimer forms as an independent module
H = QUERIES["g:i  Tif35-Tif34"]
for nm, mdl in [("declared", m_dec), ("AF3 10/10", m_af3), ("AF3 >=8/10", m_af8)]:
    vals = [support_weighted(weighted_trees(mdl, l), H) for l in (0.0, 1.0, 4.0)]
    print(f"  P(H: Tif35-Tif34 forms as a module) under {nm:<12} "
          f"lam=0:{vals[0]:.3f}  lam=1:{vals[1]:.3f}  lam=4:{vals[2]:.3f}")

print("\n  counterfactual: in-silico deletion of the Prt1-Tif35 interface (AF3 10/10)")
alt = TypedModel(V=V, C=AF3_10 - {frozenset({"Prt1", "Tif35"})}, P=PREREQ)
f0 = merger_F(m_af3)(frozenset(V))
f1 = merger_F(alt)(frozenset(V))
p0 = support_weighted(weighted_trees(m_af3, 1.0), H)
p1 = support_weighted(weighted_trees(alt, 1.0), H)
h0 = h2_neff([p for _, _, p in weighted_trees(m_af3, 1.0)])
h1 = h2_neff([p for _, _, p in weighted_trees(alt, 1.0)])
print(f"    trees {f0} -> {f1}   P(H) {p0:.3f} -> {p1:.3f} "
      f"(dPr {p1-p0:+.3f})   N_eff {h0[1]:.2f} -> {h1[1]:.2f} "
      f"(dH2 {h1[0]-h0[0]:+.3f} bits)")
