"""
verify_paper.py -- check every numerical claim in the manuscript against live
computation.

One command reproduces, from the declared encodings, each figure that appears in
the text, and reports PASS or FAIL per claim with the expected and obtained
values. Intended as the release gate and as the reproducibility entry point for
a reader: nothing here is hard-coded except the CLAIMS themselves, which are
transcribed from the manuscript.

    python3 verify_paper.py            # full suite
    python3 verify_paper.py --section eif3

Requires sj_prob.py, sj_navigate.py, sj_cryoem.py in the same directory.
"""

from __future__ import annotations

import argparse
import itertools
import math
import sys
from typing import Callable, Dict, List, Tuple

from sj_prob import (TypedModel, seeded_F, merger_F, merger_support,
                     drop_relations, decision_value_relation)
from sj_navigate import Navigator
from sj_cryoem import CryoEMAnalysis, Observation

RESULTS: List[Tuple[str, str, object, object, bool]] = []


def check(section: str, claim: str, expected, got, tol: float = 1e-9) -> None:
    if isinstance(expected, float) or isinstance(got, float):
        ok = abs(float(got) - float(expected)) <= max(tol, 5e-4)
    else:
        ok = expected == got
    RESULTS.append((section, claim, expected, got, ok))


# ---------------------------------------------------------------- encodings
V_EIF3 = ("Tif32", "Prt1", "Nip1", "Tif35", "Tif34")
C_DECL = frozenset(frozenset(e) for e in [
    ("Prt1", "Tif34"), ("Tif32", "Prt1"), ("Tif35", "Tif34"), ("Tif32", "Nip1")])
P_MERGER = (("Tif34", "Tif35"),)
P_SEEDED = (("Prt1", "Tif34"), ("Prt1", "Tif32"),
            ("Tif34", "Tif35"), ("Tif32", "Nip1"))
AF3 = {frozenset({"Tif32", "Prt1"}): (1.0, 195),
       frozenset({"Tif32", "Nip1"}): (1.0, 206),
       frozenset({"Tif35", "Tif34"}): (1.0, 84),
       frozenset({"Prt1", "Tif34"}): (1.0, 73),
       frozenset({"Prt1", "Tif35"}): (1.0, 70),
       frozenset({"Prt1", "Nip1"}): (1.0, 60),
       frozenset({"Tif32", "Tif35"}): (1.0, 40),
       frozenset({"Nip1", "Tif35"}): (0.8, 92),
       frozenset({"Nip1", "Tif34"}): (0.1, 24)}
AF10 = frozenset(e for e, (r, _) in AF3.items() if r >= 1.0)
AF8 = frozenset(e for e, (r, _) in AF3.items() if r >= 0.8)
Q_GI = frozenset({"Tif35", "Tif34"})
Q_BI = frozenset({"Prt1", "Tif34"})
Q_BGI = frozenset({"Prt1", "Tif35", "Tif34"})


# ------------------------------------------------- weighted-tree machinery
def enumerate_trees(model: TypedModel):
    def rec(S):
        S = frozenset(S)
        if len(S) == 1:
            return [(next(iter(S)), [])]
        if not model.valid_intermediate(S):
            return []
        out, anchor = [], min(S)
        for r in range(len(S)):
            for combo in itertools.combinations(sorted(S - {anchor}), r):
                A = frozenset({anchor}) | frozenset(combo)
                B = S - A
                if not B or not (model.valid_intermediate(A)
                                 and model.valid_intermediate(B)):
                    continue
                if not any(model.has_contact(u, v) for u in A for v in B):
                    continue
                for ta, ma in rec(A):
                    for tb, mb in rec(B):
                        out.append((tuple(sorted([repr(ta), repr(tb)])),
                                    ma + mb + [(A, B)]))
        return out
    return rec(frozenset(model.V))


def score(model, A, B, additive: bool) -> float:
    ks = [frozenset({u, v}) for u in A for v in B
          if frozenset({u, v}) in model.C and frozenset({u, v}) in AF3]
    if not ks:
        return 0.0
    if additive:
        return sum(AF3[e][0] + math.log(1 + AF3[e][1]) for e in ks)
    r = sum(AF3[e][0] for e in ks) / len(ks)
    return r + math.log(1 + sum(AF3[e][1] for e in ks))


def weighted(model, lam, additive=False):
    rows = [(f, m, lam * sum(score(model, A, B, additive) for A, B in m))
            for f, m in enumerate_trees(model)]
    mx = max(r[2] for r in rows)
    tot = sum(math.exp(r[2] - mx) for r in rows)
    return [(f, m, math.exp(w - mx) / tot) for f, m, w in rows]


def neff(rows) -> float:
    ps = [p for _, _, p in rows if p > 0]
    return 2 ** (-sum(p * math.log2(p) for p in ps))


def wsupport(rows, E) -> float:
    tot = 0.0
    for _, merges, p in rows:
        subs = {frozenset(A) for A, _ in merges} | {frozenset(B) for _, B in merges}
        if E in subs:
            tot += p
    return tot


# ------------------------------------------------------------- the sections
def sec_eif3() -> None:
    s = "eIF3 exact"
    seeded = TypedModel(V=V_EIF3, C=C_DECL, P=P_SEEDED)
    check(s, "seeded orders (of 24)", 6, seeded_F(seeded, "Prt1")(frozenset()))
    typed = TypedModel(V=V_EIF3, C=C_DECL, P=P_MERGER)
    check(s, "typed merger trees", 14, merger_F(typed)(frozenset(V_EIF3)))
    check(s, "contacts only, no temporal claim", 14,
          merger_F(TypedModel(V=V_EIF3, C=C_DECL, P=()))(frozenset(V_EIF3)))
    for lbl, E, exp in [("support g:i", Q_GI, 5), ("support b:i", Q_BI, 5),
                        ("support b:g:i", Q_BGI, 4)]:
        g, f = merger_support(typed, E)
        check(s, f"{lbl} (of 14)", exp, g)

    s = "eIF3 AlphaFold"
    for name, C, ntrees, gi in [("10/10", AF10, 27, 12), (">=8/10", AF8, 30, 15)]:
        m = TypedModel(V=V_EIF3, C=C, P=P_MERGER)
        check(s, f"AF3 {name} trees", ntrees, merger_F(m)(frozenset(V_EIF3)))
        g, f = merger_support(m, Q_GI)
        check(s, f"AF3 {name} g:i support", gi, g)


def sec_additivity() -> None:
    s = "additivity lemma"
    m = TypedModel(V=V_EIF3, C=AF10, P=P_MERGER)
    totals = {round(sum(score(m, A, B, True) for A, B in merges), 6)
              for _, merges in enumerate_trees(m)}
    check(s, "distinct edge-additive tree scores", 1, len(totals))
    check(s, "common total score", 38.4447, next(iter(totals)), tol=5e-4)
    check(s, "sum over all contacts equals it", round(next(iter(totals)), 4),
          round(sum(AF3[e][0] + math.log(1 + AF3[e][1]) for e in m.C), 4))
    for lam in (0.0, 1.0, 4.0):
        check(s, f"edge-additive N_eff at lambda={lam}", 27.0,
              neff(weighted(m, lam, additive=True)), tol=1e-6)

    s = "lambda sweep (interface-level)"
    for name, C, exp in [("declared", C_DECL, {0: 14.0, 1: 14.0, 2: 14.0, 4: 14.0}),
                         ("AF3 10/10", AF10, {0: 27.0, 1: 24.49, 2: 20.42, 4: 15.61}),
                         ("AF3 >=8/10", AF8, {0: 30.0, 1: 27.86, 2: 24.16, 4: 18.06})]:
        mm = TypedModel(V=V_EIF3, C=C, P=P_MERGER)
        for lam, e in exp.items():
            check(s, f"{name} N_eff at lambda={lam}", e,
                  round(neff(weighted(mm, lam)), 2), tol=0.02)

    s = "robustness band"
    lo, hi = 1.0, 0.0
    for C in (C_DECL, AF10, AF8):
        mm = TypedModel(V=V_EIF3, C=C, P=P_MERGER)
        for lam in (0, 0.25, 0.5, 1, 2, 4):
            v = wsupport(weighted(mm, lam), Q_GI)
            lo, hi = min(lo, v), max(hi, v)
    check("robustness band", "min support of Tif35-Tif34", 0.357, round(lo, 3), tol=0.002)
    check("robustness band", "max support of Tif35-Tif34", 0.500, round(hi, 3), tol=0.002)


def sec_deletion() -> None:
    s = "interface deletion"
    m = TypedModel(V=V_EIF3, C=AF10, P=P_MERGER)
    expect = {("Prt1", "Tif34"): 12, ("Tif35", "Tif34"): 15, ("Nip1", "Prt1"): 17,
              ("Prt1", "Tif32"): 20, ("Nip1", "Tif32"): 20,
              ("Tif32", "Tif35"): 23, ("Prt1", "Tif35"): 27}
    for pair, exp in expect.items():
        alt = TypedModel(V=V_EIF3, C=AF10 - {frozenset(pair)}, P=P_MERGER)
        check(s, f"delete {pair[0]}--{pair[1]}", exp,
              merger_F(alt)(frozenset(V_EIF3)))


def sec_attribution() -> None:
    s = "evidence attribution"
    seeded = TypedModel(V=V_EIF3, C=C_DECL, P=P_SEEDED)
    exp = {("Tif34", "Tif35"): 1.0, ("Tif32", "Nip1"): 1.0,
           ("Prt1", "Tif34"): 0.0, ("Prt1", "Tif32"): 0.0}
    for p, e in exp.items():
        r = decision_value_relation(seeded, "Prt1", drop_P=[p], label=str(p))
        check(s, f"delta H2 for {p[0]}->{p[1]} (bits)", e, round(r["bits_removed"], 3), tol=2e-3)
    r = decision_value_relation(seeded, "Prt1", drop_P=list(P_SEEDED), label="all")
    check(s, "all prerequisites jointly (bits)", 2.0, round(r["bits_removed"], 3), tol=2e-3)


def sec_cases() -> None:
    s = "case studies"
    V = ("mTOR", "mLST8", "RAPTOR", "S6K1", "FKBP12rap")
    C = frozenset(frozenset(e) for e in [("mTOR", "mLST8"), ("mTOR", "RAPTOR"),
                                         ("mTOR", "S6K1"), ("mTOR", "FKBP12rap")])
    X = frozenset({frozenset({"S6K1", "FKBP12rap"})})
    for state in [("mTOR", "mLST8", "RAPTOR", "S6K1"),
                  ("mTOR", "mLST8", "RAPTOR", "FKBP12rap")]:
        sub = TypedModel(V=state, C=frozenset(e for e in C if e <= frozenset(state)),
                         P=(), X=X)
        check(s, f"mTOR trees {state[-1]}", 6, merger_F(sub)(frozenset(state)))
    full = TypedModel(V=V, C=C, P=(), X=X)
    check(s, "joint FRB occupancy admissible", False,
          full.valid_intermediate(frozenset({"mTOR", "S6K1", "FKBP12rap"})))

    Vp = ("13S", "b1", "b5", "b6")
    Cp = frozenset(frozenset(e) for e in [("13S", "b1"), ("13S", "b5"), ("13S", "b6")])
    imposed = TypedModel(V=Vp, C=Cp, P=(("b5", "b6"), ("b6", "b1")))
    check(s, "proteasome under imposed order", 1,
          Navigator(imposed, "13S").n_histories())
    check(s, "proteasome with precedence removed", 6,
          Navigator(TypedModel(V=Vp, C=Cp, P=()), "13S").n_histories())


def sec_cryoem() -> None:
    s = "cryo-EM multimodality"
    V = ("13S", "b1", "b5", "b6")
    C = frozenset(frozenset({"13S", x}) for x in ("b1", "b5", "b6"))
    OBS = [Observation("13S+b1", frozenset({"13S", "b1"}), "9rla", 3.16),
           Observation("13S+b5+b6", frozenset({"13S", "b5", "b6"}), "9rm0", 3.29),
           Observation("13S+b1+b5", frozenset({"13S", "b1", "b5"}), "9rm1", 4.11),
           Observation("15S", frozenset({"13S", "b1", "b5", "b6"}), "9rlz", 3.12)]
    free = CryoEMAnalysis(TypedModel(V=V, C=C, P=()), "13S", OBS)
    cov = free.minimum_cover()
    check(s, "best single mechanism explains", 3, cov["best_single"])
    check(s, "minimum coexisting mechanisms", 2, cov["n_required"])
    check(s, "mutually exclusive observation pairs", 2, len(free.conflicts()))
    loo = {r["dropped"]: r["n_required"] for r in free.leave_one_out()}
    check(s, "conclusion carried by 9rm0 (drop -> 1)", 1, loo["13S+b5+b6"])
    check(s, "dropping 9rla still needs 2", 2, loo["13S+b1"])

    imposed = CryoEMAnalysis(
        TypedModel(V=V, C=C, P=(("b5", "b6"), ("b6", "b1"))), "13S", OBS)
    imp = [r for r in imposed.admissibility() if r["verdict"] == "IMPOSSIBLE"]
    check(s, "impossible under imposed order", 2, len(imp))


SECTIONS: Dict[str, Callable[[], None]] = {
    "eif3": sec_eif3, "additivity": sec_additivity, "deletion": sec_deletion,
    "attribution": sec_attribution, "cases": sec_cases, "cryoem": sec_cryoem,
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--section", choices=list(SECTIONS) + ["all"], default="all")
    args = ap.parse_args()
    todo = SECTIONS if args.section == "all" else {args.section: SECTIONS[args.section]}
    for fn in todo.values():
        fn()

    width = max(len(c) for _, c, _, _, _ in RESULTS) + 2
    cur = None
    for sec, claim, exp, got, ok in RESULTS:
        if sec != cur:
            print(f"\n{sec.upper()}")
            print("-" * (width + 34))
            cur = sec
        mark = "PASS" if ok else "FAIL"
        print(f"  {mark}  {claim:<{width}} expected {str(exp):>10}   got {str(got):>10}")
    n = len(RESULTS)
    bad = [r for r in RESULTS if not r[4]]
    print("\n" + "=" * (width + 36))
    print(f"  {n - len(bad)}/{n} checks passed"
          + ("" if not bad else f"   -- {len(bad)} FAILED"))
    print("=" * (width + 36))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
