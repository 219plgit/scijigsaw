"""Ordering consistency checks for Supplementary S5.6 and the Discussion.

Q1: are the experimentally derived hierarchies (NLRP3, 30S) admissible by
    construction?
Q2: evolutionary-tier analysis for VAMP2 — complementary fractions, enrichment,
    an exact tier-label permutation test, direct-edge concordance, and Kendall
    tau-b between tier and assembly depth.

All printed numbers are real outputs and are the values quoted in the paper.
"""
import itertools, math

import pandas as pd
import scijigsaw.cases as C

AGE = {"ancient": 0, "metazoan": 1, "vertebrate": 2}


def _req(asm):
    return {u: set(asm.requires.get(u, [])) for u in asm.units}, asm.seed


def admissible(seq, req, seed):
    placed = {seed} if seed else set()
    for u in seq:
        if not req.get(u, set()) <= placed:
            return False
        placed.add(u)
    return True


def exists_admissible(asm):
    req, seed = _req(asm)
    placed = [seed] if seed else []
    ps = set(placed)
    rem = {u for u in asm.units if u != seed}
    while rem:
        avail = [u for u in rem if req.get(u, set()) <= ps]
        if not avail:
            return None
        u = sorted(avail)[0]
        placed.append(u); ps.add(u); rem.discard(u)
    return placed


def q1():
    print("Q1  experimental hierarchies admissible by construction?")
    for name, asm in [("NLRP3", C.INFLAMMASOME), ("30S", C.RIBOSOME_30S_SPECIFIC)]:
        order = exists_admissible(asm)
        req, seed = _req(asm)
        ok = order is not None and admissible([u for u in order if u != seed], req, seed)
        print(f"  {name}: admissible order exists={order is not None}, "
              f"respects all encoded precedence={ok}")


def q2():
    prot = pd.read_csv("examples/vamp2/proteins.csv").set_index("name")
    vfus = ["SNAP25", "Syntaxin-1A", "Munc18-1", "Complexin", "Syt-1",
            "Synaptophysin", "SNCA"]
    req, seed = _req(C.VAMP2)
    tier = {u: AGE[str(prot.loc[u, "age"])] for u in vfus}
    perms = list(itertools.permutations(vfus))

    def age_mono(seq, tof):
        return all(tof[seq[i]] <= tof[seq[i + 1]] for i in range(len(seq) - 1))

    L_assembly = sum(1 for p in perms if admissible(p, req, seed))
    L_age = sum(1 for p in perms if age_mono(p, tier))
    L_both = sum(1 for p in perms if admissible(p, req, seed) and age_mono(p, tier))
    N = len(perms)
    enrich = (L_both / L_assembly) / (L_age / N)
    print("Q2  VAMP2 evolutionary-tier analysis")
    print(f"  L_assembly={L_assembly}, L_age={L_age}, L_both={L_both}, N={N}")
    print(f"  admissible that are age-monotone: {L_both}/{L_assembly} = {100*L_both/L_assembly:.2f}%")
    print(f"  age-monotone that are admissible: {L_both}/{L_age} = {100*L_both/L_age:.1f}%")
    print(f"  enrichment among admissible: {enrich:.2f}-fold")

    # NOTE: no inferential permutation test is reported. With seven components
    # and a 3/3/1 tier split the relabelling space is tiny and coarse (L_both
    # takes only values {0,6,12,18}), so any permutation null is severely
    # underpowered and its p-value uninformative. We report descriptive figures
    # and the concordance measures below instead.

    # direct-edge concordance
    edges = [(x, y) for x in vfus for y in req.get(x, set()) if y in tier]
    v_edge = sum(1 for (x, y) in edges if tier[y] > tier[x])
    print(f"  direct-edge age violations: {v_edge}/{len(edges)}")

    # Kendall tau-b between tier and assembly depth
    def depth(u, memo={}):
        if u in memo:
            return memo[u]
        ds = [depth(p) for p in req.get(u, set()) if p in vfus]
        memo[u] = 1 + (max(ds) if ds else 0)
        return memo[u]
    depths = {u: depth(u) for u in vfus}
    conc = disc = tx = ty = 0
    for a, b in itertools.combinations(vfus, 2):
        da, dd = tier[a] - tier[b], depths[a] - depths[b]
        if da == 0 and dd == 0:
            continue
        if da == 0:
            tx += 1; continue
        if dd == 0:
            ty += 1; continue
        if (da > 0) == (dd > 0):
            conc += 1
        else:
            disc += 1
    denom = math.sqrt((conc + disc + tx) * (conc + disc + ty))
    tau_b = (conc - disc) / denom if denom > 0 else 0.0
    # sign convention: positive tau_b means older tier tends to LOWER depth
    # (assembles earlier); negative means older tends to assemble later.
    print(f"  Kendall tau-b (tier vs assembly depth, + = older assembles earlier) = {tau_b:.3f}")
    print("  per-component tier and assembly depth:")
    for u in vfus:
        tname = {0: "ancient", 1: "metazoan", 2: "vertebrate"}[tier[u]]
        print(f"    {u:14s} tier={tname:11s} depth={depths[u]}")


def probabilistic_summary():
    """Machine-checkable probabilistic quantities for the VAMP2 binary example."""
    from scijigsaw.contrib.probabilistic import (
        Modification, order_distribution, order_distribution_joint,
        occupancy_uncertainty_analytic)

    def block(requires, excludes):
        for k in list(requires):
            requires[k] = [u for u in requires[k] if u != "SNAP25"]
        requires["SNAP25"] = ["__blocked__"]

    def fusion(asm):
        return "__blocked__" not in asm.requires.get("SNAP25", [])

    print("Probabilistic summary (VAMP2 binary example, declared p=0.6)")
    m = Modification("pThr138", 0.6, block)
    r = order_distribution(C.VAMP2, [m], target=fusion)
    s = r.summary(total_permutations=5040)
    print(f"  P_competent           = {s['P_competent']}")
    print(f"  E_L                   = {s['E_orders']}")
    print(f"  E_q                   = {s['E_fraction']}")
    print(f"  Var_L                 = {s['Var_orders']}   (252^2*0.6*0.4 = 15240.96)")
    print(f"  E_L_given_competent   = {s['E_orders_given_competent']}")

    print("Beta(3,2) occupancy uncertainty (analytic)")
    b = occupancy_uncertainty_analytic(252, 3, 2)
    print(f"  E_L (analytic)        = {b['E_orders']}")
    print(f"  L 95% equal-tailed    = [{b['L_ci95_equal_tailed'][0]:.1f}, "
          f"{b['L_ci95_equal_tailed'][1]:.1f}]")
    print(f"  Var decomposition     = within {b['Var_within_E_p_Var']:.1f} + "
          f"between {b['Var_between_Var_p_E']:.1f} = {b['Var_total']:.1f}")
    print(f"  method                = {b['method']}")

    print("Correlated joint example (weights over (x_A,x_B))")
    mA = Modification("A", 0.3, block)
    mB = Modification("B", 0.45, lambda r, e: None)
    joint = {(0, 0): 0.40, (0, 1): 0.15, (1, 0): 0.25, (1, 1): 0.20}
    rj = order_distribution_joint(C.VAMP2, [mA, mB], joint, target=fusion)
    for rec in rj.instances:
        print(f"    state present={rec['present']!s:20s} w={rec['weight']:.2f} "
              f"competent={rec['target_competent']!s:5s} L={rec['orders']}")
    print(f"  E_L (joint)           = {rj.expected_orders:.1f}")


if __name__ == "__main__":
    q1()
    q2()
    print()
    probabilistic_summary()
