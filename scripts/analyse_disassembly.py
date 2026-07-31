"""Disassembly-order counting (Supplementary S5.7).

Demonstrates that, when disassembly is the microscopic reverse of assembly, the
admissible disassembly orders are the linear extensions of the DUAL of the
assembly precedence poset, and that a poset and its dual have equally many
linear extensions. Active disassembly (e.g. NSF/alpha-SNAP-driven SNARE
separation) is not a reversal and would require a separately declared
disassembly precedence relation, which the same enumerator would then process.
"""
import itertools
import random


def linext_count(n, precedes):
    """Count linear extensions given a set of (a, b) meaning a must precede b."""
    req = {i: set() for i in range(n)}
    for a, b in precedes:
        req[b].add(a)

    def ok(seq):
        placed = set()
        for u in seq:
            if not req[u] <= placed:
                return False
            placed.add(u)
        return True

    return sum(1 for p in itertools.permutations(range(n)) if ok(p))


def check_dual_equality(trials=200, seed=0):
    rng = random.Random(seed)
    mismatches = 0
    for _ in range(trials):
        n = rng.randint(3, 6)
        edges = {(a, b) for a in range(n) for b in range(n)
                 if a < b and rng.random() < 0.4}
        fwd = linext_count(n, edges)
        dual = linext_count(n, {(b, a) for (a, b) in edges})
        if fwd != dual:
            mismatches += 1
    return trials, mismatches


def main():
    trials, mismatches = check_dual_equality()
    print("Disassembly = reverse of assembly (passive, reversible case):")
    print("  assembly orders = linear extensions of the precedence poset P")
    print("  disassembly orders = linear extensions of the dual poset P*")
    print(f"  tested {trials} random posets; count(P) != count(P*) in "
          f"{mismatches} cases")
    print("  -> equal by bijection (reverse each extension); no algorithm change.")
    print()
    print("Active disassembly (e.g. NSF/alpha-SNAP SNARE separation) is not a")
    print("reversal; it needs a separately declared disassembly precedence")
    print("relation, which the same enumerator then counts.")


if __name__ == "__main__":
    main()


def counterfactual_demo():
    """Counterfactual histories: factorization and probabilistic treatment (S5.8)."""
    from scijigsaw.contrib.probabilistic import (
        order_distribution_counterfactual, linext_factorizes_disjoint)
    from scijigsaw.assembly import Assembly

    # illustrative declared histories (abstract board; not an empirical VAMP2 claim)
    H1 = Assembly(requires={"B": ["A"], "C": ["A"], "D": ["B", "C"]}, seed="s")
    H0 = Assembly(requires={"C": ["A"], "D": ["B", "C"]}, seed="s")
    print("Counterfactual histories (declared alternative boards):")
    print(f"  H1 (event occurred): L = {H1.n_orders_permitted()}")
    print(f"  H0 (event absent):   L = {H0.n_orders_permitted()}")
    print("  probabilistic counterfactual (declared scenario prior q on H1):")
    for q in (0.0, 0.5, 1.0):
        r = order_distribution_counterfactual({"H1": H1, "H0": H0},
                                              {"H1": q, "H0": 1 - q})
        print(f"    q={q}: E[L]={r.expected_orders:.2f}  Var={r.variance():.3f}")
    print("  (q is a declared scenario weight, not an inferred evolutionary probability)")
    print("  factorization over INDEPENDENT events: "
          f"L(P u Q)=C(5,2)*L(P)*L(Q)={linext_factorizes_disjoint(2,1,3,1)}")


if __name__ == "__main__":
    main()
    print()
    counterfactual_demo()
