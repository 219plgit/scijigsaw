"""Assembly-order space and molecular crowding (Supplementary S5.9).

The admissible fraction q = L / n! is the probability that a uniformly random
arrival order is assembly-competent under the declared constraints -- a
combinatorial proxy for robustness to the out-of-order encounters that
macromolecular crowding produces. This is NOT a kinetic or thermodynamic model
of crowding; it reports which orders remain viable, not rates or yields.
"""
import math
import scijigsaw.cases as C


def q_random_order(asm):
    """Exact q = L / n! for an assembly (n = non-seed component count)."""
    n = len([u for u in asm.units if u != asm.seed])
    L = asm.n_orders_permitted()
    return L, math.factorial(n), (L / math.factorial(n))


def main():
    print("q = P(uniformly random arrival order is admissible) = L / n!")
    for name, asm in [("VAMP2 (shallow/wide)", C.VAMP2),
                      ("NLRP3 (deep/narrow)", C.INFLAMMASOME)]:
        L, N, q = q_random_order(asm)
        print(f"  {name}: L={L}, n!={N}, q={q:.3g}")
    print("  -> wide, shallow boards tolerate stochastic out-of-order encounter;")
    print("     deep, narrow boards are fragile and need chaperoning/retry.")
    print("  (combinatorial robustness proxy, not a kinetic crowding model)")


if __name__ == "__main__":
    main()
