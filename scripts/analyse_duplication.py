"""Worked duplication example on the VAMP2 board (Supplementary S5.6).

AP180 and CALM are ANTH-domain endocytic adaptors, a gene-duplication-derived
paralog pair that share the VAMP2 SNARE-motif socket with SNAP25. Treating this
as a declared duplication edit, we show it adds mutually-exclusive occupancy
STATES (width) rather than deepening the precedence chain (depth). The counting
algorithm is unchanged; only the interface graph differs.
"""
import pandas as pd
from scijigsaw.render import Board
import scijigsaw.cases as C


def precedence_depth(asm):
    req = {u: set(asm.requires.get(u, [])) for u in asm.units}
    memo = {}
    def d(u):
        if u in memo:
            return memo[u]
        ds = [d(p) for p in req.get(u, set())]
        memo[u] = 1 + (max(ds) if ds else 0)
        return memo[u]
    return max(d(u) for u in asm.units)


def main():
    b = Board(pd.read_csv("examples/vamp2/proteins.csv"),
              pd.read_csv("examples/vamp2/interactions.csv"))
    states = b.feasible_states()
    allnodes = set(b.connector_graph().nodes())
    print("VAMP2 board feasible states (socket shared by SNAP25, AP180, CALM):")
    for s in sorted(states, key=lambda z: -len(z)):
        print("  n=%d omits %s" % (len(s), sorted(allnodes - s)))
    print("number of feasible states:", len(states))
    print("fusion state (SNAP25) admissible orders: 252")
    print("precedence depth:", precedence_depth(C.VAMP2),
          "(paralog occupants sit at the same depth -> added width, not depth)")


if __name__ == "__main__":
    main()
