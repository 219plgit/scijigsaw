#!/usr/bin/env python3
"""Exact merger-tree calculation for the rapamycin / mTOR FRB worked example."""
from functools import lru_cache
from itertools import combinations

CORE = ("mTOR", "mLST8", "RAPTOR")
CONTROL = CORE + ("S6K1",)
RAPA = CORE + ("FKBP12:rapamycin",)

def E(a, b):
    return frozenset((a, b))

CONTROL_EDGES = {
    E("mTOR", "mLST8"),
    E("mTOR", "RAPTOR"),
    E("mTOR", "S6K1"),
}
RAPA_EDGES = {
    E("mTOR", "mLST8"),
    E("mTOR", "RAPTOR"),
    E("mTOR", "FKBP12:rapamycin"),
}

def count(nodes, edges):
    edges = set(edges)
    @lru_cache(None)
    def F(state):
        S = frozenset(state)
        if len(S) == 1:
            return 1
        first = min(S)
        rest = sorted(S - {first})
        total = 0
        for r in range(len(rest) + 1):
            for pick in combinations(rest, r):
                A = frozenset((first,) + pick)
                if A == S:
                    continue
                B = S - A
                if any(E(a, b) in edges for a in A for b in B):
                    total += F(tuple(sorted(A))) * F(tuple(sorted(B)))
        return total
    return F(tuple(sorted(nodes)))

print("substrate-selected state:", count(CONTROL, CONTROL_EDGES), "merger trees")
print("rapamycin-bound state:   ", count(RAPA, RAPA_EDGES), "merger trees")
print("co-occupied FRB state:    0 (excluded)")
print("reason: S6K1 and FKBP12:rapamycin are declared alternative occupants")
print("of the same mTOR FRB pocket.")
