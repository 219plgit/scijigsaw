"""Run the weighted merger-tree layer, and verify each tree is counted exactly once."""
from itertools import combinations
from sj_prob import (TypedModel, merger_F, merger_Z, merger_branch_probs,
                     merger_support, merger_weighted_support)

# ---- 1. encoding (EDIT THIS BLOCK) ------------------------------------------
V = ("P", "T34", "T35", "N", "T32")
C = frozenset(frozenset(e) for e in [
    ("P", "T34"), ("P", "T35"), ("T35", "T34"), ("P", "N"), ("N", "T32")])
M = TypedModel(V=V, C=C, P=(), X=frozenset())
E = frozenset({"T35", "T34"})          # queried subcomplex

# ---- 2. exact counts and support --------------------------------------------
print("trees:", merger_F(M)(frozenset(V)))
g, f = merger_support(M, E)
print(f"support({sorted(E)}) = {g}/{f}")

# ---- 3. weighted version -----------------------------------------------------
def phi(A, B):                          # evidence weight per merge (edit freely)
    return float(sum(1 for u in A for v in B if M.has_contact(u, v)))

ze, z = merger_weighted_support(M, E, phi)
print(f"weighted support = {ze/z:.3f}   (Z = {z:.1f})")
for A, B, p in sorted(merger_branch_probs(M, frozenset(V), phi), key=lambda t: -t[2]):
    print(f"  root split {sorted(A)} | {sorted(B)}  P = {p:.3f}")

# ---- 4. uniqueness check: build every tree explicitly and count distinct -----
def all_trees(S):
    """Every admissible tree on S as a canonical hashable object."""
    S = frozenset(S)
    if len(S) == 1:
        return [next(iter(S))]
    if not M.valid_intermediate(S):
        return []
    out = []
    anchor = min(S)
    for r in range(len(S)):
        for combo in combinations(sorted(S - {anchor}), r):
            A = frozenset({anchor}) | frozenset(combo)
            B = S - A
            if not B or not (M.valid_intermediate(A) and M.valid_intermediate(B)):
                continue
            if not any(M.has_contact(u, v) for u in A for v in B):
                continue
            for ta in all_trees(A):
                for tb in all_trees(B):
                    out.append(frozenset({("L", repr(ta)), ("R", repr(tb))})
                               if False else (min(repr(ta), repr(tb)),
                                              max(repr(ta), repr(tb))))
    return out

trees = all_trees(V)
print("\nuniqueness check:")
print("  explicit trees generated:", len(trees))
print("  distinct trees:          ", len(set(trees)))
print("  DP count:                ", merger_F(M)(frozenset(V)))
assert len(trees) == len(set(trees)) == merger_F(M)(frozenset(V)), "DOUBLE COUNTING"
print("  each admissible tree counted exactly once")
