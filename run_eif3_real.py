"""CENTREPIECE: weighted merger-tree inference on the REAL eIF3 encoding.

Encoding transcribed from Scientific Jigsaw Supplementary S9.10 (seeded) and
S8.7 (typed merger trees). Interface-size potentials from Supplementary
Table (AF3 recurrent contacts, median qualifying interface residues).
"""
import math
from sj_prob import (TypedModel, seeded_F, merger_F, merger_support,
                     merger_Z, merger_branch_probs, merger_weighted_support,
                     decision_value_relation, decision_table, entropy_bits)

V = ("Tif32", "Prt1", "Nip1", "Tif35", "Tif34")
C = frozenset(frozenset(e) for e in [
    ("Prt1", "Tif34"), ("Tif32", "Prt1"), ("Tif35", "Tif34"), ("Tif32", "Nip1")])
TYPED = TypedModel(V=V, C=C, P=(("Tif34", "Tif35"),))
QUERIES = {"g:i Tif35-Tif34": {"Tif35", "Tif34"},
           "b:i Prt1-Tif34": {"Prt1", "Tif34"},
           "b:g:i Prt1-Tif35-Tif34": {"Prt1", "Tif35", "Tif34"}}

# AF3 median qualifying interface residues (Supplementary AF3 table)
IFACE = {frozenset({"Tif32", "Prt1"}): 195, frozenset({"Tif32", "Nip1"}): 206,
         frozenset({"Tif35", "Tif34"}): 84, frozenset({"Prt1", "Tif34"}): 73}

def phi(A, B, lam=1.0):
    """Interface-size potential: sum of qualifying residues across the merge
    interface, raised to lambda. Strictly positive (floor 1.0)."""
    s = sum(IFACE.get(frozenset({u, v}), 0)
            for u in A for v in B if TYPED.has_contact(u, v))
    return max(float(s), 1.0) ** lam

print("=" * 72)
print("1. EXACT LAYER (reproduces the published numbers)")
print("=" * 72)
SEEDED = TypedModel(V=V, C=C, P=(("Prt1","Tif34"), ("Prt1","Tif32"),
                                 ("Tif34","Tif35"), ("Tif32","Nip1")))
print(f"  seeded orders (full precedence): {seeded_F(SEEDED,'Prt1')(frozenset())}  [paper: 6 of 24]")
print(f"  typed merger trees: {merger_F(TYPED)(frozenset(V))}  [paper: 14]")
for lbl, E in QUERIES.items():
    g, f = merger_support(TYPED, frozenset(E))
    print(f"  support {lbl:<24} {g}/{f} = {g/f:.3f}")

print()
print("=" * 72)
print("2. WEIGHTED LAYER (interface-size potentials; space unchanged)")
print("=" * 72)
for lam in [0.0, 0.5, 1.0, 2.0]:
    Zf = merger_Z(TYPED, lambda A, B, l=lam: phi(A, B, l))
    z = Zf(frozenset(V))
    row = []
    for lbl, E in QUERIES.items():
        ze, zz = merger_weighted_support(TYPED, frozenset(E),
                                         lambda A, B, l=lam: phi(A, B, l))
        row.append(f"{lbl.split()[0]}={ze/zz:.3f}")
    print(f"  lambda={lam:<4} Z={z:10.3e}   " + "  ".join(row))
print("  (lambda=0 is the uniform case: supports return to 0.357/0.357/0.286)")

print()
print("=" * 72)
print("3. PROBABILISTIC ASSEMBLY TREE (root splits, lambda=1)")
print("=" * 72)
bp = merger_branch_probs(TYPED, frozenset(V), lambda A, B: phi(A, B, 1.0))
for A, B, p in sorted(bp, key=lambda t: -t[2]):
    print(f"  {sorted(A)} | {sorted(B)}   P = {p:.4f}")
h = -sum(p * math.log2(p) for _, _, p in bp if p > 0)
print(f"  root branch entropy: {h:.3f} bits over {len(bp)} splits")

print()
print("=" * 72)
print("4. EVIDENCE ATTRIBUTION (Delta H2, seeded mode)")
print("=" * 72)
seeded = SEEDED
rows = [decision_value_relation(seeded, "Prt1", drop_P=[p],
                                label=f"prerequisite {p[0]}->{p[1]}")
        for p in seeded.P]
rows.append(decision_value_relation(seeded, "Prt1", drop_P=list(seeded.P),
                                    label="all declared prerequisites"))
print(decision_table(rows))


print()
print("=" * 72)
print("5. WHY INTERFACE WEIGHTING IS INERT HERE: TREE-TOPOLOGY LEMMA")
print("=" * 72)
print("  The declared eIF3 contact graph is a TREE (5 nodes, 4 edges: "
      "Nip1-Tif32-Prt1-Tif34-Tif35).")
print("  Every binary merger tree on it cuts each contact edge exactly once, so a")
print("  potential depending only on the cut interface gives every tree the same")
print("  weight: 195*206*84*73 = %d, times 14 trees = Z." % (195*206*84*73))
print("  Interface weighting is therefore ORDER-INVARIANT on tree-topology contact")
print("  graphs -- the merger-mode analogue of the Plackett-Luce lemma. Weighting")
print("  discriminates only when the contact graph contains a cycle, as the denser")
print("  AlphaFold3 graphs do (see run_eif3_af3.py).")
