"""Template: plug your own encoding into sj_prob.py."""
from sj_prob import (TypedModel, seeded_F, merger_F, merger_support,
                     calibrated_tree, prefix_contains, prob_event, n_eff,
                     plackett_luce, support_interval, rank_histories,
                     model_distribution)

# 1. Encode the typed model -------------------------------------------------
V = ("Prt1", "Tif34", "Tif35", "Nip1", "Tif32")
C = frozenset(frozenset(e) for e in [
    ("Prt1", "Tif34"), ("Prt1", "Tif35"), ("Tif35", "Tif34"),
    ("Prt1", "Nip1"), ("Nip1", "Tif32"),
])
P = ()                      # directed prerequisites (p, q): p precedes q
X = frozenset()             # exclusions, e.g. frozenset({frozenset({"A","B"})})
eta = {"P:Tif34": "native-MS, Politis 2015"}
M = TypedModel(V=V, C=C, P=P, X=X, eta=eta)

# 2. Exact enumeration ------------------------------------------------------
print("merger trees:", merger_F(M)(frozenset(V)))
print("seeded orders:", seeded_F(M, "Prt1")(frozenset()))

# 3. Exact subcomplex support ----------------------------------------------
for E in [{"Tif35", "Tif34"}, {"Prt1", "Tif34"}, {"Prt1", "Tif35", "Tif34"}]:
    g, f = merger_support(M, frozenset(E))
    print(f"support({sorted(E)}) = {g}/{f} = {g/f:.3f}" if f else "no trees")

# 4. Weighted ranking from an abundance vector per tissue -------------------
tissues = [
    {"Tif34": 1.0, "Tif35": 0.6, "Nip1": 1.4, "Tif32": 0.9},
    {"Tif34": 1.2, "Tif35": 0.5, "Nip1": 1.1, "Tif32": 1.3},
]
res = support_interval(M, "Prt1", frozenset({"Prt1", "Tif34"}), tissues)
print("support interval:", {k: round(v, 3) for k, v in res.items()})

for r in rank_histories(M, "Prt1", tissues, top=3):
    print("  ", "->".join(r["history"]), round(r["mass_mean"], 4),
          "stable" if r["rank_stable"] else "unstable")
