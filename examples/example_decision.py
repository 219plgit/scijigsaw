"""Rank evidence sources by bits of entropy removed. Edit the model, then run."""
from sj_prob import (TypedModel, entropy_bits, decision_table,
                     decision_value_relation, decision_value_model_choice,
                     decision_value_weights)

# ---- 1. your encoding ------------------------------------------------------
V = ("Prt1", "Tif34", "Tif35", "Nip1", "Tif32")
C = frozenset(frozenset(e) for e in [
    ("Prt1", "Tif34"), ("Prt1", "Tif35"), ("Tif35", "Tif34"),
    ("Prt1", "Nip1"), ("Nip1", "Tif32"),
])
P = (("Nip1", "Tif32"),)          # declared temporal prerequisites
X = frozenset()                   # exclusions
SEED = "Prt1"

M = TypedModel(V=V, C=C, P=P, X=X)
print("baseline entropy:", round(entropy_bits(M, SEED), 3), "bits\n")

rows = []

# ---- 2. hard relations: one row per declared prerequisite/exclusion --------
for (p, q) in P:
    rows.append(decision_value_relation(M, SEED, drop_P=[(p, q)],
                                        label=f"prerequisite {p}->{q}"))
if len(P) > 1:
    rows.append(decision_value_relation(M, SEED, drop_P=list(P),
                                        label="all temporal prerequisites"))
for e in X:
    u, v = sorted(e)
    rows.append(decision_value_relation(M, SEED, drop_X=[e],
                                        label=f"exclusion {u}--{v}"))

# ---- 3. structural ensemble: one TypedModel per AlphaFold/Boltz variant ----
alt = TypedModel(V=V, C=C - {frozenset({"Tif35", "Tif34"})}, P=P, X=X)
rows.append(decision_value_model_choice([(0.8, M), (0.2, alt)], SEED,
                                        label="AlphaFold contact variant"))

# ---- 4. context selection: needs a DISTINCT model per Tahoe/LINCS state ----
state_on = M
state_off = TypedModel(V=V, C=C, P=P + (("Tif34", "Tif35"),), X=X)
rows.append(decision_value_model_choice([(0.5, state_on), (0.5, state_off)], SEED,
                                        label="Tahoe/LINCS state selection"))

# ---- 5. soft abundance weights --------------------------------------------
abundance = {"Tif34": 1.0, "Tif35": 0.6, "Nip1": 1.4, "Tif32": 0.9}
rows.append(decision_value_weights(M, SEED, abundance, label="expression weights"))

print(decision_table(rows))
