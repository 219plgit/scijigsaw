"""Counterfactual model comparison. Edit section 1, then: python3 example_counterfactual.py"""
from sj_prob import TypedModel, drop_relations, counterfactual_compare

# ---- 1. factual encoding (EDIT THIS BLOCK) ---------------------------------
V = ("Prt1", "Tif34", "Tif35", "Nip1", "Tif32")
C = frozenset(frozenset(e) for e in [
    ("Prt1", "Tif34"), ("Prt1", "Tif35"), ("Tif35", "Tif34"),
    ("Prt1", "Nip1"), ("Nip1", "Tif32")])
P = (("Nip1", "Tif32"),)              # declared temporal prerequisites
SEED = "Prt1"
QUERY = frozenset({"Prt1", "Tif34"})  # subcomplex whose support change we track

factual = TypedModel(V=V, C=C, P=P, X=frozenset())

# ---- 2. the counterfactual edit --------------------------------------------
# relax one declared prerequisite ("what if this temporal relation were absent?")
alt = drop_relations(factual, drop_P=[("Nip1", "Tif32")])
# other edits: drop_C=[("Prt1","Tif35")] is an in-silico interface mutation;
#              drop_X=[...] removes an exclusion (e.g. a glue-blocked site).

# ---- 3. compare -------------------------------------------------------------
r = counterfactual_compare(factual, alt, SEED, query=QUERY)
print(f"histories: {r['n_factual']} -> {r['n_alternative']}")
print(f"killed: {len(r['killed'])}  (mass {r['killed_mass_factual']:.3f})")
print(f"resurrected: {len(r['resurrected'])}  (mass {r['resurrected_mass_alt']:.3f})")
print(f"TV on shared: {r['tv_on_shared']:.3f}")
print(f"dH2: {r['delta_H2_bits']:+.3f} bits   dN_eff: {r['delta_N_eff']:+.2f}")
print(f"dPr(query {sorted(QUERY)}): {r['query_support_delta']:+.3f}")
print(f"edited relations: {r['edited_relations']}")
