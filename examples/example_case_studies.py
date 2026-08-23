"""The v9 case studies as navigator sessions. Encodings for mTOR and the
proteasome late stage are rebuilt from the manuscript's own description;
eIF3/SNARE contact graphs are placeholders until the release encodings load."""
from sj_prob import TypedModel, merger_F, merger_support
from sj_navigate import Navigator, before

print("=" * 70)
print("CASE 1: rapamycin-mTOR (molecular glue; shared FRB site exclusion)")
print("=" * 70)
V = ("mTOR", "mLST8", "RAPTOR", "S6K1", "FKBP12rap")
C = frozenset(frozenset(e) for e in [
    ("mTOR", "mLST8"), ("mTOR", "RAPTOR"),
    ("mTOR", "S6K1"), ("mTOR", "FKBP12rap")])       # both engage the FRB site
X = frozenset({frozenset({"S6K1", "FKBP12rap"})})   # structure-supported exclusion
mtor = TypedModel(V=V, C=C, P=(), X=X)

for state, name in [(("mTOR","mLST8","RAPTOR","S6K1"), "substrate state"),
                    (("mTOR","mLST8","RAPTOR","FKBP12rap"), "drug-bound state")]:
    sub = TypedModel(V=state, C=frozenset(e for e in C if e <= frozenset(state)),
                     P=(), X=X)
    print(f"  {name}: merger trees = {merger_F(sub)(frozenset(state))}   "
          f"(paper: exactly six)")
joint = frozenset({"mTOR", "S6K1", "FKBP12rap"})
print(f"  joint FRB occupancy admissible intermediate? "
      f"{mtor.valid_intermediate(joint)}   "
      f"cause: {mtor.blocking_reason(joint)}")

print()
print("=" * 70)
print("CASE 2: proteasome late stage (assumption failure -> do() repair)")
print("=" * 70)
V = ("13S", "b1", "b5", "b6")
C = frozenset(frozenset(e) for e in [("13S","b1"), ("13S","b5"), ("13S","b6")])
imposed = TypedModel(V=V, C=C, P=(("b5","b6"), ("b6","b1")))  # unsupported order
nav = Navigator(imposed, "13S")
print(f"  under imposed b5->b6->b1: histories = {nav.n_histories()} of 6")
r = nav.ask(before("b1", "b5"))
print(f"  ask before(b1,b5): p={r['prob']:.3f}  {r['verdict']}"
      + (f"  [{r['cause']}]" if r['cause'] else ""))
repaired = nav.do(drop_P=[("b5","b6"), ("b6","b1")])
cmp_ = nav.compare_with(repaired)
print(f"  do(remove precedence): {cmp_['n_factual']} -> {cmp_['n_alternative']} "
      f"histories (paper: all six restored), dH2 {cmp_['delta_H2_bits']:+.3f} bits, "
      f"cause {cmp_['edited_relations']['P_removed']}")

print()
print("=" * 70)
print("CASE 3: eIF3 representation failure (placeholder contact graph)")
print("=" * 70)
V = ("Prt1", "Tif34", "Tif35", "Nip1", "Tif32")
C = frozenset(frozenset(e) for e in [
    ("Prt1","Tif34"), ("Prt1","Tif35"), ("Tif35","Tif34"),
    ("Prt1","Nip1"), ("Nip1","Tif32")])
eif3 = TypedModel(V=V, C=C, P=())
nav = Navigator(eif3, "Prt1")
r = nav.ask(nav.forms_exact({"Tif35", "Tif34"}))   # session-aware: seed auto-included
print(f"  seeded forms_exact(Tif35-Tif34): p={r['prob']:.3f}  [{r['cause']}]")
g, f = merger_support(eif3, frozenset({"Tif35", "Tif34"}))
print(f"  merger-tree support of the same dimer: {g}/{f} = {g/f:.3f}")
print("  -> same query, two representations: the paper's Table 1 in two calls.")
print("  (swap in the release encoding to reconcile against 14 trees, 5/14)")
