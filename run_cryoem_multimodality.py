"""Cryo-EM multimodality: Mark et al. 2026 proteasome precursor structures.

Six deposited yeast PC intermediates (EMD-54029/PDB 9rl3 13S; EMD-54032/9rla
13S+b1; EMD-54047/9rm0 13S+b5+b6; EMD-54048/9rm1 13S+b1+b5; EMD-54045/9rlt
dimerised 13S-13S+b5; EMD-54046/9rlz 15S). Each resolved intermediate enters as
a QUERY, not as a contact or a precedence.

Question: can any single assembly mechanism explain all observed intermediates?
"""
from itertools import permutations
from sj_prob import TypedModel, seeded_F
from sj_navigate import Navigator

SEED = "13S"
V = (SEED, "b1", "b5", "b6")
C = frozenset(frozenset(e) for e in [(SEED,"b1"), (SEED,"b5"), (SEED,"b6")])

OBSERVED = {                      # exact intermediates resolved by cryo-EM
 "13S+b1        (9rla, 3.16 A)": {SEED,"b1"},
 "13S+b5+b6     (9rm0, 3.29 A)": {SEED,"b5","b6"},
 "13S+b1+b5     (9rm1, 4.11 A)": {SEED,"b1","b5"},
 "15S           (9rlz, 3.12 A)": {SEED,"b1","b5","b6"},
}

def report(model, title):
    nav = Navigator(model, SEED)
    print(f"\n{title}: {nav.n_histories()} admissible mechanisms")
    for lbl, E in OBSERVED.items():
        r = nav.ask(nav.forms_exact(E))
        n = round(r["prob"]*nav.n_histories())
        tag = "IMPOSSIBLE" if r["prob"]==0 else f"{n}/{nav.n_histories()} = {r['prob']:.3f}"
        print(f"   {lbl:<30} {tag}" + (f"   [{r['cause']}]" if r["cause"] else ""))
    return nav

print("="*74)
print("1. UNDER THE PREVIOUSLY ASSUMED ORDER  b5 -> b6 -> b1")
print("="*74)
imposed = TypedModel(V=V, C=C, P=(("b5","b6"),("b6","b1")))
report(imposed, "imposed precedence")

print()
print("="*74)
print("2. WITH THE UNSUPPORTED PRECEDENCE REMOVED (counterfactual)")
print("="*74)
free = TypedModel(V=V, C=C, P=())
nav = report(free, "no imposed precedence")

print()
print("="*74)
print("3. MULTIMODALITY TEST: does ONE mechanism explain ALL observations?")
print("="*74)
hist = [h for h,_ in nav.masses]
covered = {}
for h in hist:
    got = set()
    for lbl, E in OBSERVED.items():
        cur = {SEED}
        if cur == E: got.add(lbl)
        for p in h:
            cur.add(p)
            if cur == E: got.add(lbl)
    covered[h] = got
    print(f"   {SEED}->{'->'.join(h):<14} explains {len(got)}/{len(OBSERVED)}: "
          + ", ".join(sorted(x.split()[0] for x in got)))

best = max(len(v) for v in covered.values())
print(f"\n   best single mechanism explains {best} of {len(OBSERVED)} observed intermediates")

# minimum number of mechanisms needed to cover every observation (exact set cover)
import itertools
need = None
for k in range(1, len(hist)+1):
    for combo in itertools.combinations(hist, k):
        if set().union(*(covered[c] for c in combo)) == set(OBSERVED):
            need = (k, combo); break
    if need: break
k, combo = need
print(f"   MINIMUM number of distinct mechanisms required: {k}")
for c in combo:
    print(f"      {SEED}->{'->'.join(c)}   covers "
          + ", ".join(sorted(x.split()[0] for x in covered[c])))
print(f"\n   => the cryo-EM intermediates are MULTIMODAL: no single assembly")
print(f"      mechanism accounts for all of them; at least {k} coexisting routes")
print(f"      are required. This is an exact statement about the declared model,")
print(f"      derived from structures alone, with no kinetic assumption.")
