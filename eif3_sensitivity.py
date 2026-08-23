#!/usr/bin/env python3
"""eIF3 Comparison 2 — sensitivity of the admissible-set result to the encoding.

Runs the endorsed baseline plus several reasonable variations, and checks in each
case whether the native-MS-observed seed-anchored subcomplexes (b:i and b:g:i)
remain admissible intermediates. Robustness across variations means the result
is not an artifact of one particular encoding choice.

Run from inside your scijigsaw working copy:
    python eif3_sensitivity.py
"""
from scijigsaw.assembly import Assembly
from itertools import permutations


def analyze(name, requires, seed):
    A = Assembly(requires=requires, seed=seed)
    units = [u for u in A.units if u != seed]
    req = A.requires

    def permitted(order):
        placed = {seed}
        for u in order:
            if not req[u] <= placed:
                return False
            placed.add(u)
        return True

    perms = [p for p in permutations(units) if permitted(p)]

    def reachable(sub):
        S = sub | {seed}
        k = len(S) - 1
        return any(set(p[:k]) | {seed} == S for p in perms)

    bi = reachable({"Prt1", "Tif34"})
    bgi = reachable({"Prt1", "Tif35", "Tif34"})
    print(f"{name}")
    print(f"   seed={seed:6}  permitted={len(perms):3}/120   "
          f"b:i={'OK' if bi else 'X'}   b:g:i={'OK' if bgi else 'X'}")
    print()


# baseline endorsed encoding
base = {
    "Prt1": set(),
    "Tif34": {"Prt1"},
    "Tif35": {"Tif34"},
    "Tif32": {"Prt1"},
    "Nip1": {"Tif32"},
}

print("=== eIF3 encoding sensitivity ===\n")

# 1. baseline
analyze("1. baseline (Prt1 seed; endorsed)", base, "Prt1")

# 2. Nip1 also contacts Prt1 (Fig 5B shows c near the b module)
v2 = {
    "Prt1": set(),
    "Tif34": {"Prt1"},
    "Tif35": {"Tif34"},
    "Tif32": {"Prt1"},
    "Nip1": {"Tif32", "Prt1"},
}
analyze("2. Nip1 requires Tif32 + Prt1", v2, "Prt1")

# 3. stricter g attachment: Tif35 requires Tif34 AND Prt1
v3 = {
    "Prt1": set(),
    "Tif34": {"Prt1"},
    "Tif35": {"Tif34", "Prt1"},
    "Tif32": {"Prt1"},
    "Nip1": {"Tif32"},
}
analyze("3. Tif35 requires Tif34 + Prt1", v3, "Prt1")

# 4. alternative seed: Tif32 (a) as nucleus
v4 = {
    "Tif32": set(),
    "Prt1": {"Tif32"},
    "Tif34": {"Prt1"},
    "Tif35": {"Tif34"},
    "Nip1": {"Tif32"},
}
analyze("4. Tif32 as seed instead of Prt1", v4, "Tif32")

# 5. Tif32 as a bridge (requires Prt1, and Nip1 hangs off it) - same topology check
v5 = {
    "Prt1": set(),
    "Tif34": {"Prt1"},
    "Tif35": {"Tif34"},
    "Tif32": {"Prt1"},
    "Nip1": {"Tif32"},
}
analyze("5. control (= baseline)", v5, "Prt1")

print("Interpretation: if b:i and b:g:i stay OK across variations, the")
print("native-MS validation is robust to the encoding choice, not an artifact.")
