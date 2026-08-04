"""Typed-edge Assembly: contacts vs prerequisites, with legacy compatibility.

RELATION TYPES
  contact(u,v)      undirected; u and v can form an interface. Says nothing
                    about order.
  prerequisite u<v  directed; the state containing u must exist before v is
                    introduced.
  exclusion u _|_ v symmetric; u and v cannot coexist in the selected state.

LEGACY COMPATIBILITY
  Assembly(requires={...}) keeps its current meaning: every listed partner is a
  PREREQUISITE (and, implicitly, also a contact). All existing counts are
  therefore unchanged.

NEW FORM
  Assembly(contacts=[(u,v),...], prerequisites=[(u,v),...], ...)
  where (u,v) in prerequisites means u must precede v.
"""
from __future__ import annotations
import math
from typing import Dict, Iterable, Set, Tuple


# provenance -> suggested relation kind (SUGGESTION ONLY, never applied silently)
PROVENANCE_SUGGESTION = {
    "structure":            "contact",
    "pisa":                 "contact",
    "predicted_structure":  "contact",
    "crosslinking":         "contact",
    "pulldown":             "contact",
    "native_ms":            "observed_subcomplex",
    "time_resolved":        "prerequisite",
    "nucleated_cascade":    "prerequisite",
    "curated":              None,      # encoder must declare
}


class TypedAssembly:
    def __init__(self, requires=None, contacts=(), prerequisites=(),
                 excludes=(), seed=None, provenance=None,
                 observed_subcomplexes=()):
        self.seed = seed
        self.provenance = dict(provenance or {})
        self.observed = [set(e) for e in observed_subcomplexes]

        # --- legacy: requires => prerequisite AND contact -------------------
        self.contacts: Set[frozenset] = set()
        self.pre: Dict[str, Set[str]] = {}
        units: Set[str] = set()

        if requires:
            for u, ps in requires.items():
                units.add(u); self.pre.setdefault(u, set())
                for p in ps:
                    units.add(p)
                    self.pre[u].add(p)
                    self.contacts.add(frozenset((u, p)))
        for u, v in contacts:
            units.update((u, v)); self.contacts.add(frozenset((u, v)))
        for u, v in prerequisites:                    # u must precede v
            units.update((u, v))
            self.pre.setdefault(v, set()).add(u)
            self.contacts.add(frozenset((u, v)))      # a prerequisite implies contact

        if seed: units.add(seed)
        self.units = sorted(units)
        for u in self.units: self.pre.setdefault(u, set())
        self.excludes = [(set(a), set(b)) for a, b in excludes]
        self._warnings = self._consistency_check()

    # -------------------------------------------------------------- checks
    def _consistency_check(self):
        """Flag prerequisites whose declared provenance suggests 'contact'.
        Warnings only: nothing is changed automatically."""
        w = []
        for v, ps in self.pre.items():
            for u in ps:
                key = tuple(sorted((u, v)))
                prov = self.provenance.get(key) or self.provenance.get(f"{u}-{v}")
                if prov and PROVENANCE_SUGGESTION.get(prov) == "contact":
                    w.append(f"relation {u}-{v} has provenance '{prov}' "
                             f"(suggests contact) but is encoded as a temporal "
                             f"prerequisite")
        return w

    def warnings(self): return list(self._warnings)

    # ----------------------------------------------------- legacy interface
    def n_orders_total(self):
        ns = [u for u in self.units if u != self.seed]
        return math.factorial(len(ns))

    def n_orders_permitted(self):
        """Seed-anchored singleton addition under prerequisite semantics.
        This is the ORIGINAL algorithm and must reproduce published counts."""
        ns = [u for u in self.units if u != self.seed]
        n = len(ns)
        if n > 25: raise ValueError(f"n = {n}: exact counting impractical beyond ~25")
        idx = {u: i for i, u in enumerate(ns)}
        req = [0]*n; exc = [0]*n
        for u, i in idx.items():
            req[i] = sum(1 << idx[q] for q in self.pre[u] if q in idx)
            m = 0
            for a, b in self.excludes:
                if u in a: m |= sum(1 << idx[q] for q in b if q in idx)
                if u in b: m |= sum(1 << idx[q] for q in a if q in idx)
            exc[i] = m
        dep = [0]*n
        for i in range(n):
            for j in range(n):
                if req[j] >> i & 1: dep[i] |= 1 << j
        memo = {0: 1}
        def f(S):
            v = memo.get(S)
            if v is not None: return v
            tot = 0
            for i in range(n):
                if not (S >> i) & 1: continue
                rest = S ^ (1 << i)
                if req[i] & ~rest: continue
                if exc[i] & rest: continue
                if dep[i] & rest: continue
                tot += f(rest)
            memo[S] = tot; return tot
        return f((1 << n) - 1)

    # ------------------------------------------------- assembly-tree layer
    def _masks(self):
        idx = {u: i for i, u in enumerate(self.units)}
        n = len(self.units)
        con = [0]*n; pre = [0]*n; exc = [0]*n
        for fs in self.contacts:
            a, b = tuple(fs) if len(fs) == 2 else (next(iter(fs)),)*2
            if a in idx and b in idx and a != b:
                con[idx[a]] |= 1 << idx[b]; con[idx[b]] |= 1 << idx[a]
        for v, ps in self.pre.items():
            if v in idx:
                pre[idx[v]] = sum(1 << idx[p] for p in ps if p in idx)
        for x, y in self.excludes:
            for u in x:
                if u in idx: exc[idx[u]] |= sum(1 << idx[q] for q in y if q in idx)
            for q in y:
                if q in idx: exc[idx[q]] |= sum(1 << idx[u] for u in x if u in idx)
        return idx, n, con, pre, exc

    def _connected(self, S, con):
        if S == 0: return False
        f = (S & -S).bit_length()-1; seen = 1 << f; st = [f]
        while st:
            i = st.pop(); nb = con[i] & S & ~seen
            while nb:
                j = (nb & -nb).bit_length()-1; seen |= 1 << j; st.append(j); nb &= nb-1
        return seen == S

    def _merge_ok(self, A, B, con, pre, exc):
        C = A | B; linked = False; m = A
        while m:
            i = (m & -m).bit_length()-1
            if con[i] & B: linked = True
            if exc[i] & B: return False
            m &= m-1
        if not linked: return False
        m = C
        while m:
            i = (m & -m).bit_length()-1
            if pre[i] & ~C: return False
            m &= m-1
        return True

    def n_trees(self, subset=None):
        """Exact count of admissible binary assembly trees (brute-force
        validated). F(S) = sum over splits of F(A)F(B)."""
        idx, n, con, pre, exc = self._masks()
        S0 = ((1 << n)-1) if subset is None else sum(1 << idx[u] for u in subset)
        memo = {}
        def F(S):
            if S & (S-1) == 0: return 1
            v = memo.get(S)
            if v is not None: return v
            tot = 0; low = S & -S; sub = (S-1) & S
            while sub:
                if sub & low:
                    A, B = sub, S ^ sub
                    if B and self._connected(A, con) and self._connected(B, con) \
                       and self._merge_ok(A, B, con, pre, exc):
                        tot += F(A)*F(B)
                sub = (sub-1) & S
            memo[S] = tot; return tot
        r = F(S0); self._states = len(memo); return r

    def formable(self, subset):
        """Can this subset form as an admissible subcomplex?"""
        idx, n, con, pre, exc = self._masks()
        S = sum(1 << idx[u] for u in subset if u in idx)
        return self._connected(S, con) and self.n_trees(subset) > 0

    def coverage(self):
        """Which declared observed subcomplexes are formable?"""
        return {frozenset(e): self.formable(e) for e in self.observed}
