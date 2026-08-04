"""Typed-edge Assembly: contacts vs prerequisites, with legacy compatibility."""
from __future__ import annotations
import math
from typing import Dict, Set

PROVENANCE_SUGGESTION = {
    "structure": "contact", "pisa": "contact", "predicted_structure": "contact",
    "crosslinking": "contact", "pulldown": "contact",
    "native_ms": "observed_subcomplex",
    "time_resolved": "prerequisite", "nucleated_cascade": "prerequisite",
    "curated": None,
}


class TypedAssembly:
    def __init__(self, requires=None, contacts=(), prerequisites=(),
                 excludes=(), seed=None, provenance=None,
                 observed_subcomplexes=(), units=None):
        self.seed = seed
        self.provenance = dict(provenance or {})
        self.observed = [set(e) for e in observed_subcomplexes]
        declared_units = list(units) if units else None
        units = set()
        self.contacts: Set[frozenset] = set()
        self.pre: Dict[str, Set[str]] = {}

        if requires:
            for u, ps in requires.items():
                units.add(u); self.pre.setdefault(u, set())
                for p in ps:
                    units.add(p); self.pre[u].add(p)
                    self.contacts.add(frozenset((u, p)))
        for u, v in contacts:
            units.update((u, v)); self.contacts.add(frozenset((u, v)))
        for u, v in prerequisites:            # u must precede v
            units.update((u, v))
            self.pre.setdefault(v, set()).add(u)
            self.contacts.add(frozenset((u, v)))

        if seed: units.add(seed)
        if declared_units: units |= set(declared_units)
        self.units = sorted(units)
        for u in self.units: self.pre.setdefault(u, set())
        self.excludes = [(set(a), set(b)) for a, b in excludes]
        self._check_acyclic()
        self._warnings = self._consistency_check()

    # ---------------------------------------------------------------- checks
    def _check_acyclic(self):
        colour = {u: 0 for u in self.units}
        def visit(u, path):
            if colour[u] == 1:
                raise ValueError("cyclic prerequisite: " +
                                 " -> ".join(path[path.index(u):] + [u]))
            if colour[u] == 2: return
            colour[u] = 1
            for p in self.pre.get(u, ()):
                if p in colour: visit(p, path + [u])
            colour[u] = 2
        for u in self.units: visit(u, [])

    def _consistency_check(self):
        w = []
        for v, ps in self.pre.items():
            for u in ps:
                key = tuple(sorted((u, v)))
                prov = self.provenance.get(key) or self.provenance.get(f"{u}-{v}")
                if prov and PROVENANCE_SUGGESTION.get(prov) == "contact":
                    w.append(f"relation {u}-{v} has provenance '{prov}' "
                             f"(suggests contact) but is encoded as a prerequisite")
        return w

    def warnings(self): return list(self._warnings)

    def _isolated_units(self):
        idx, n, con, pre, exc = self._masks()
        return [u for u in self.units if con[idx[u]] == 0]

    def check_contact_graph(self, strict=False):
        iso = self._isolated_units()
        idx, n, con, pre, exc = self._masks()
        connected = self._connected((1 << n) - 1, con) if n else False
        msg = None
        if iso:
            msg = (f"contact graph has isolated units {iso}: no merger can "
                   f"incorporate them, so tree enumeration yields 0.")
        elif not connected:
            msg = "contact graph is disconnected: tree enumeration yields 0."
        if msg and strict: raise ValueError(msg)
        return msg

    def diagnose(self):
        m = self.check_contact_graph()
        if m: return m
        if self.n_trees() > 0: return None
        if self.excludes:
            return "no admissible tree: exclusions leave no compatible merger"
        return "no admissible assembly tree under the declared relations"

    # ------------------------------------------------------ legacy interface
    def n_orders_total(self):
        return math.factorial(len([u for u in self.units if u != self.seed]))

    def n_orders_permitted(self):
        ns = [u for u in self.units if u != self.seed]
        n = len(ns)
        if n > 25: raise ValueError(f"n = {n}: impractical beyond ~25")
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

    # -------------------------------------------------- assembly-tree layer
    def _masks(self):
        idx = {u: i for i, u in enumerate(self.units)}
        n = len(self.units)
        con = [0]*n; pre = [0]*n; exc = [0]*n
        for fs in self.contacts:
            t = tuple(fs)
            if len(t) == 2 and t[0] in idx and t[1] in idx:
                con[idx[t[0]]] |= 1 << idx[t[1]]
                con[idx[t[1]]] |= 1 << idx[t[0]]
        for v, ps in self.pre.items():
            if v in idx: pre[idx[v]] = sum(1 << idx[p] for p in ps if p in idx)
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
                j = (nb & -nb).bit_length()-1
                seen |= 1 << j; st.append(j); nb &= nb-1
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
        idx, n, con, pre, exc = self._masks()
        S = sum(1 << idx[u] for u in subset if u in idx)
        return self._connected(S, con) and self.n_trees(subset) > 0

    def coverage(self):
        return {frozenset(e): self.formable(e) for e in self.observed}
