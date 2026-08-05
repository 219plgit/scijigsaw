"""Exact enumeration of the assembly orders permitted by an encoded seating model.

THE MODEL, STATED PRECISELY
---------------------------
Let U be the units of an assembly and s a seed already present.

  PRECEDENCE  For each u, R(u) is the set of units u binds. u may be seated only
              in a state where R(u) is present. A *bridge* is any u with
              |R(u)| >= 2: it cannot be seated until the seam beneath it closes.
  EXCLUSION   Pairs {u, v} whose interfaces on a common partner overlap by more
              than a threshold cannot both be on the board.

A seating order is a permutation of U; it is PERMITTED if, read left to right,
every unit finds its requirements met and no excluded partner already present.

WHAT THIS IS NOT
----------------
The precedence rule requires ALL partners of a unit to be present before it is
seated. That is exact for seating a piece into a completed board. It is not
universally equivalent to biological assembly: a multivalent protein may bind one
partner, remain partially engaged, and only then encounter a second. The count
returned is therefore of orders permitted under this seating model -- a LOWER
BOUND on the biologically accessible set, not a characterisation of it.

WHY EXACT COUNTING, NOT SAMPLING
--------------------------------
Once the permitted fraction falls below ~1e-6 a uniformly random permutation
essentially never lands inside it; the estimate collapses to zero and the
reported reduction becomes spuriously enormous. We hit exactly this failure in
development, obtaining a fabricated six-billion-fold reduction from zero sampled
hits. Count; do not sample.
"""
from __future__ import annotations

import math
from functools import lru_cache
from typing import Dict, Iterable, Set, Tuple


class Assembly:
    """A constraint poset: precedence plus exclusion."""

    def __init__(self, requires: Dict[str, Iterable[str]] | None = None,
                 excludes: Iterable[Tuple[Set[str], Set[str]]] = (),
                 seed: str | None = None,
                 contacts: Iterable[Tuple[str, str]] = (),
                 prerequisites: Iterable[Tuple[str, str]] = (),
                 units: Iterable[str] | None = None,
                 provenance: Dict[Tuple[str, str], str] | None = None,
                 observed_subcomplexes: Iterable[Iterable[str]] = ()):
        """A constraint poset: precedence plus exclusion.

        `requires` retains its original meaning: every listed partner is a
        PREREQUISITE, and implicitly also a contact. Encodings and counts
        published under that model are unchanged.

        Typed relations may be declared instead of, or in addition to,
        `requires`:
          contacts       undirected (u, v): u and v can form an interface;
                         no claim about order.
          prerequisites  directed (u, v): u must precede v.
          units          components that must exist even if they appear in no
                         relation. Without this a declared component that has
                         no relation is silently absent from the assembly.
          provenance     {(u, v): evidence_kind}; raises warnings only, and
                         never reassigns a relation type.
        """
        self.requires = {k: set(v) for k, v in (requires or {}).items()}
        self.excludes = [(set(a), set(b)) for a, b in excludes]
        self.seed = seed
        self.provenance = dict(provenance or {})
        self.observed = [set(e) for e in observed_subcomplexes]

        # contacts implied by `requires`, plus any declared directly
        self._contacts = {frozenset((u, p))
                          for u, ps in self.requires.items() for p in ps}
        for u, v in contacts:
            self._contacts.add(frozenset((u, v)))
        for u, v in prerequisites:            # u must precede v
            # A temporal prerequisite is NOT evidence of a physical interface:
            # u may enable v indirectly. Declared prerequisites therefore do not
            # create contacts. (A legacy `requires` entry does imply both, for
            # backward compatibility.) If a prerequisite-only encoding is passed
            # to the tree layer, check_contact_graph() reports the missing
            # contacts rather than silently inventing them.
            self.requires.setdefault(v, set()).add(u)
            self.requires.setdefault(u, set())

        # components declared explicitly, or introduced by a typed relation,
        # become units. A legacy call (requires only) is unaffected.
        declared = set(units or ())
        if contacts or prerequisites:
            for fs in self._contacts:
                declared |= set(fs)
        for u in declared:
            self.requires.setdefault(u, set())

        self.units = list(self.requires)
        self._check_acyclic()
        self._warnings = self._provenance_warnings()

    def _check_acyclic(self):
        colour = {u: 0 for u in self.units}

        def visit(u):
            if colour[u] == 1:
                raise ValueError(f"cyclic precedence at {u!r}")
            if colour[u] == 2:
                return
            colour[u] = 1
            for p in self.requires[u]:
                if p in colour:
                    visit(p)
            colour[u] = 2

        for u in self.units:
            visit(u)

    _PROVENANCE_SUGGESTION = {
        "structure": "contact", "pisa": "contact",
        "predicted_structure": "contact", "crosslinking": "contact",
        "pulldown": "contact", "native_ms": "observed_subcomplex",
        "time_resolved": "prerequisite", "nucleated_cascade": "prerequisite",
        "curated": None,
    }

    def _provenance_warnings(self) -> list:
        """Flag prerequisites whose declared provenance suggests a contact.

        Warnings only: no relation is ever reassigned automatically."""
        out = []
        for v, ps in self.requires.items():
            for u in ps:
                kind = (self.provenance.get(tuple(sorted((u, v))))
                        or self.provenance.get((u, v)))
                if kind and self._PROVENANCE_SUGGESTION.get(kind) == "contact":
                    out.append(f"relation {u}-{v} has provenance {kind!r} "
                               f"(suggests contact) but is encoded as a "
                               f"temporal prerequisite")
        return out

    def warnings(self) -> list:
        return list(self._warnings)

    # ---------------------------------------------------------------- counts
    def n_orders_total(self) -> int:
        return math.factorial(len(self.units))

    def dependency_graph(self, include_seed=False):
        """Undirected graph of the precedence relation.

        For assemblies that are enumerated but not rendered as a tile board this
        is the graph a physical layout would have to realise. The seed is
        excluded by default because it is present before assembly begins;
        `include_seed=True` adds it as a physical component."""
        import networkx as nx
        G = nx.Graph()
        for unit, reqs in self.requires.items():
            G.add_node(unit)
            for r in reqs:
                if r == self.seed and not include_seed:
                    continue
                G.add_edge(unit, r)
        if include_seed and self.seed is not None:
            for unit, reqs in self.requires.items():
                if self.seed in reqs or not (set(reqs) - {self.seed}):
                    G.add_edge(unit, self.seed)
        return G

    def n_orders_permitted(self) -> int:
        """Exact count of linear extensions, restricted by exclusion.

        f(S) = number of permitted orders seating exactly S. A unit u may be
        seated LAST in S iff R(u) is contained in S\\{u} (plus the seed) and no
        unit excluded with u lies in S\\{u}. Every permitted order of S ends in
        exactly one such u, so

            f(S) = sum over admissible-last u of f(S \\ {u}),   f({}) = 1

        partitions the permitted orders of S without overlap.

        The recursion is evaluated lazily (memoised, top-down), so only states
        that are actually reachable under the precedence and exclusion rules are
        ever allocated. Reachable states are the order ideals of the poset: their
        number falls sharply as the poset deepens. The worst case (a poset with
        few relations, where every subset is reachable) is unchanged at O(2^n n)
        time and O(2^n) memory, which is what the guard below reflects.
        """
        units = self.units
        n = len(units)
        if n > 25:
            raise ValueError(
                f"n = {n}: exact counting is worst-case O(2^n n) and impractical "
                "beyond ~25 units. Lazy evaluation reduces memory on constrained "
                "posets but does not change the worst case, which is a poset with "
                "few relations. Use decomposition or bounded approximation.")
        idx = {u: i for i, u in enumerate(units)}
        req = [0] * n
        exc = [0] * n
        for u, i in idx.items():
            req[i] = sum(1 << idx[q] for q in self.requires[u] if q in idx)
            m = 0
            for a, b in self.excludes:
                if u in a:
                    m |= sum(1 << idx[q] for q in b if q in idx)
                if u in b:
                    m |= sum(1 << idx[q] for q in a if q in idx)
            exc[i] = m

        # dep[i] = units that REQUIRE i. If any of them is still present, i is not
        # maximal in S, and every order of S\\{i} would be unsatisfiable (f = 0).
        # Skipping those branches confines the recursion to the order ideals.
        dep = [0] * n
        for i in range(n):
            for j in range(n):
                if req[j] >> i & 1:
                    dep[i] |= 1 << j

        # Lazy (memoised) evaluation: allocate only reachable states. Recursion
        # depth is bounded by n (one unit is removed per level).
        memo = {0: 1}

        def f(S):
            v = memo.get(S)
            if v is not None:
                return v
            tot = 0
            for i in range(n):
                if not (S >> i) & 1:
                    continue
                rest = S ^ (1 << i)
                if req[i] & ~rest:          # a requirement is absent
                    continue
                if exc[i] & rest:           # an excluded partner is present
                    continue
                if dep[i] & rest:           # something left still requires i
                    continue
                tot += f(rest)
            memo[S] = tot
            return tot

        total = f((1 << n) - 1)
        self._states_visited = len(memo)
        return total

    def reduction(self) -> float:
        p = self.n_orders_permitted()
        if p == 0:
            raise ValueError("no permitted order: the constraints are unsatisfiable")
        return self.n_orders_total() / p

    # ------------------------------------------------------------- topology
    def depth(self) -> int:
        """Length of the longest chain. This, not n, governs the pruning."""
        @lru_cache(maxsize=None)
        def d(u):
            return 1 + max((d(p) for p in self.requires[u] if p in self.requires),
                           default=0)
        return max((d(u) for u in self.units), default=0)

    def width(self) -> int:
        @lru_cache(maxsize=None)
        def lvl(u):
            return 1 + max((lvl(p) for p in self.requires[u] if p in self.requires),
                           default=0)
        counts: Dict[int, int] = {}
        for u in self.units:
            counts[lvl(u)] = counts.get(lvl(u), 0) + 1
        return max(counts.values(), default=0)

    def bridges(self) -> list:
        """Units binding two or more partners: they need a closed seam."""
        return [u for u, r in self.requires.items() if len(r) >= 2]

    def summary(self) -> dict:
        return dict(n=len(self.units), depth=self.depth(), width=self.width(),
                    bridges=len(self.bridges()),
                    total=self.n_orders_total(),
                    permitted=self.n_orders_permitted(),
                    reduction=self.reduction())

    # ------------------------------------------------- assembly-tree layer
    def _tree_units(self):
        """Units for the tree layer. The seed IS a component here: a merger
        model needs the scaffold's interfaces, whereas the seed-anchored count
        treats it as present from the start and excludes it from the units."""
        u = list(self.units)
        if self.seed is not None and self.seed not in u:
            u.append(self.seed)
        return u

    def _tree_masks(self):
        tu = self._tree_units()
        idx = {u: i for i, u in enumerate(tu)}
        n = len(tu)
        con = [0] * n
        pre = [0] * n
        exc = [0] * n
        for fs in self._contacts:
            t = tuple(fs)
            if len(t) == 2 and t[0] in idx and t[1] in idx:
                con[idx[t[0]]] |= 1 << idx[t[1]]
                con[idx[t[1]]] |= 1 << idx[t[0]]
        for v, ps in self.requires.items():
            if v in idx:
                pre[idx[v]] = sum(1 << idx[p] for p in ps if p in idx)
        for x, y in self.excludes:
            for u in x:
                if u in idx:
                    exc[idx[u]] |= sum(1 << idx[q] for q in y if q in idx)
            for q in y:
                if q in idx:
                    exc[idx[q]] |= sum(1 << idx[u] for u in x if u in idx)
        return idx, n, con, pre, exc

    @staticmethod
    def _connected(S, con):
        if S == 0:
            return False
        f = (S & -S).bit_length() - 1
        seen = 1 << f
        stack = [f]
        while stack:
            i = stack.pop()
            nb = con[i] & S & ~seen
            while nb:
                j = (nb & -nb).bit_length() - 1
                seen |= 1 << j
                stack.append(j)
                nb &= nb - 1
        return seen == S

    @staticmethod
    def _merge_ok(A, B, con, pre, exc):
        C = A | B
        linked = False
        m = A
        while m:
            i = (m & -m).bit_length() - 1
            if con[i] & B:
                linked = True
            if exc[i] & B:
                return False
            m &= m - 1
        if not linked:
            return False
        m = C
        while m:
            i = (m & -m).bit_length() - 1
            if pre[i] & ~C:
                return False
            m &= m - 1
        return True

    def n_trees(self, subset: Iterable[str] | None = None) -> int:
        """Exact count of admissible binary assembly trees.

        A subcomplex is a connected subset of the contact graph; an event is a
        binary merger, permitted when a contact links the two parts, no
        exclusion is violated between them, and every prerequisite of every
        component of the merged subcomplex lies within it:

            F({i}) = 1,   F(S) = sum over unordered splits of F(A) F(B)

        The seed-anchored count of `n_orders_permitted` is recovered when every
        relation is a prerequisite and every merger adds one component to the
        subcomplex containing the seed. Verified against exhaustive enumeration
        on all connected contact graphs up to five components.
        """
        idx, n, con, pre, exc = self._tree_masks()
        S0 = ((1 << n) - 1) if subset is None \
            else sum(1 << idx[u] for u in subset if u in idx)
        memo = {}

        def F(S):
            if S & (S - 1) == 0:
                return 1
            v = memo.get(S)
            if v is not None:
                return v
            tot = 0
            low = S & -S
            sub = (S - 1) & S
            while sub:
                if sub & low:
                    A, B = sub, S ^ sub
                    if B and self._connected(A, con) \
                       and self._connected(B, con) \
                       and self._merge_ok(A, B, con, pre, exc):
                        tot += F(A) * F(B)
                sub = (sub - 1) & S
            memo[S] = tot
            return tot

        r = F(S0)
        self._tree_states = len(memo)
        return r

    def formable(self, subset: Iterable[str]) -> bool:
        """Can this subset form as an admissible subcomplex?"""
        idx, n, con, pre, exc = self._tree_masks()
        S = sum(1 << idx[u] for u in subset if u in idx)
        return self._connected(S, con) and self.n_trees(subset) > 0

    # ------------------------------------------- experimental support queries
    def _trees_containing(self, E_mask, S, con, pre, exc, memoF, memoG):
        """G(S) = number of admissible trees producing S in which the subset E
        appears as a subcomplex (an internal node or leaf of the tree).

        A proper subset E can lie on at most one side of any root split, because
        the two sides are disjoint. Hence

            G(S) = sum over splits {A,B} of
                     G(A)*F(B)  if E subset A
                   + F(A)*G(B)  if E subset B
                   + F(A)*F(B)  if S == E   (the root node itself is E)

        F is the unrestricted tree count computed by the same recursion."""
        if S == E_mask:
            return self._F(S, con, pre, exc, memoF)
        v = memoG.get(S)
        if v is not None:
            return v
        if E_mask & ~S:                       # E not contained in S at all
            memoG[S] = 0
            return 0
        if S & (S - 1) == 0:                  # singleton, and S != E
            memoG[S] = 0
            return 0
        tot = 0
        low = S & -S
        sub = (S - 1) & S
        while sub:
            if sub & low:
                A, B = sub, S ^ sub
                if B and self._connected(A, con) and self._connected(B, con) \
                   and self._merge_ok(A, B, con, pre, exc):
                    if not (E_mask & ~A):     # E subset of A
                        tot += self._trees_containing(E_mask, A, con, pre, exc,
                                                      memoF, memoG) \
                               * self._F(B, con, pre, exc, memoF)
                    elif not (E_mask & ~B):   # E subset of B
                        tot += self._F(A, con, pre, exc, memoF) \
                               * self._trees_containing(E_mask, B, con, pre,
                                                        exc, memoF, memoG)
            sub = (sub - 1) & S
        memoG[S] = tot
        return tot

    def _F(self, S, con, pre, exc, memo):
        if S & (S - 1) == 0:
            return 1
        v = memo.get(S)
        if v is not None:
            return v
        tot = 0
        low = S & -S
        sub = (S - 1) & S
        while sub:
            if sub & low:
                A, B = sub, S ^ sub
                if B and self._connected(A, con) and self._connected(B, con) \
                   and self._merge_ok(A, B, con, pre, exc):
                    tot += self._F(A, con, pre, exc, memo) \
                           * self._F(B, con, pre, exc, memo)
            sub = (sub - 1) & S
        memo[S] = tot
        return tot

    def support(self, subset) -> dict:
        """Exact fraction of admissible assembly trees in which `subset` appears
        as a subcomplex.

        Returns n_containing, n_total and support = n_containing / n_total, with
        a categorical reading: support 0 means the subcomplex is impossible under
        the declared model; 0 < support < 1 means it is optional, occurring on
        some admissible routes and not others; support 1 means every admissible
        assembly passes through it, so it is necessary.
        """
        idx, n, con, pre, exc = self._tree_masks()
        full = (1 << n) - 1
        E = sum(1 << idx[u] for u in subset if u in idx)
        memoF, memoG = {}, {}
        total = self._F(full, con, pre, exc, memoF)
        if total == 0:
            return dict(subset=sorted(subset), n_containing=0, n_total=0,
                        support=float("nan"), reading="no admissible tree")
        if not self._connected(E, con):
            k = 0
        else:
            k = self._trees_containing(E, full, con, pre, exc, memoF, memoG)
        f = k / total
        reading = ("impossible" if k == 0 else
                   "necessary" if k == total else "optional")
        return dict(subset=sorted(subset), n_containing=k, n_total=total,
                    support=f, reading=reading)

    def support_table(self) -> list:
        """Exact support for every declared observed subcomplex."""
        return [self.support(e) for e in self.observed]

    def coverage(self) -> dict:
        """Which declared observed subcomplexes are formable?"""
        return {frozenset(e): self.formable(e) for e in self.observed}

    def check_contact_graph(self, strict: bool = False):
        """Tree enumeration needs a connected contact graph; the seed-anchored
        count does not. Returns a message, or None if the graph is usable."""
        idx, n, con, pre, exc = self._tree_masks()
        iso = [u for u in self._tree_units() if con[idx[u]] == 0]
        msg = None
        if iso:
            msg = (f"contact graph has isolated units {sorted(iso)}: no merger "
                   f"can incorporate them, so tree enumeration yields 0. "
                   f"Supply contact evidence for these units.")
        elif n and not self._connected((1 << n) - 1, con):
            msg = ("contact graph is disconnected: tree enumeration yields 0. "
                   "Supply contacts linking the components.")
        if msg and strict:
            raise ValueError(msg)
        return msg

    def diagnose(self):
        """Explain a zero tree count, or return None if trees exist."""
        m = self.check_contact_graph()
        if m:
            return m
        if self.n_trees() > 0:
            return None
        if self.excludes:
            return ("no admissible assembly tree: the declared exclusions "
                    "leave no compatible merger sequence")
        return "no admissible assembly tree under the declared relations"
