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

    def __init__(self, requires: Dict[str, Iterable[str]],
                 excludes: Iterable[Tuple[Set[str], Set[str]]] = (),
                 seed: str | None = None):
        self.requires = {k: set(v) for k, v in requires.items()}
        self.excludes = [(set(a), set(b)) for a, b in excludes]
        self.seed = seed
        self.units = list(self.requires)
        self._check_acyclic()

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

    # ---------------------------------------------------------------- counts
    def n_orders_total(self) -> int:
        return math.factorial(len(self.units))

    def n_orders_permitted(self) -> int:
        """Exact count of linear extensions, restricted by exclusion.

        f(S) = number of permitted orders seating exactly S. A unit u may be
        seated LAST in S iff R(u) is contained in S\\{u} (plus the seed) and no
        unit excluded with u lies in S\\{u}. Every permitted order of S ends in
        exactly one such u, so

            f(S) = sum over admissible-last u of f(S \\ {u}),   f({}) = 1

        partitions the permitted orders of S without overlap. O(2^n n).
        """
        units = self.units
        n = len(units)
        if n > 25:
            raise ValueError(
                f"n = {n}: exact counting is O(2^n n) and impractical beyond ~25 "
                "units. Use decomposition or bounded approximation.")
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

        f = [0] * (1 << n)
        f[0] = 1
        for S in range(1, 1 << n):
            tot = 0
            for i in range(n):
                if not (S >> i) & 1:
                    continue
                rest = S ^ (1 << i)
                if req[i] & ~rest:          # a requirement is absent
                    continue
                if exc[i] & rest:           # an excluded partner is present
                    continue
                tot += f[rest]
            f[S] = tot
        return f[(1 << n) - 1]

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
