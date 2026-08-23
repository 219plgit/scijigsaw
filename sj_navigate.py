"""
sj_navigate.py -- hypothesis testing and navigation on top of sj_prob.py.

Wraps the validated machinery in the four operations the paper promises:
  see   : condition on an observation, renormalise, keep querying
  do    : intervene on the ENCODING (edit typed relations), get a new session
  ask   : probability of a hypothesis, with verdict (impossible / possible /
          necessary / weighted), and named cause when impossible
  walk  : stand at an assembly state and inspect branch probabilities,
          remaining entropy, and most probable continuations

Hypotheses are built from readable constructors instead of lambdas:
  forms({"a","b"}, seed) -- subcomplex ever assembled
  before("a","b")        -- a joins strictly before b
  at_step(2,"b")         -- b is the 2nd non-seed addition
  last("d")              -- d joins last
  AND(h1,h2), OR(h1,h2), NOT(h)

Seeded mode only (navigation is over histories). Spaces are materialised as
leaf lists, which is exact and fine for constrained boards (hundreds to a few
thousand histories); for weakly constrained large boards, prune first.

Usage: put this file next to sj_prob.py, then `python3 sj_navigate.py` runs a
demo session on a placeholder encoding.
"""

from __future__ import annotations

import math
from typing import Callable, Dict, FrozenSet, List, Optional, Sequence, Tuple

from sj_prob import (TypedModel, calibrated_tree, leaf_masses, drop_relations,
                     counterfactual_compare)

Component = str
Subset = FrozenSet[Component]
History = Tuple[Component, ...]
Hypothesis = Callable[[History], bool]


# ----------------------------------------------------------------- hypotheses
def contains_components(components, seed: Optional[Component] = None) -> Hypothesis:
    """All queried components have joined by some point (CONTAINMENT). NOTE: in a
    complete seeded history every component eventually joins, so this is trivially
    true for any subset of V; it is useful mainly inside compound hypotheses. It
    does NOT test exact intermediate formation -- use forms_exact for that."""
    E = frozenset(components)

    def h(path: History, _E=E, _s=seed) -> bool:
        cur = {_s} if _s else set()
        if _E <= cur:
            return True
        for p in path:
            cur.add(p)
            if _E <= cur:
                return True
        return False
    h.label = f"contains({sorted(E)})"
    return h


def forms_exact(components, seed: Optional[Component] = None) -> Hypothesis:
    """The EXACT component set existed as an intermediate: at some step the
    assembled state equals the queried set. In seed-anchored mode every
    intermediate contains the seed, so an exact intermediate omitting the seed is
    impossible by construction -- this predicate expresses the eIF3/SNARE
    representation-failure claim."""
    E = frozenset(components)

    def h(path: History, _E=E, _s=seed) -> bool:
        cur = {_s} if _s else set()
        if cur == _E:
            return True
        for p in path:
            cur.add(p)
            if cur == _E:
                return True
        return False
    h.label = f"forms_exact({sorted(E)})"
    return h


forms = contains_components  # backward-compatible alias (containment semantics)


def before(a: Component, b: Component) -> Hypothesis:
    def h(path: History) -> bool:
        return a in path and b in path and path.index(a) < path.index(b)
    h.label = f"before({a},{b})"
    return h


def at_step(k: int, p: Component) -> Hypothesis:
    def h(path: History) -> bool:
        return len(path) >= k and path[k - 1] == p
    h.label = f"at_step({k},{p})"
    return h


def last(p: Component) -> Hypothesis:
    def h(path: History) -> bool:
        return bool(path) and path[-1] == p
    h.label = f"last({p})"
    return h


def AND(*hs: Hypothesis) -> Hypothesis:
    def h(path: History) -> bool:
        return all(x(path) for x in hs)
    h.label = " & ".join(getattr(x, "label", "?") for x in hs)
    return h


def OR(*hs: Hypothesis) -> Hypothesis:
    def h(path: History) -> bool:
        return any(x(path) for x in hs)
    h.label = " | ".join(getattr(x, "label", "?") for x in hs)
    return h


def NOT(x: Hypothesis) -> Hypothesis:
    def h(path: History) -> bool:
        return not x(path)
    h.label = f"not({getattr(x, 'label', '?')})"
    return h


# ------------------------------------------------------------------ session
class Navigator:
    """A session over one declared model: a measure on admissible histories,
    plus a record of every conditioning step taken to reach it."""

    def __init__(self, model: TypedModel, seed: Component,
                 step_weight=None, _masses=None, _trail=None):
        self.model = model
        self.seed = seed
        self.step_weight = step_weight
        if _masses is None:
            tree = calibrated_tree(model, seed, step_weight=step_weight)
            _masses = [(p, m) for p, m in leaf_masses(tree) if m > 0]
        self.masses: List[Tuple[History, float]] = _masses
        self.trail: List[str] = _trail or []

    # ---- basic quantities ---------------------------------------------------
    def total(self) -> float:
        return sum(m for _, m in self.masses)

    def n_histories(self) -> int:
        return len(self.masses)

    def entropy_bits(self) -> float:
        z = self.total()
        ms = [m / z for _, m in self.masses if m > 0]
        return -sum(m * math.log2(m) for m in ms) if ms else float("nan")

    def n_eff(self) -> float:
        return 2 ** self.entropy_bits()

    def contains(self, components) -> Hypothesis:
        return contains_components(components, seed=self.seed)

    def forms(self, components) -> Hypothesis:  # alias of contains (see docstring)
        return contains_components(components, seed=self.seed)

    def forms_exact(self, components) -> Hypothesis:
        """Seed-aware exact-intermediate predicate."""
        return forms_exact(components, seed=self.seed)

    # ---- ask: hypothesis probability with verdict ---------------------------
    def ask(self, h: Hypothesis) -> Dict[str, object]:
        z = self.total()
        p = sum(m for path, m in self.masses if h(path)) / z if z else 0.0
        if p == 0.0:
            verdict = "impossible (under the declared model and conditioning)"
            cause = self._diagnose(h)
        elif p == 1.0:
            verdict, cause = "necessary", None
        else:
            verdict, cause = f"possible; probability mass = {p:.3f}", None
        return {"hypothesis": getattr(h, "label", "?"), "prob": p,
                "verdict": verdict, "cause": cause, "given": list(self.trail)}

    def _diagnose(self, h: Hypothesis) -> str:
        """Name why a hypothesis has zero mass: conditioning, or a typed relation."""
        base = Navigator(self.model, self.seed, self.step_weight)
        if any(h(path) for path, _ in base.masses):
            return ("excluded by conditioning: " + "; ".join(self.trail)
                    if self.trail else "excluded (no admissible history)")
        # not in the unconditioned space either: try single-relation ablations
        for (a, b) in self.model.P:
            relaxed = drop_relations(self.model, drop_P=[(a, b)])
            if any(h(path) for path, _ in
                   Navigator(relaxed, self.seed).masses):
                return f"blocked by prerequisite {a}->{b}"
        for e in self.model.X:
            relaxed = drop_relations(self.model, drop_X=[e])
            if any(h(path) for path, _ in
                   Navigator(relaxed, self.seed).masses):
                u, v = sorted(e)
                return f"blocked by exclusion {u}--{v}"
        return ("not expressible under the seed-anchored representation "
                "(consider merger-tree mode) or blocked by relation combination")

    # ---- see: condition and continue ---------------------------------------
    def see(self, h: Hypothesis) -> "Navigator":
        kept = [(p, m) for p, m in self.masses if h(p)]
        if not kept:
            raise ValueError(f"conditioning on zero-mass event: "
                             f"{getattr(h, 'label', '?')}")
        return Navigator(self.model, self.seed, self.step_weight,
                         _masses=kept,
                         _trail=self.trail + [getattr(h, "label", "?")])

    # ---- do: intervene on the encoding --------------------------------------
    def do(self, *, drop_P=(), drop_X=(), drop_C=(),
           add_P=(), add_X=(), add_C=()) -> "Navigator":
        alt = drop_relations(self.model, drop_P=drop_P, drop_X=drop_X,
                             drop_C=drop_C)
        alt = TypedModel(
            V=alt.V,
            C=alt.C | frozenset(frozenset(e) for e in add_C),
            P=tuple(alt.P) + tuple(add_P),
            X=alt.X | frozenset(frozenset(e) for e in add_X),
            eta=alt.eta)
        return Navigator(alt, self.seed, self.step_weight,
                         _trail=[f"do: -P{list(drop_P)} +P{list(add_P)} "
                                 f"-X{[sorted(e) for e in drop_X]} +X{[sorted(e) for e in add_X]} "
                                 f"-C{[sorted(e) for e in drop_C]} +C{[sorted(e) for e in add_C]}"])

    def compare_with(self, other: "Navigator",
                     query: Optional[Subset] = None) -> Dict[str, object]:
        return counterfactual_compare(self.model, other.model, self.seed,
                                      step_weight=self.step_weight,
                                      query=query)

    # ---- walk: stand at a state, look around --------------------------------
    def walk(self, added: Sequence[Component] = ()) -> Dict[str, object]:
        prefix = tuple(added)
        sub = [(p, m) for p, m in self.masses if p[:len(prefix)] == prefix]
        z = sum(m for _, m in sub)
        if z == 0:
            return {"state": [self.seed, *prefix], "reachable": False}
        nxt: Dict[Component, float] = {}
        for p, m in sub:
            if len(p) > len(prefix):
                nxt[p[len(prefix)]] = nxt.get(p[len(prefix)], 0.0) + m / z
        ms = [m / z for _, m in sub]
        h2 = -sum(m * math.log2(m) for m in ms if m > 0)
        hb = -sum(p * math.log2(p) for p in nxt.values() if p > 0)
        return {"state": [self.seed, *prefix], "reachable": True,
                "mass_here": z / self.total(),
                "branches": dict(sorted(nxt.items(), key=lambda kv: -kv[1])),
                "branch_entropy_bits": hb,
                "remaining_entropy_bits": h2,
                "remaining_n_eff": 2 ** h2}

    def compare(self, hs: Sequence[Hypothesis]) -> List[Dict[str, object]]:
        """Side-by-side hypothesis comparison, sorted by probability mass."""
        rows = [self.ask(h) for h in hs]
        return sorted(rows, key=lambda r: -r["prob"])

    def uncertain_branches(self, k: int = 5, min_mass: float = 0.02
                           ) -> List[Dict[str, object]]:
        """Rank partial assembly states by residual uncertainty (bits), weighted
        by the mass reaching them: where competing mechanism families diverge,
        and hence where a candidate discriminating experiment would look."""
        prefixes = {}
        z = self.total()
        for path, m in self.masses:
            for i in range(len(path)):
                prefixes.setdefault(path[:i], 0.0)
        rows = []
        for pre in prefixes:
            w = self.walk(pre)
            if not w["reachable"] or w["mass_here"] < min_mass:
                continue
            if len(w["branches"]) < 2:
                continue
            bp = list(w["branches"].values())
            h2b = -sum(p * math.log2(p) for p in bp if p > 0)  # branch-level entropy
            rows.append({"state": w["state"], "mass_here": w["mass_here"],
                         "branch_entropy_bits": h2b,
                         "weighted_uncertainty": w["mass_here"] * h2b,
                         "branches": w["branches"]})
        rows.sort(key=lambda r: -r["weighted_uncertainty"])
        return rows[:k]

    def top(self, k: int = 5) -> List[Tuple[History, float]]:
        z = self.total()
        return sorted(((p, m / z) for p, m in self.masses),
                      key=lambda t: -t[1])[:k]

    # ---- report -------------------------------------------------------------
    def report(self) -> str:
        lines = [f"model: |V|={len(self.model.V)}, seed={self.seed}, "
                 f"|P|={len(self.model.P)}, |X|={len(self.model.X)}",
                 f"conditioning trail: {self.trail or '(none)'}",
                 f"histories: {self.n_histories()}   "
                 f"H2={self.entropy_bits():.3f} bits   "
                 f"N_eff={self.n_eff():.2f}",
                 "most probable histories:"]
        for p, m in self.top(3):
            lines.append(f"  {self.seed}->" + "->".join(p) + f"   {m:.4f}")
        return "\n".join(lines)


# ------------------------------------------------------------------ demo
if __name__ == "__main__":
    V = ("Prt1", "Tif34", "Tif35", "Nip1", "Tif32")
    C = frozenset(frozenset(e) for e in [
        ("Prt1", "Tif34"), ("Prt1", "Tif35"), ("Tif35", "Tif34"),
        ("Prt1", "Nip1"), ("Nip1", "Tif32")])
    model = TypedModel(V=V, C=C, P=(("Nip1", "Tif32"),))
    nav = Navigator(model, "Prt1")

    print("== session ==")
    print(nav.report())
    print("\n== ask ==")
    for h in [nav.forms_exact({"Prt1", "Tif34"}),
              nav.forms_exact({"Tif35", "Tif34"}),   # off-seed: impossible by representation
              before("Tif34", "Tif35"),
              last("Tif32"), at_step(1, "Tif32")]:
        r = nav.ask(h)
        print(f"  {r['hypothesis']:<28} p={r['prob']:.3f}  {r['verdict']}"
              + (f"  [{r['cause']}]" if r["cause"] else ""))
    print("\n== see: condition on 'Tif34 before Tif35' ==")
    nav2 = nav.see(before("Tif34", "Tif35"))
    print(nav2.report())
    r = nav2.ask(last("Tif32"))
    print(f"  {r['hypothesis']:<28} p={r['prob']:.3f}  {r['verdict']}")
    print("\n== walk: after adding Nip1 ==")
    w = nav.walk(["Nip1"])
    print(f"  state {w['state']}  mass here {w['mass_here']:.3f}  "
          f"remaining N_eff {w['remaining_n_eff']:.2f}")
    print("  branches: " + ", ".join(f"{k}:{v:.3f}"
                                     for k, v in w["branches"].items()))
    print("\n== compare hypotheses ==")
    for r in nav.compare([before("Tif34", "Tif35"), before("Tif35", "Tif34"),
                          last("Tif32"), last("Tif35")]):
        print(f"  {r['hypothesis']:<24} p={r['prob']:.3f}")
    print("\n== uncertain branches (candidate discriminating experiments) ==")
    for r in nav.uncertain_branches(3):
        print(f"  state {'->'.join(r['state']):<28} mass {r['mass_here']:.3f}  "
              f"branch H2 {r['branch_entropy_bits']:.3f} bits  "
              f"-> {', '.join(f'{c}:{p:.2f}' for c, p in r['branches'].items())}")
    print("\n== do: relax the Nip1->Tif32 prerequisite, compare ==")
    nav3 = nav.do(drop_P=[("Nip1", "Tif32")])
    cmp_ = nav.compare_with(nav3, query=frozenset({"Prt1", "Tif34"}))
    print(f"  histories {cmp_['n_factual']} -> {cmp_['n_alternative']}, "
          f"dH2 {cmp_['delta_H2_bits']:+.3f} bits, "
          f"dPr(query) {cmp_['query_support_delta']:+.3f}")
    r = nav3.ask(at_step(1, "Tif32"))
    print(f"  under do(): {r['hypothesis']:<20} p={r['prob']:.3f}  {r['verdict']}")
