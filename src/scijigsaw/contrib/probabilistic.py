"""Probabilistic layer over declared post-translational modifications.

A modification is a *declared* edit to the interface constraints (add/remove a
prerequisite edge, or add/remove an exclusion), carrying an occupancy: the
fraction of molecules in which it is present. Occupancies are evidence-declared
inputs, never inferred.

The deterministic kernel is reused unchanged: for each on/off assignment of the
modifications we build one Assembly and call its exact ``n_orders_permitted``.
Weighting the resulting counts by the occupancies turns a single order count
into a distribution. This is the implementation of Supplementary Section S5.4.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from itertools import product
from typing import Callable, Dict, Iterable, List, Tuple
import random

from ..assembly import Assembly


@dataclass
class Modification:
    """A declared modification and its effect on the constraint graph.

    name       : identifier (e.g. "SNAP25:pThr138").
    occupancy  : declared probability the site is modified, in [0, 1].
    apply      : function taking (requires, excludes) copies and mutating them
                 to encode the modification's declared effect when PRESENT.
    """
    name: str
    occupancy: float
    apply: Callable[[Dict[str, set], List[Tuple[set, set]]], None]

    def __post_init__(self):
        if not 0.0 <= self.occupancy <= 1.0:
            raise ValueError(f"occupancy for {self.name} must be in [0,1]")


def remove_prerequisite(a: str, b: str):
    """Modification effect: remove the requirement that `a` needs `b`."""
    def _apply(requires, excludes):
        if a in requires:
            requires[a] = [u for u in requires[a] if u != b]
    return _apply


def add_exclusion(x: str, y: str):
    """Modification effect: make `x` and `y` mutually exclusive."""
    def _apply(requires, excludes):
        excludes.append(({x}, {y}))
    return _apply


def _rebuild(base: Assembly, active: Iterable[Modification]) -> Assembly:
    """Construct one Assembly with the active modifications applied."""
    requires = {k: list(v) for k, v in base.requires.items()}
    excludes = [ (set(a), set(b)) for (a, b) in base.excludes ]
    for m in active:
        m.apply(requires, excludes)
    return Assembly(requires=requires, excludes=excludes, seed=base.seed)


@dataclass
class MixtureResult:
    instances: List[dict] = field(default_factory=list)  # per-assignment records
    expected_orders: float = 0.0
    exact: bool = True
    n_samples: int = 0

    def probability_nonzero(self) -> float:
        """Probability the assembly admits at least one permitted order."""
        return sum(r["weight"] for r in self.instances if r["orders"] > 0)


def order_distribution(base: Assembly,
                       mods: List[Modification],
                       target: Callable[[Assembly], bool] | None = None,
                       max_exact: int = 15,
                       n_samples: int = 4000,
                       rng_seed: int = 0) -> MixtureResult:
    """Distribution of permitted-order counts under partial occupancy.

    For k <= max_exact, enumerate all 2**k assignments exactly. For larger k,
    sample assignments from the declared occupancies (legitimate here: the draws
    come from a known distribution over declared states, not from an attempt to
    estimate an unknown fraction). Returns a MixtureResult with the expected
    order count and the per-instance records.

    If ``target`` is given, an instance contributes its order count only when the
    predicate holds for the rebuilt Assembly (e.g. a chosen target state remains
    feasible); otherwise it contributes zero. This yields the expected number of
    orders competent for that target, matching Supplementary Section S5.4.
    """
    k = len(mods)
    res = MixtureResult()

    if k <= max_exact:
        for bits in product((0, 1), repeat=k):
            active = [m for m, b in zip(mods, bits) if b]
            w = 1.0
            for m, b in zip(mods, bits):
                w *= m.occupancy if b else (1.0 - m.occupancy)
            if w == 0.0:
                continue
            inst = _rebuild(base, active)
            competent = (target is None) or bool(target(inst))
            orders = inst.n_orders_permitted() if competent else 0
            res.instances.append(
                {"present": [m.name for m in active], "weight": w,
                 "orders": orders, "target_competent": competent})
            res.expected_orders += w * orders
        res.exact = True
        return res

    # sampling for large k
    rng = random.Random(rng_seed)
    acc = 0.0
    for _ in range(n_samples):
        active = [m for m in mods if rng.random() < m.occupancy]
        inst = _rebuild(base, active)
        competent = (target is None) or bool(target(inst))
        orders = inst.n_orders_permitted() if competent else 0
        acc += orders
        res.instances.append(
            {"present": [m.name for m in active], "weight": 1.0 / n_samples,
             "orders": orders, "target_competent": competent})
    res.expected_orders = acc / n_samples
    res.exact = False
    res.n_samples = n_samples
    return res
