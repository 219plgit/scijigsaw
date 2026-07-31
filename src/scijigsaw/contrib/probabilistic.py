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
        """P_competent: probability the target admits at least one order."""
        return sum(r["weight"] for r in self.instances if r["orders"] > 0)

    # --- summary quantities (Supplementary S5.5) ---
    def p_competent(self) -> float:
        """Alias for probability_nonzero, in the reviewer's notation."""
        return self.probability_nonzero()

    def expected_fraction(self, total_permutations: int) -> float:
        """E[q_T] = sum_x Pr(x) L_T(G_x) / N, the expected admissible fraction.

        `total_permutations` is N = n! for the target component set (e.g. 5040
        for the 7-component VAMP2 fusion state). More comparable across graphs
        than the raw expected count.
        """
        return self.expected_orders / total_permutations

    def variance(self) -> float:
        """Var(L_T) = sum_x Pr(x) (L_T - E[L_T])^2.

        Nonzero variance makes explicit that no single instance has the
        expected (e.g. fractional) count: each instance is in one graph state.
        """
        mu = self.expected_orders
        return sum(r["weight"] * (r["orders"] - mu) ** 2 for r in self.instances)

    def std(self) -> float:
        return self.variance() ** 0.5

    def expected_given_competent(self) -> float:
        """E[L_T | T competent] = E[L_T] / P_competent (0 if never competent)."""
        pc = self.probability_nonzero()
        return (self.expected_orders / pc) if pc > 0 else 0.0

    def summary(self, total_permutations: int | None = None) -> dict:
        """All summary quantities in one dict (for reporting/reproduction)."""
        out = {
            "P_competent": self.p_competent(),
            "E_orders": self.expected_orders,
            "Var_orders": self.variance(),
            "SD_orders": self.std(),
            "E_orders_given_competent": self.expected_given_competent(),
            "exact": self.exact,
        }
        if total_permutations:
            out["E_fraction"] = self.expected_fraction(total_permutations)
        return out


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


def order_distribution_joint(base: Assembly,
                             mods: List[Modification],
                             joint: Dict[Tuple[int, ...], float],
                             target: Callable[[Assembly], bool] | None = None
                             ) -> MixtureResult:
    """Mixture over graphs with a *declared joint* distribution over states.

    `mods` gives the k modifications (their `occupancy` fields are ignored here);
    `joint` maps an assignment tuple x in {0,1}^k to a declared probability w_x,
    which need not factorise. This implements correlated modifications
    (Supplementary S5.5, item F): E[L_T] = sum_x w_x L_T(G_x). Weights must be
    non-negative and sum to 1 (within tolerance).
    """
    total = sum(joint.values())
    if abs(total - 1.0) > 1e-6:
        raise ValueError(f"joint weights must sum to 1 (got {total})")
    res = MixtureResult()
    for x, w in joint.items():
        if len(x) != len(mods):
            raise ValueError("assignment length must equal number of mods")
        if w < 0:
            raise ValueError("joint weights must be non-negative")
        if w == 0.0:
            continue
        active = [m for m, b in zip(mods, x) if b]
        inst = _rebuild(base, active)
        competent = (target is None) or bool(target(inst))
        orders = inst.n_orders_permitted() if competent else 0
        res.instances.append(
            {"present": [m.name for m in active], "weight": w,
             "orders": orders, "target_competent": competent})
        res.expected_orders += w * orders
    res.exact = True
    return res


def occupancy_uncertainty_analytic(total_orders: int,
                                   alpha: float,
                                   beta: float) -> dict:
    """Analytic uncertainty for the single-modification binary example.

    For the binary switch L = total_orders * (1 - p) with p ~ Beta(alpha, beta):
    because E is linear in p, E[L] = total_orders * (1 - alpha/(alpha+beta))
    exactly, with no simulation. The 95% equal-tailed interval maps the Beta
    quantiles of p through the (decreasing) transform. P_competent = 1 - p, so
    E[P_competent] = beta/(alpha+beta).

    Returns exact mean plus an equal-tailed 95% interval computed from Beta
    quantiles (bisection on the regularised incomplete Beta; no SciPy needed).

    Also returns the law-of-total-variance decomposition
    Var(L) = E_p[Var(L|p)] + Var_p[E(L|p)]. Here Var(L|p) is the population
    (Bernoulli-mixture) variance total_orders^2 * p(1-p) at fixed p, and
    E(L|p) = total_orders*(1-p); both terms are given in closed form.
    """
    import math

    def _betacdf(x, a, b, steps=20000):
        if x <= 0:
            return 0.0
        if x >= 1:
            return 1.0
        logB = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
        B = math.exp(logB)
        h = x / steps
        s = 0.0
        for i in range(steps):
            t = (i + 0.5) * h
            s += t ** (a - 1) * (1 - t) ** (b - 1)
        return s * h / B

    def _betaquantile(q, a, b):
        lo, hi = 0.0, 1.0
        for _ in range(80):
            m = 0.5 * (lo + hi)
            if _betacdf(m, a, b) < q:
                lo = m
            else:
                hi = m
        return 0.5 * (lo + hi)

    mean_p = alpha / (alpha + beta)
    E_L = total_orders * (1 - mean_p)                 # exact (E linear in p)
    E_Pcomp = beta / (alpha + beta)
    q025, q975 = _betaquantile(0.025, alpha, beta), _betaquantile(0.975, alpha, beta)
    L_lo, L_hi = total_orders * (1 - q975), total_orders * (1 - q025)

    # law of total variance for L = N*Bernoulli-mixture, treating each instance
    # as 252-or-0: Var(L|p) = N^2 p(1-p); E(L|p) = N(1-p).
    N = total_orders
    var_p = alpha * beta / ((alpha + beta) ** 2 * (alpha + beta + 1))
    E_pp = mean_p
    # E_p[p(1-p)] = E[p] - E[p^2] = mean_p - (var_p + mean_p^2)
    E_p_p1mp = E_pp - (var_p + E_pp ** 2)
    within = N ** 2 * E_p_p1mp            # E_p[Var(L|p)]
    between = N ** 2 * var_p              # Var_p[E(L|p)] = N^2 Var(p)
    return {
        "alpha": alpha, "beta": beta,
        "E_orders": E_L, "E_P_competent": E_Pcomp,
        "L_ci95_equal_tailed": (L_lo, L_hi),
        "Var_within_E_p_Var": within,
        "Var_between_Var_p_E": between,
        "Var_total": within + between,
        "method": "analytic (Beta quantiles by bisection); no simulation",
    }


def occupancy_uncertainty(base: Assembly,
                          make_mod: Callable[[float], Modification],
                          alpha: float,
                          beta: float,
                          target: Callable[[Assembly], bool] | None = None,
                          n_draws: int = 2000,
                          rng_seed: int = 0) -> dict:
    """Monte Carlo propagation of Beta(alpha, beta) occupancy uncertainty.

    General fallback for cases without a closed form (e.g. several modifications,
    or a nonlinear target). For the single binary switch prefer
    ``occupancy_uncertainty_analytic``, which is exact. Reports the seed, number
    of draws, and Monte Carlo standard error so the estimate is not mistaken for
    the analytic value.
    """
    rng = random.Random(rng_seed)
    E_draws, Pc_draws = [], []
    for _ in range(n_draws):
        p = rng.betavariate(alpha, beta)
        m = make_mod(p)
        r = order_distribution(base, [m], target=target)
        E_draws.append(r.expected_orders)
        Pc_draws.append(r.probability_nonzero())

    def _stats(v):
        v_sorted = sorted(v)
        mean = sum(v) / len(v)
        var = sum((x - mean) ** 2 for x in v) / (len(v) - 1)
        mcse = (var / len(v)) ** 0.5
        lo = v_sorted[int(0.025 * len(v))]
        hi = v_sorted[min(len(v) - 1, int(0.975 * len(v)))]
        return mean, (lo, hi), mcse

    Em, Eci, Emcse = _stats(E_draws)
    Pm, Pci, Pmcse = _stats(Pc_draws)
    return {
        "alpha": alpha, "beta": beta, "n_draws": n_draws, "rng_seed": rng_seed,
        "method": "Monte Carlo",
        "E_orders_mean": Em, "E_orders_ci95_percentile": Eci, "E_orders_mcse": Emcse,
        "P_competent_mean": Pm, "P_competent_ci95_percentile": Pci,
        "P_competent_mcse": Pmcse,
    }


def order_distribution_counterfactual(
        histories: Dict[str, Assembly],
        prior: Dict[str, float],
        target: Callable[[Assembly], bool] | None = None) -> MixtureResult:
    """Distribution over admissible order-space under alternative *histories*.

    A counterfactual hypothesis is a declared alternative constraint graph --- a
    duplication that did not happen, an interface repurposed differently, a
    socket shared rather than distinct. `histories` maps a label to the Assembly
    that would obtain under that history; `prior` gives a declared weight per
    label. The result is the S5.5 mixture with histories in place of modification
    states: E[L], Var, P_competent and the conditional quantities all apply.

    IMPORTANT: `prior` is a *declared scenario weight* expressing a degree of
    belief in, or interest in, an alternative history. It is NOT an evolutionary
    probability the model infers; unlike a modification occupancy it generally
    cannot be measured. The output is therefore a sensitivity distribution over
    what assembly-order space *would* have been admissible, not a prediction of
    which history occurred.
    """
    total = sum(prior.values())
    if abs(total - 1.0) > 1e-6:
        raise ValueError(f"prior weights must sum to 1 (got {total})")
    res = MixtureResult()
    for label, asm in histories.items():
        w = prior.get(label, 0.0)
        if w < 0:
            raise ValueError("prior weights must be non-negative")
        if w == 0.0:
            continue
        competent = (target is None) or bool(target(asm))
        orders = asm.n_orders_permitted() if competent else 0
        res.instances.append(
            {"history": label, "weight": w, "orders": orders,
             "target_competent": competent})
        res.expected_orders += w * orders
    res.exact = True
    return res


def linext_factorizes_disjoint(n_P, L_P, n_Q, L_Q):
    """Exact count for a disjoint union of two independent sub-posets.

    L(P u Q) = C(n_P + n_Q, n_P) * L(P) * L(Q). Counterfactual events confined to
    independent subassemblies therefore act on separate multiplicative factors;
    events at a shared locus are mutually exclusive and combine additively (a
    mixture). Coupled events that share components do not separate and require the
    joint recount.
    """
    from math import comb
    return comb(n_P + n_Q, n_P) * L_P * L_Q
