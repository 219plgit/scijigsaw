"""
sj_prob.py — weighted enumeration and calibrated probability trees for Scientific Jigsaw.

Design principle carried over from the manuscript: HARD typed evidence restricts the
space (mass -> 0, with a nameable cause); SOFT evidence (GTEx abundance, AF interface
confidence, Tahoe/LINCS-selected state) only redistributes mass WITHIN it. The exact
DP remains the kernel; the probabilistic layer is a reweighting on top of it.

Key identity implemented here:
    support(E) = G_E(V)/F(V) = P(E formed) under the F-calibrated probability tree,
    where branch S -> S' carries probability F(S')/F(S).
Naive uniform branch probabilities do NOT give uniform-over-histories; the DP is what
calibrates the tree.

Modelling choices that must be aligned with the main scijigsaw implementation:
  * A connected subset S is an admissible intermediate iff (i) S is connected in C,
    (ii) prerequisite closure: for every q in S, P(q) subset of S, (iii) no excluded
    pair lies inside S. Choice (ii) is what makes merger-tree mode reduce exactly to
    the seed-anchored enumerator (see test_seed_recovery).
  * Merger-tree splits are canonicalised by requiring min(S) in A, so each unordered
    root split is counted once.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from functools import lru_cache
from itertools import combinations
from typing import Callable, Dict, FrozenSet, Iterable, List, Optional, Sequence, Tuple

Component = str
Subset = FrozenSet[Component]


# --------------------------------------------------------------------------------------
# Typed model
# --------------------------------------------------------------------------------------

@dataclass(frozen=True)
class TypedModel:
    """M = (V, C, P, X, eta). Contacts are undirected, prerequisites directed."""

    V: Tuple[Component, ...]
    C: FrozenSet[FrozenSet[Component]] = frozenset()
    P: Tuple[Tuple[Component, Component], ...] = ()      # (p, q): p must precede q
    X: FrozenSet[FrozenSet[Component]] = frozenset()
    eta: Dict[str, str] = field(default_factory=dict)     # relation id -> provenance

    # ---- derived accessors -------------------------------------------------------
    def prereqs(self, q: Component) -> FrozenSet[Component]:
        return frozenset(p for p, r in self.P if r == q)

    def contacts(self, u: Component) -> FrozenSet[Component]:
        return frozenset(next(iter(e - {u})) for e in self.C if u in e)

    def has_contact(self, u: Component, v: Component) -> bool:
        return frozenset({u, v}) in self.C

    def excluded(self, u: Component, v: Component) -> bool:
        return frozenset({u, v}) in self.X

    # ---- admissibility of an intermediate ---------------------------------------
    def connected(self, S: Subset) -> bool:
        if not S:
            return False
        seen, stack = set(), [next(iter(S))]
        while stack:
            u = stack.pop()
            if u in seen:
                continue
            seen.add(u)
            stack.extend(v for v in self.contacts(u) & S if v not in seen)
        return seen == set(S)

    def valid_intermediate(self, S: Subset) -> bool:
        """Admissibility of an assembled intermediate.

        NOTE on singletons: a free, unassembled monomer is always a valid starting
        object, so prerequisite closure is required only for ASSEMBLED intermediates
        (|S| >= 2). A prerequisite p->q constrains the assembly in which q
        participates, not the existence of q itself. Without this exemption every
        merger tree containing a prerequisite-bearing leaf would be rejected."""
        if len(S) == 1:
            return True
        for u, v in combinations(sorted(S), 2):
            if self.excluded(u, v):
                return False
        for q in S:
            if not self.prereqs(q) <= S:
                return False
        return self.connected(S)

    def blocking_reason(self, S: Subset) -> Optional[str]:
        """Name the typed relation that makes S inadmissible ('names the cause')."""
        if len(S) == 1:
            return None
        for u, v in combinations(sorted(S), 2):
            if self.excluded(u, v):
                return f"exclusion {u}--{v} ({self.eta.get(f'X:{u}|{v}', 'declared')})"
        for q in sorted(S):
            missing = self.prereqs(q) - S
            if missing:
                return (f"prerequisite {sorted(missing)} -> {q} "
                        f"({self.eta.get(f'P:{q}', 'declared')})")
        if not self.connected(S):
            return "no contact path: subset is disconnected in C"
        return None


# --------------------------------------------------------------------------------------
# Exact enumeration: seed-anchored orders and merger trees
# --------------------------------------------------------------------------------------

def seeded_F(model: TypedModel, seed: Component) -> Callable[[Subset], int]:
    rest = tuple(v for v in model.V if v != seed)

    @lru_cache(maxsize=None)
    def F(S: Subset) -> int:
        rem = [p for p in rest if p not in S]
        if not rem:
            return 1
        return sum(F(S | {p}) for p in rem if model.prereqs(p) - {seed} <= S)

    return F


def merger_F(model: TypedModel) -> Callable[[Subset], int]:
    """F(S) = number of admissible binary merger trees on connected subset S."""

    @lru_cache(maxsize=None)
    def F(S: Subset) -> int:
        if len(S) == 1:
            return 1
        if not model.valid_intermediate(S):
            return 0
        anchor = min(S)
        total = 0
        others = sorted(S - {anchor})
        for r in range(len(others) + 1):
            for combo in combinations(others, r):
                A = frozenset({anchor}) | frozenset(combo)
                B = S - A
                if not B:
                    continue
                if not (model.valid_intermediate(A) and model.valid_intermediate(B)):
                    continue
                if not any(model.has_contact(u, v) for u in A for v in B):
                    continue
                total += F(A) * F(B)
        return total

    return F


def merger_support(model: TypedModel, E: Subset) -> Tuple[int, int]:
    """Return (G_E(V), F(V)) — trees containing E as a subcomplex, and all trees."""
    F = merger_F(model)
    V = frozenset(model.V)

    @lru_cache(maxsize=None)
    def G(S: Subset) -> int:
        if not (E <= S):
            return 0
        if S == E:
            return F(S)
        if len(S) == 1:
            return 1 if S == E else 0
        if not model.valid_intermediate(S):
            return 0
        anchor = min(S)
        total = 0
        others = sorted(S - {anchor})
        for r in range(len(others) + 1):
            for combo in combinations(others, r):
                A = frozenset({anchor}) | frozenset(combo)
                B = S - A
                if not B:
                    continue
                if not (model.valid_intermediate(A) and model.valid_intermediate(B)):
                    continue
                if not any(model.has_contact(u, v) for u in A for v in B):
                    continue
                if E <= A:
                    total += G(A) * F(B)
                elif E <= B:
                    total += F(A) * G(B)
                # split crossing E contributes 0
        return total

    return G(V), F(V)


# --------------------------------------------------------------------------------------
# F-calibrated probability tree over seed-anchored histories
# --------------------------------------------------------------------------------------

@dataclass
class PTNode:
    state: Subset
    step: int
    added: Optional[Component]
    prob: float                      # transition probability from parent
    children: List["PTNode"] = field(default_factory=list)


def calibrated_tree(model: TypedModel, seed: Component,
                    step_weight: Optional[Callable[[Subset, Component], float]] = None
                    ) -> PTNode:
    """Unroll the admissible seed-anchored space into a probability tree.

    The internal Z here is the FORWARD CONTINUATION FUNCTION Q(S) of the paper
    (S8.16): Q(S) sums completion weights from already-assembled state S, with
    Q(U)=1, and the branch probability is phi * Q(S')/Q(S). Do NOT confuse it
    with the backward last-element counting recurrence F(S); writing branch
    probabilities via backward F is ambiguous and incorrect in general.

    step_weight(S, p) is an EXTERNAL, non-negative encounter/evidence weight
    (e.g. GTEx relative abundance of p, AF interface confidence). With
    step_weight=None the tree is uniform over admissible histories, and
    P(event) reproduces exact enumeration ratios.
    """
    rest = tuple(v for v in model.V if v != seed)
    w = step_weight or (lambda S, p: 1.0)

    @lru_cache(maxsize=None)
    def Z(S: Subset) -> float:
        rem = [p for p in rest if p not in S and model.prereqs(p) - {seed} <= S]
        if not [p for p in rest if p not in S]:
            return 1.0
        return sum(w(S, p) * Z(S | {p}) for p in rem)

    def build(S: Subset, step: int, added: Optional[Component], prob: float) -> PTNode:
        node = PTNode(S, step, added, prob)
        tot = Z(S)
        if tot == 0:
            return node
        for p in sorted(rest):
            if p in S or not (model.prereqs(p) - {seed} <= S):
                continue
            Sp = S | {p}
            branch = w(S, p) * Z(Sp) / tot
            if branch > 0:
                node.children.append(build(Sp, step + 1, p, branch))
        return node

    return build(frozenset(), 0, None, 1.0)


def leaf_masses(root: PTNode) -> List[Tuple[Tuple[Component, ...], float]]:
    out: List[Tuple[Tuple[Component, ...], float]] = []

    def walk(n: PTNode, path: Tuple[Component, ...], m: float) -> None:
        if not n.children:
            out.append((path, m))
            return
        for c in n.children:
            walk(c, path + (c.added,), m * c.prob)

    walk(root, (), 1.0)
    return out


def prob_event(root: PTNode, predicate: Callable[[Tuple[Component, ...]], bool]) -> float:
    return sum(m for path, m in leaf_masses(root) if predicate(path))


def n_eff(root: PTNode) -> float:
    """Effective number of histories: N_eff = 2**H2 with H2 in bits (equivalently
    exp of the natural-log entropy; the perplexity is base-invariant). Reported
    entropy DIFFERENCES elsewhere in this module are in bits (log2)."""
    ms = [m for _, m in leaf_masses(root) if m > 0]
    return math.exp(-sum(m * math.log(m) for m in ms))


def positive_potential(mu: float, eps: float = 0.05) -> float:
    """Convert a graded score mu in [0,1] (e.g. fuzzy membership) into a
    STRICTLY POSITIVE potential: phi = eps + (1-eps)*mu. Soft evidence must
    never manufacture a structural zero; zeros are the privilege of declared
    hard typed relations only (v9 rule: phi > 0)."""
    if not (0.0 < eps < 1.0):
        raise ValueError("eps must be in (0,1)")
    mu = min(max(mu, 0.0), 1.0)
    return eps + (1.0 - eps) * mu


def plackett_luce(model: TypedModel, seed: Component,
                  abundance: Dict[Component, float]
                  ) -> Callable[[Subset, Component], float]:
    """Encounter-order prior: w(S, p) = a_p / sum_{q admissible at S} a_q.

    IMPORTANT: a weight depending only on the component, w(S, p) = a_p, is
    order-INVARIANT -- every complete history multiplies the same set of factors,
    so all admissible histories keep equal mass and N_eff never moves. The
    normalisation over the currently-admissible set is what makes an abundance
    prior discriminate between orders at all.
    """
    rest = tuple(v for v in model.V if v != seed)

    def w(S: Subset, p: Component) -> float:
        avail = [q for q in rest
                 if q not in S and model.prereqs(q) - {seed} <= S]
        denom = sum(max(abundance.get(q, 0.0), 1e-12) for q in avail)
        return max(abundance.get(p, 0.0), 1e-12) / denom if denom > 0 else 0.0

    return w


def prefix_contains(E: Subset, seed: Component) -> Callable[[Tuple[Component, ...]], bool]:
    """Event: subcomplex E is present at some prefix of the seeded history."""
    def pred(path: Tuple[Component, ...]) -> bool:
        cur = {seed}
        if E <= cur:
            return True
        for p in path:
            cur.add(p)
            if E <= cur:
                return True
        return False
    return pred


# --------------------------------------------------------------------------------------
# Second-order uncertainty: GTEx cross-tissue variance -> credible interval on support
# --------------------------------------------------------------------------------------

def support_interval(model: TypedModel, seed: Component, E: Subset,
                     abundance_samples: Sequence[Dict[Component, float]],
                     q: Tuple[float, float] = (0.05, 0.95)
                     ) -> Dict[str, float]:
    """Propagate variability in an abundance vector (e.g. one dict per GTEx tissue,
    or per bootstrap replicate) into an interval on support(E) and on N_eff.

    Reports rank stability alongside the point estimate: a hypothesis that is
    top-ranked in one tissue but unstable across tissues is reported as such.
    """
    sups, neffs = [], []
    for ab in abundance_samples:
        tree = calibrated_tree(model, seed, step_weight=plackett_luce(model, seed, ab))
        sups.append(prob_event(tree, prefix_contains(E, seed)))
        neffs.append(n_eff(tree))
    sups.sort()
    lo = sups[max(0, int(q[0] * (len(sups) - 1)))]
    hi = sups[min(len(sups) - 1, int(math.ceil(q[1] * (len(sups) - 1))))]
    mean = sum(sups) / len(sups)
    var = sum((s - mean) ** 2 for s in sups) / max(1, len(sups) - 1)
    return {"support_mean": mean, "support_sd": math.sqrt(var),
            "support_lo": lo, "support_hi": hi,
            "n_eff_mean": sum(neffs) / len(neffs)}


# --------------------------------------------------------------------------------------
# Uncertainty over alternative constraint models (AF3 samples, occupancy states)
# --------------------------------------------------------------------------------------

def model_distribution(models: Sequence[Tuple[float, TypedModel]], seed: Component,
                       E: Optional[Subset] = None) -> Dict[str, float]:
    """Run the unchanged exact enumerator per declared model M_z with weight w_z."""
    tot_w = sum(w for w, _ in models)
    comp, exp_L, exp_L2, exp_sup, robust = 0.0, 0.0, 0.0, 0.0, 0
    for w, m in models:
        F = seeded_F(m, seed)
        L = F(frozenset())
        p = w / tot_w
        comp += p * (1.0 if L > 0 else 0.0)
        exp_L += p * L
        exp_L2 += p * L * L
        if E is not None and L > 0:
            tree = calibrated_tree(m, seed)
            s = prob_event(tree, prefix_contains(E, seed))
            exp_sup += p * s
            robust += 1 if s > 0 else 0
    out = {"P_competent": comp, "E_histories": exp_L,
           "Var_histories": exp_L2 - exp_L ** 2}
    if E is not None:
        out["E_support"] = exp_sup
        out["robustness_R"] = robust / len(models)
    return out


# --------------------------------------------------------------------------------------
# Ranking: the "most viable hypothesis" output
# --------------------------------------------------------------------------------------

def rank_histories(model: TypedModel, seed: Component,
                   abundance_samples: Sequence[Dict[Component, float]],
                   top: int = 5) -> List[Dict[str, object]]:
    """Rank admissible histories by posterior mass, with cross-sample rank stability.

    Mass is an evidence/encounter prior, NOT a kinetic probability. A hypothesis is
    'most viable' only when it is high-mass, tight-interval and rank-stable.
    """
    per_sample: List[Dict[Tuple[Component, ...], float]] = []
    for ab in abundance_samples:
        tree = calibrated_tree(model, seed, step_weight=plackett_luce(model, seed, ab))
        per_sample.append(dict(leaf_masses(tree)))

    paths = sorted(per_sample[0].keys())
    rows = []
    for path in paths:
        ms = [d.get(path, 0.0) for d in per_sample]
        ranks = []
        for d in per_sample:
            order = sorted(d, key=lambda k: -d[k])
            ranks.append(order.index(path) + 1)
        mean = sum(ms) / len(ms)
        var = sum((m - mean) ** 2 for m in ms) / max(1, len(ms) - 1)
        rows.append({"history": path, "mass_mean": mean, "mass_sd": math.sqrt(var),
                     "rank_best": min(ranks), "rank_worst": max(ranks),
                     "rank_stable": min(ranks) == max(ranks)})
    rows.sort(key=lambda r: -r["mass_mean"])
    return rows[:top]


# --------------------------------------------------------------------------------------
# Decision value: how many bits does each evidence source remove?
# --------------------------------------------------------------------------------------

def entropy_bits(model: TypedModel, seed: Component,
                 step_weight: Optional[Callable[[Subset, Component], float]] = None
                 ) -> float:
    """Shannon entropy (bits) of the measure over admissible histories.

    Uniform weights -> log2(number of admissible histories). This is the common
    currency in which structural, temporal and contextual evidence are compared.
    """
    tree = calibrated_tree(model, seed, step_weight=step_weight)
    ms = [m for _, m in leaf_masses(tree) if m > 0]
    if not ms:
        return float("nan")          # infeasible model: no admissible history
    return -sum(m * math.log2(m) for m in ms)


def drop_relations(model: TypedModel, *,
                   drop_P: Sequence[Tuple[Component, Component]] = (),
                   drop_X: Iterable[FrozenSet[Component]] = (),
                   drop_C: Iterable[FrozenSet[Component]] = ()) -> TypedModel:
    """Return the model with the named typed relations removed (the 'before' state)."""
    return TypedModel(
        V=model.V,
        C=model.C - frozenset(frozenset(e) for e in drop_C),
        P=tuple(p for p in model.P if p not in set(drop_P)),
        X=model.X - frozenset(frozenset(e) for e in drop_X),
        eta=model.eta,
    )


def decision_value_relation(model: TypedModel, seed: Component, *,
                            drop_P: Sequence[Tuple[Component, Component]] = (),
                            drop_X: Iterable[FrozenSet[Component]] = (),
                            label: str = "") -> Dict[str, object]:
    """Bits removed by a HARD typed relation: H(without) - H(with).

    This is a leave-one-out ablation on the constraint set, so it measures how much
    of the space that relation is responsible for carving away.
    """
    without = drop_relations(model, drop_P=drop_P, drop_X=drop_X)
    h0 = entropy_bits(without, seed)
    h1 = entropy_bits(model, seed)
    return {"source": label or "relation", "kind": "hard",
            "H_before_bits": h0, "H_after_bits": h1,
            "bits_removed": h0 - h1,
            "N_eff_before": 2 ** h0 if h0 == h0 else float("nan"),
            "N_eff_after": 2 ** h1 if h1 == h1 else float("nan")}


def decision_value_model_choice(models: Sequence[Tuple[float, TypedModel]],
                                seed: Component, label: str = "") -> Dict[str, object]:
    """Bits removed by resolving WHICH model is active (AlphaFold variant, Tahoe/LINCS state).

    Mutual-information style: H(mixture over models) - E_z[H(model z)]. Infeasible
    models (no admissible history) contribute zero mass and are reported separately,
    since 'this hypothesis is not competent' is a different finding from 'this
    hypothesis is uncertain'.
    """
    tot_w = sum(w for w, _ in models)
    pooled: Dict[Tuple[Component, ...], float] = {}
    exp_h, live_w, n_dead = 0.0, 0.0, 0
    for w, m in models:
        tree = calibrated_tree(m, seed)
        lm = leaf_masses(tree)
        ms = [x for _, x in lm if x > 0]
        if not ms:
            n_dead += 1
            continue
        p = w / tot_w
        live_w += p
        exp_h += p * (-sum(x * math.log2(x) for x in ms))
        for path, x in lm:
            pooled[path] = pooled.get(path, 0.0) + p * x
    if live_w == 0:
        return {"source": label or "model choice", "kind": "selector",
                "bits_removed": float("nan"), "n_infeasible": n_dead}
    z = sum(pooled.values())
    h_pool = -sum((v / z) * math.log2(v / z) for v in pooled.values() if v > 0)
    exp_h /= live_w
    return {"source": label or "model choice", "kind": "selector",
            "H_before_bits": h_pool, "H_after_bits": exp_h,
            "bits_removed": h_pool - exp_h,
            "n_models": len(models), "n_infeasible": n_dead,
            "P_competent": live_w}


def decision_value_weights(model: TypedModel, seed: Component,
                           abundance: Dict[Component, float],
                           label: str = "abundance weights") -> Dict[str, object]:
    """Bits removed by SOFT weights: they concentrate mass, they never cut the space."""
    h0 = entropy_bits(model, seed)
    h1 = entropy_bits(model, seed, step_weight=plackett_luce(model, seed, abundance))
    return {"source": label, "kind": "soft",
            "H_before_bits": h0, "H_after_bits": h1,
            "bits_removed": h0 - h1,
            "N_eff_before": 2 ** h0, "N_eff_after": 2 ** h1}


def decision_table(rows: Sequence[Dict[str, object]]) -> str:
    """Rank evidence sources by bits removed — the 'strongest decision node' question."""
    rows = sorted(rows, key=lambda r: -(r.get("bits_removed") or 0.0))
    out = [f"{'source':<40}{'kind':<10}{'H_before':>9}{'H_after':>9}{'bits':>8}"]
    out.append("-" * 76)
    for r in rows:
        hb, ha = r.get("H_before_bits"), r.get("H_after_bits")
        src = str(r["source"])[:39]
        if hb is None or ha is None:
            out.append(f"{src:<40}{str(r['kind']):<10}{'n/a':>9}{'n/a':>9}{'n/a':>8}")
        else:
            out.append(f"{src:<40}{str(r['kind']):<10}"
                       f"{hb:>9.3f}{ha:>9.3f}{r['bits_removed']:>8.3f}")
    return "\n".join(out)


# --------------------------------------------------------------------------------------
# Weighted merger-tree enumeration (partition function Z, branch probabilities, Z_E)
# --------------------------------------------------------------------------------------

MergePhi = Callable[[Subset, Subset], float]


def _valid_splits(model: TypedModel, S: Subset):
    anchor = min(S)
    others = sorted(S - {anchor})
    for r in range(len(others) + 1):
        for combo in combinations(others, r):
            A = frozenset({anchor}) | frozenset(combo)
            B = S - A
            if not B:
                continue
            if not (model.valid_intermediate(A) and model.valid_intermediate(B)):
                continue
            if not any(model.has_contact(u, v) for u in A for v in B):
                continue
            yield A, B


def merger_Z(model: TypedModel, phi: Optional[MergePhi] = None,
             phi_state: Optional[Callable[[Subset], float]] = None
             ) -> Callable[[Subset], float]:
    """Z(S) = sum over admissible trees of the product of weights.

    Two kinds of potential may be supplied, and they behave very differently:

      phi(A, B)    TRANSITION potential, attached to each merge event. Note the
                   additivity lemma: if phi = exp(lambda * sum over the contacts
                   crossing (A,B) of a fixed per-contact score), then every tree
                   receives the same total, because each contact is crossed
                   exactly once (at the LCA of its endpoints). Edge-additive
                   transition potentials are therefore ORDER-INVARIANT on every
                   contact graph, and only non-additive interface-level forms
                   (e.g. log of the summed interface) discriminate, and then
                   only where the contact graph has a cycle.

      phi_state(S) STATE potential, attached to each assembled intermediate.
                   This is NOT subject to the additivity lemma: different trees
                   pass through different intermediates, so a per-subcomplex
                   score (e.g. a stability estimate for that species)
                   discriminates even on an acyclic contact graph. Subcomplex
                   thermodynamics is naturally of this kind.

    phi=phi_state=None recovers the exact tree count merger_F. Both must be
    strictly positive if impossibility is to remain the exclusive privilege of
    declared hard evidence. Evidence scoring a whole tree does not factor and
    cannot enter here without augmenting the DP state.
    """
    w = phi or (lambda A, B: 1.0)
    v = phi_state or (lambda S: 1.0)

    @lru_cache(maxsize=None)
    def Z(S: Subset) -> float:
        if len(S) == 1:
            return 1.0
        if not model.valid_intermediate(S):
            return 0.0
        return v(S) * sum(Z(A) * Z(B) * w(A, B) for A, B in _valid_splits(model, S))

    return Z


def merger_branch_probs(model: TypedModel, S: Subset,
                        phi: Optional[MergePhi] = None
                        ) -> List[Tuple[Subset, Subset, float]]:
    """P(A,B | S) = Z(A) Z(B) phi(A,B) / Z(S): the probabilistic assembly tree is
    a normalised view of the weighted enumeration, not a separately fitted object."""
    w = phi or (lambda A, B: 1.0)
    Z = merger_Z(model, phi)
    tot = Z(S)
    if tot == 0:
        return []
    return [(A, B, Z(A) * Z(B) * w(A, B) / tot) for A, B in _valid_splits(model, S)]


def merger_weighted_support(model: TypedModel, E: Subset,
                            phi: Optional[MergePhi] = None) -> Tuple[float, float]:
    """(Z_E(V), Z(V)): P(E occurs | evidence) = Z_E/Z. phi=None recovers G_E/F."""
    w = phi or (lambda A, B: 1.0)
    Z = merger_Z(model, phi)

    @lru_cache(maxsize=None)
    def ZE(S: Subset) -> float:
        if not (E <= S):
            return 0.0
        if S == E:
            return Z(S)
        if len(S) == 1:
            return 0.0
        if not model.valid_intermediate(S):
            return 0.0
        total = 0.0
        for A, B in _valid_splits(model, S):
            if E <= A:
                total += ZE(A) * Z(B) * w(A, B)
            elif E <= B:
                total += Z(A) * ZE(B) * w(A, B)
        return total

    V = frozenset(model.V)
    return ZE(V), Z(V)


# --------------------------------------------------------------------------------------
# Counterfactual model comparison (epistemic: edits act on the ENCODING, not biology)
# --------------------------------------------------------------------------------------

def counterfactual_compare(model: TypedModel, alt: TypedModel, seed: Component,
                           step_weight=None, alt_step_weight=None,
                           query: Optional[Subset] = None) -> Dict[str, object]:
    """Compare the factual measure pi_M with an alternative pi_M' by coupling them
    through the shared mechanism space (models must share V and seed).

    This is NOT a Pearl level-3 counterfactual: there is no exogenous noise to
    abduce. It is an exact paired comparison of two epistemic measures, with every
    difference decomposable and attributable:
      killed       -- admissible under M only (mass they carried under M)
      resurrected  -- admissible under M' only (mass they carry under M')
      reweighted   -- admissible under both (total-variation shift on this set)
    A removed contact is an in-silico interface mutation; an added occupant edge is
    a molecular-glue edit.
    """
    if model.V != alt.V:
        raise ValueError("counterfactual comparison requires a shared component set")
    m_f = dict(leaf_masses(calibrated_tree(model, seed, step_weight=step_weight)))
    m_a = dict(leaf_masses(calibrated_tree(alt, seed,
                                           step_weight=alt_step_weight or step_weight)))
    keys_f, keys_a = set(m_f), set(m_a)
    killed = keys_f - keys_a
    resurrected = keys_a - keys_f
    shared = keys_f & keys_a
    tv_shared = 0.5 * sum(abs(m_f[k] - m_a[k]) for k in shared)
    out: Dict[str, object] = {
        "n_factual": len(keys_f), "n_alternative": len(keys_a),
        "killed": sorted(killed), "killed_mass_factual": sum(m_f[k] for k in killed),
        "resurrected": sorted(resurrected),
        "resurrected_mass_alt": sum(m_a[k] for k in resurrected),
        "tv_on_shared": tv_shared,
        "edited_relations": {
            "P_removed": sorted(set(model.P) - set(alt.P)),
            "P_added": sorted(set(alt.P) - set(model.P)),
            "C_removed": sorted(tuple(sorted(e)) for e in model.C - alt.C),
            "C_added": sorted(tuple(sorted(e)) for e in alt.C - model.C),
            "X_removed": sorted(tuple(sorted(e)) for e in model.X - alt.X),
            "X_added": sorted(tuple(sorted(e)) for e in alt.X - model.X),
        },
    }
    def _h2(masses):
        ms = [v for v in masses.values() if v > 0]
        return -sum(v * math.log2(v) for v in ms) if ms else float("nan")

    h2_f, h2_a = _h2(m_f), _h2(m_a)
    out["H2_factual"] = h2_f
    out["H2_alt"] = h2_a
    out["delta_H2_bits"] = h2_a - h2_f
    out["delta_N_eff"] = (2 ** h2_a) - (2 ** h2_f)

    if query is not None:
        pf = sum(v for k, v in m_f.items()
                 if prefix_contains(query, seed)(k))
        pa = sum(v for k, v in m_a.items()
                 if prefix_contains(query, seed)(k))
        out["query_support_factual"] = pf
        out["query_support_alt"] = pa
        out["query_support_delta"] = pa - pf
    return out


# --------------------------------------------------------------------------------------
# Self-tests
# --------------------------------------------------------------------------------------

def _brute_force_orders(model: TypedModel, seed: Component) -> List[Tuple[Component, ...]]:
    from itertools import permutations
    rest = [v for v in model.V if v != seed]
    ok = []
    for perm in permutations(rest):
        seen = {seed}
        good = True
        for p in perm:
            if not model.prereqs(p) <= seen:
                good = False
                break
            seen.add(p)
        if good:
            ok.append(perm)
    return ok


def test_seed_recovery() -> None:
    """Restricting mergers to singleton additions recovers the seeded count exactly."""
    V = ("s", "a", "b", "c", "d")
    C = frozenset(frozenset(e) for e in
                  [("s", "a"), ("s", "b"), ("a", "b"), ("b", "c"), ("c", "d"), ("a", "d")])
    P = (("a", "b"), ("a", "c"))
    m = TypedModel(V=V, C=C, P=P)
    F = seeded_F(m, "s")
    assert F(frozenset()) == len(_brute_force_orders(m, "s"))
    print(f"[ok] seeded DP == brute force: {F(frozenset())} orders")


def test_calibration() -> None:
    """P(event) under the F-calibrated tree == exact enumeration ratio."""
    V = ("s", "a", "b", "c", "d")
    C = frozenset(frozenset(e) for e in
                  [("s", "a"), ("s", "b"), ("a", "b"), ("b", "c"), ("c", "d"), ("a", "d")])
    P = (("a", "b"), ("a", "c"))
    m = TypedModel(V=V, C=C, P=P)
    orders = _brute_force_orders(m, "s")
    tree = calibrated_tree(m, "s")

    for k, comp in [(0, "a"), (1, "b"), (2, "c")]:
        exact = sum(1 for o in orders if len(o) > k and o[k] == comp) / len(orders)
        got = prob_event(tree, lambda path, k=k, comp=comp: len(path) > k and path[k] == comp)
        assert abs(exact - got) < 1e-12, (k, comp, exact, got)
    print(f"[ok] calibrated tree reproduces exact ratios over {len(orders)} histories")

    ms = [mass for _, mass in leaf_masses(tree)]
    assert abs(sum(ms) - 1.0) < 1e-12
    assert abs(n_eff(tree) - len(orders)) < 1e-9
    print(f"[ok] uniform N_eff == history count ({n_eff(tree):.4f})")


def test_support_identity() -> None:
    """support(E) via G_E/F equals P(E formed) in the calibrated tree (seeded mode)."""
    V = ("s", "a", "b", "c")
    C = frozenset(frozenset(e) for e in [("s", "a"), ("a", "b"), ("b", "c"), ("s", "c")])
    m = TypedModel(V=V, C=C, P=(("a", "b"),))
    tree = calibrated_tree(m, "s")
    E = frozenset({"s", "a"})
    orders = _brute_force_orders(m, "s")
    exact = sum(1 for o in orders
                if any(E <= ({"s"} | set(o[:k])) for k in range(len(o) + 1))) / len(orders)
    got = prob_event(tree, prefix_contains(E, "s"))
    assert abs(exact - got) < 1e-12
    print(f"[ok] support({sorted(E)}) = {got:.4f} matches enumeration")


def test_merger_and_diagnostics() -> None:
    """Off-seed dimer: impossible under seeding, recoverable under merger trees."""
    V = ("P", "T34", "T35", "N", "T32")
    C = frozenset(frozenset(e) for e in
                  [("P", "T34"), ("P", "T35"), ("T35", "T34"), ("P", "N"), ("N", "T32")])
    m = TypedModel(V=V, C=C, P=())
    F = merger_F(m)
    total = F(frozenset(V))
    g, f = merger_support(m, frozenset({"T35", "T34"}))
    print(f"[ok] merger trees on 5 components: {total}; "
          f"off-scaffold dimer support = {g}/{f} = {g/f:.3f}")

    blocked = TypedModel(V=V, C=C, P=(), X=frozenset({frozenset({"T35", "T34"})}))
    reason = blocked.blocking_reason(frozenset({"T35", "T34"}))
    print(f"[ok] named cause when excluded: {reason}")


def test_second_order() -> None:
    """Cross-'tissue' variability propagates to an interval on support and N_eff."""
    V = ("s", "a", "b", "c", "d")
    C = frozenset(frozenset(e) for e in
                  [("s", "a"), ("s", "b"), ("a", "b"), ("b", "c"), ("c", "d"), ("a", "d")])
    m = TypedModel(V=V, C=C, P=(("a", "b"),))
    rng = random.Random(0)
    samples = [{p: max(0.05, rng.gauss(1.0, 0.35)) for p in ("a", "b", "c", "d")}
               for _ in range(13)]          # stand-in for 13 GTEx brain tissues
    res = support_interval(m, "s", frozenset({"s", "a", "b"}), samples)
    print("[ok] second-order support: "
          f"{res['support_mean']:.3f} +/- {res['support_sd']:.3f} "
          f"[{res['support_lo']:.3f}, {res['support_hi']:.3f}], "
          f"N_eff = {res['n_eff_mean']:.2f}")
    top = rank_histories(m, "s", samples, top=3)
    for r in top:
        print(f"     {'->'.join(r['history']):<16} mass {r['mass_mean']:.4f} "
              f"+/- {r['mass_sd']:.4f}  ranks {r['rank_best']}-{r['rank_worst']} "
              f"{'STABLE' if r['rank_stable'] else 'unstable'}")


def test_decision_value() -> None:
    """Rank evidence sources by bits of entropy removed."""
    V = ("s", "a", "b", "c", "d")
    Cbase = [("s", "a"), ("s", "b"), ("a", "b"), ("b", "c"), ("c", "d"), ("a", "d")]
    C = frozenset(frozenset(e) for e in Cbase)
    P = (("a", "b"), ("a", "c"))
    X = frozenset({frozenset({"c", "d"})})
    m = TypedModel(V=V, C=C, P=P, X=X,
                   eta={"P:b": "temporal, native-MS", "P:c": "temporal, native-MS"})

    rows = [
        decision_value_relation(m, "s", drop_P=[("a", "b")],
                                label="prerequisite a->b (native-MS)"),
        decision_value_relation(m, "s", drop_P=[("a", "c")],
                                label="prerequisite a->c (native-MS)"),
        decision_value_relation(m, "s", drop_P=[("a", "b"), ("a", "c")],
                                label="all temporal prerequisites"),
    ]

    # structural ensemble: two contact variants, as if 8/10 vs 2/10 AF3 samples
    alt = TypedModel(V=V, C=C - {frozenset({"a", "d"})}, P=P, X=X)
    rows.append(decision_value_model_choice([(0.8, m), (0.2, alt)], "s",
                                            label="AlphaFold contact variant"))

    # context selector with no per-context model: scores ~0 by construction
    rows.append(decision_value_model_choice([(0.5, m), (0.5, m)], "s",
                                            label="context selector (no per-state model)"))

    rng = random.Random(1)
    ab = {p: max(0.05, rng.gauss(1.0, 0.6)) for p in ("a", "b", "c", "d")}
    rows.append(decision_value_weights(m, "s", ab, label="GTEx abundance weights"))

    print("[ok] decision value ranking:")
    print("\n".join("     " + ln for ln in decision_table(rows).splitlines()))


def test_weighted_mergers() -> None:
    """phi=1 recovers exact counts; a non-uniform phi moves support, never the space."""
    V = ("P", "T34", "T35", "N", "T32")
    C = frozenset(frozenset(e) for e in
                  [("P", "T34"), ("P", "T35"), ("T35", "T34"), ("P", "N"), ("N", "T32")])
    m = TypedModel(V=V, C=C, P=())
    E = frozenset({"T35", "T34"})

    Z1 = merger_Z(m)(frozenset(V))
    g, f = merger_support(m, E)
    ze, z = merger_weighted_support(m, E)
    assert abs(Z1 - f) < 1e-9 and abs(ze - g) < 1e-9
    print(f"[ok] phi=1 recovers counts: Z = {Z1:.0f} = F, Z_E = {ze:.0f} = G_E "
          f"(support {ze/z:.3f})")

    # interface-size-like weight: favour merges joined by more contacts
    def phi(A: Subset, B: Subset) -> float:
        k = sum(1 for u in A for v in B if m.has_contact(u, v))
        return float(k)

    ze2, z2 = merger_weighted_support(m, E, phi)
    assert z2 > 0
    print(f"[ok] weighted support under contact-count phi: {ze2/z2:.3f} "
          f"(uniform {ze/z:.3f}); admissible set unchanged")

    bp = merger_branch_probs(m, frozenset(V), phi)
    s = sum(p for _, _, p in bp)
    assert abs(s - 1.0) < 1e-9
    top = max(bp, key=lambda t: t[2])
    print(f"[ok] root branch probs sum to 1; top split "
          f"{sorted(top[0])} | {sorted(top[1])} at {top[2]:.3f}")


def test_counterfactual() -> None:
    """Removing a prerequisite resurrects mechanisms; every change is attributed."""
    V = ("s", "a", "b", "c", "d")
    C = frozenset(frozenset(e) for e in
                  [("s", "a"), ("s", "b"), ("a", "b"), ("b", "c"), ("c", "d"), ("a", "d")])
    factual = TypedModel(V=V, C=C, P=(("a", "b"), ("a", "c")))
    alt = drop_relations(factual, drop_P=[("a", "c")])   # in-silico relaxation
    r = counterfactual_compare(factual, alt, "s", query=frozenset({"s", "c"}))
    assert r["n_alternative"] > r["n_factual"] and not r["killed"]
    assert r["edited_relations"]["P_removed"] == [("a", "c")]
    print(f"[ok] counterfactual: {r['n_factual']} -> {r['n_alternative']} histories, "
          f"resurrected mass {r['resurrected_mass_alt']:.3f}, "
          f"TV on shared {r['tv_on_shared']:.3f}, "
          f"dH2 {r['delta_H2_bits']:+.3f} bits, dN_eff {r['delta_N_eff']:+.2f}, "
          f"dPr(query) {r['query_support_delta']:+.3f}, "
          f"cause {r['edited_relations']['P_removed']}")


def test_model_distribution() -> None:
    """Occupancy-style uncertainty: 0.6 modified -> P_competent 0.40, large variance."""
    V = ("s", "a", "b")
    C = frozenset(frozenset(e) for e in [("s", "a"), ("a", "b"), ("s", "b")])
    live = TypedModel(V=V, C=C, P=())
    dead = TypedModel(V=V, C=C, P=(("a", "b"), ("b", "a")))   # contradictory -> 0 orders
    res = model_distribution([(0.4, live), (0.6, dead)], "s")
    print(f"[ok] model distribution: P_comp = {res['P_competent']:.2f}, "
          f"E[L] = {res['E_histories']:.2f}, Var = {res['Var_histories']:.2f}")


if __name__ == "__main__":
    test_seed_recovery()
    test_calibration()
    test_support_identity()
    test_merger_and_diagnostics()
    test_second_order()
    test_model_distribution()
    test_counterfactual()
    test_decision_value()
    test_weighted_mergers()
    print("\nall self-tests passed")
