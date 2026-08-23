"""
sj_cryoem.py -- reasoning from a set of structurally resolved intermediates.

Cryo-EM (and any other technique that resolves partial complexes) produces a set
of OBSERVED COMPOSITIONS. Each is a query, never a constraint: the framework
asks what the set as a whole implies about assembly, rather than forcing every
observed species into every pathway.

The module answers five questions, in the order a reader would ask them:

  1. ADMISSIBILITY   Is each observed composition reachable under the declared
                     model? If not, which typed relation blocks it?
  2. CONFLICT        Which PAIRS of observations cannot occur in the same
                     mechanism? This is the direct evidence for branching, and
                     is more informative than any single count.
  3. COVER           What is the minimum number of coexisting mechanisms needed
                     to explain every observation, and which mechanisms are they?
  4. ROBUSTNESS      Does that conclusion survive if any single observation is
                     misassigned or withdrawn? (Leave-one-out.)
  5. IMPLICATION     What precedence relations are entailed, and which are ruled
                     out, by the observation set?

Everything is exact: no sampling, no kinetics, no assumed pathway. The minimum
cover is solved by exhaustive search, which is correct for the small mechanism
spaces these analyses produce; the module refuses rather than approximates if
the space is too large.

Depends on sj_prob.py and sj_navigate.py.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, Iterable, List, Optional, Sequence, Set, Tuple

from sj_prob import TypedModel
from sj_navigate import Navigator

Component = str
Composition = FrozenSet[Component]
History = Tuple[Component, ...]

MAX_MECHANISMS_FOR_EXACT_COVER = 5000


@dataclass(frozen=True)
class Observation:
    """One structurally resolved species."""
    label: str
    composition: Composition
    accession: str = ""
    resolution: Optional[float] = None
    note: str = ""

    def describe(self) -> str:
        res = f", {self.resolution} \u00c5" if self.resolution else ""
        acc = f" ({self.accession}{res})" if self.accession else ""
        return f"{self.label}{acc}"


@dataclass
class CryoEMAnalysis:
    """Assembly reasoning from a set of observed compositions."""
    model: TypedModel
    seed: Component
    observations: Sequence[Observation]
    _nav: Navigator = field(init=False)
    _hist: List[History] = field(init=False)
    _explains: Dict[History, Set[str]] = field(init=False)

    def __post_init__(self) -> None:
        self._nav = Navigator(self.model, self.seed)
        self._hist = [h for h, _ in self._nav.masses]
        if len(self._hist) > MAX_MECHANISMS_FOR_EXACT_COVER:
            raise ValueError(
                f"{len(self._hist)} mechanisms exceeds the exact-cover limit "
                f"({MAX_MECHANISMS_FOR_EXACT_COVER}); constrain the model first")
        self._explains = {h: self._explained_by(h) for h in self._hist}

    # ---------------------------------------------------------------- basics
    def _states(self, h: History) -> List[Composition]:
        """Every intermediate composition passed through by mechanism h."""
        cur, out = {self.seed}, [frozenset({self.seed})]
        for p in h:
            cur = cur | {p}
            out.append(frozenset(cur))
        return out

    def _explained_by(self, h: History) -> Set[str]:
        states = set(self._states(h))
        return {o.label for o in self.observations if o.composition in states}

    @property
    def n_mechanisms(self) -> int:
        return len(self._hist)

    # ------------------------------------------------- 1. admissibility
    def admissibility(self) -> List[Dict[str, object]]:
        rows = []
        for o in self.observations:
            n = sum(1 for h in self._hist if o.label in self._explains[h])
            row: Dict[str, object] = {
                "observation": o.describe(),
                "composition": sorted(o.composition),
                "mechanisms": n,
                "fraction": n / self.n_mechanisms if self.n_mechanisms else 0.0,
            }
            if n == 0:
                r = self._nav.ask(self._nav.forms_exact(o.composition))
                row["verdict"] = "IMPOSSIBLE"
                row["cause"] = r["cause"]
            elif n == self.n_mechanisms:
                row["verdict"] = "necessary"
            else:
                row["verdict"] = "optional"
            rows.append(row)
        return rows

    # ------------------------------------------------------ 2. conflicts
    def conflicts(self) -> List[Tuple[str, str, str]]:
        """Pairs of observations that cannot co-occur in one mechanism, with the
        structural reason. These pairs are the direct evidence of branching."""
        out = []
        for a, b in itertools.combinations(self.observations, 2):
            together = any(a.label in self._explains[h] and b.label in self._explains[h]
                           for h in self._hist)
            if together:
                continue
            reason = self._conflict_reason(a, b)
            out.append((a.label, b.label, reason))
        return out

    def _conflict_reason(self, a: Observation, b: Observation) -> str:
        """Nested compositions are compatible in principle; disjoint additions are
        not, because each exact intermediate fixes what had NOT yet arrived."""
        ca, cb = a.composition, b.composition
        if ca <= cb or cb <= ca:
            return "nested compositions excluded by declared relations"
        only_a = sorted(ca - cb)
        only_b = sorted(cb - ca)
        return (f"{a.label} requires {only_a} present while {only_b} absent; "
                f"{b.label} requires the reverse")

    # ---------------------------------------------------------- 3. cover
    def minimum_cover(self) -> Dict[str, object]:
        """Fewest mechanisms jointly explaining every explainable observation."""
        target = {o.label for o in self.observations
                  if any(o.label in self._explains[h] for h in self._hist)}
        unexplainable = {o.label for o in self.observations} - target
        best_single = max((len(self._explains[h]) for h in self._hist), default=0)
        chosen: Optional[Tuple[History, ...]] = None
        for k in range(1, len(self._hist) + 1):
            for combo in itertools.combinations(self._hist, k):
                if set().union(*(self._explains[c] for c in combo)) >= target:
                    chosen = combo
                    break
            if chosen:
                break
        return {
            "n_required": len(chosen) if chosen else None,
            "mechanisms": [tuple(c) for c in chosen] if chosen else [],
            "covers": [sorted(self._explains[c]) for c in chosen] if chosen else [],
            "best_single": best_single,
            "n_observations": len(target),
            "unexplainable": sorted(unexplainable),
        }

    # ----------------------------------------------------- 4. robustness
    def leave_one_out(self) -> List[Dict[str, object]]:
        """Would the multimodality conclusion survive if one observation were
        misassigned or withdrawn? A conclusion resting on a single 4 A map is
        weaker than one that survives dropping any observation."""
        rows = []
        for o in self.observations:
            sub = CryoEMAnalysis(self.model, self.seed,
                                 [x for x in self.observations if x is not o])
            rows.append({"dropped": o.label,
                         "n_required": sub.minimum_cover()["n_required"]})
        return rows

    # ---------------------------------------------------- 5. implications
    def precedence_implications(self) -> Dict[str, List[str]]:
        """Precedence relations entailed by, or ruled out by, the observations.

        A relation p->q is ENTAILED if it holds in every mechanism that explains
        at least one observation; it is RULED OUT if no such mechanism has it.
        Both are statements about the declared model plus the observation set,
        not about kinetics.
        """
        relevant = [h for h in self._hist if self._explains[h]]
        comps = [c for c in self.model.V if c != self.seed]
        entailed, ruled_out = [], []
        for p, q in itertools.permutations(comps, 2):
            holds = [h.index(p) < h.index(q) for h in relevant
                     if p in h and q in h]
            if holds and all(holds):
                entailed.append(f"{p} -> {q}")
            elif holds and not any(holds):
                ruled_out.append(f"{p} -> {q}")
        return {"entailed": entailed, "ruled_out": ruled_out}

    # ------------------------------------------------------- 6. reporting
    def report(self) -> str:
        L = []
        add = L.append
        add("=" * 76)
        add(f"CRYO-EM ASSEMBLY ANALYSIS   seed={self.seed}   "
            f"{self.n_mechanisms} admissible mechanisms")
        add("=" * 76)

        add("\n1. ADMISSIBILITY OF EACH OBSERVED SPECIES")
        for r in self.admissibility():
            line = (f"   {str(r['observation']):<34} {r['verdict']:<10} "
                    f"{r['mechanisms']}/{self.n_mechanisms} = {r['fraction']:.3f}")
            add(line)
            if r.get("cause"):
                add(f"       cause: {r['cause']}")

        add("\n2. MUTUALLY EXCLUSIVE OBSERVATIONS (evidence for branching)")
        cf = self.conflicts()
        if not cf:
            add("   none: all observations are compatible with a single mechanism")
        for a, b, why in cf:
            add(f"   {a}  vs  {b}")
            add(f"       {why}")

        add("\n3. MINIMUM NUMBER OF COEXISTING MECHANISMS")
        cov = self.minimum_cover()
        add(f"   best single mechanism explains "
            f"{cov['best_single']}/{cov['n_observations']} observations")
        add(f"   MINIMUM REQUIRED: {cov['n_required']}")
        for m, c in zip(cov["mechanisms"], cov["covers"]):
            add(f"     {self.seed}->{'->'.join(m)}   covers {', '.join(c)}")
        if cov["unexplainable"]:
            add(f"   NOT explainable under this model: "
                f"{', '.join(cov['unexplainable'])}")

        add("\n4. ROBUSTNESS (leave one observation out)")
        for r in self.leave_one_out():
            add(f"   without {str(r['dropped']):<26} minimum = {r['n_required']}")
        vals = [r["n_required"] for r in self.leave_one_out()]
        if vals and min(v for v in vals if v) >= 2:
            add("   -> branching conclusion does NOT depend on any single structure")
        else:
            add("   -> branching conclusion depends on at least one structure "
                "(see which line drops to 1)")

        add("\n5. PRECEDENCE IMPLICATIONS OF THE OBSERVATION SET")
        imp = self.precedence_implications()
        add(f"   entailed:  {', '.join(imp['entailed']) or '(none)'}")
        add(f"   ruled out: {', '.join(imp['ruled_out']) or '(none)'}")
        add("\n   Note: these follow from the declared model and the observed")
        add("   compositions. No kinetic or population assumption is used, and")
        add("   class populations are never read as pathway probabilities.")
        return "\n".join(L)


# ------------------------------------------------------------------- demo
if __name__ == "__main__":
    # Mark E, Ramos PC, Nunes MM, Matias AC, Dohmen RJ, Wendler P.
    # Nat Commun 17 (2026). doi:10.1038/s41467-026-70525-w
    SEED = "13S"
    V = (SEED, "b1", "b5", "b6")
    C = frozenset(frozenset({SEED, x}) for x in ("b1", "b5", "b6"))

    OBS = [
        Observation("13S+b1",    frozenset({SEED, "b1"}),             "9rla", 3.16),
        Observation("13S+b5+b6", frozenset({SEED, "b5", "b6"}),       "9rm0", 3.29),
        Observation("13S+b1+b5", frozenset({SEED, "b1", "b5"}),       "9rm1", 4.11),
        Observation("15S",       frozenset({SEED, "b1", "b5", "b6"}), "9rlz", 3.12),
    ]

    print("### A. under the previously assumed order  b5 -> b6 -> b1\n")
    imposed = TypedModel(V=V, C=C, P=(("b5", "b6"), ("b6", "b1")))
    print(CryoEMAnalysis(imposed, SEED, OBS).report())

    print("\n\n### B. with the unsupported precedence removed\n")
    free = TypedModel(V=V, C=C, P=())
    print(CryoEMAnalysis(free, SEED, OBS).report())
