#!/usr/bin/env python3
"""Tests for the model-level layer (declared alternatives at an occupancy).

Kept separate from verify_paper.py on purpose: that script is the manuscript
regression check and its count (53/53) is quoted in the paper, so it is not
extended here. This file covers the v3.1 additions only.

Reference numbers are the worked SNAP25 phospho-Thr138 example of
Supplementary Section S8.2.
"""
from __future__ import annotations

from sj_prob import (TypedModel, edit_relations, target_orders,
                     occupancy_mixture, seeded_F)
from sj_session import BOARDS
import sj_visual

PASS, FAIL = [], []


def check(name, got, want, tol=1e-9):
    ok = abs(got - want) <= tol if isinstance(want, float) else got == want
    (PASS if ok else FAIL).append(name)
    print(f"  {'PASS' if ok else 'FAIL'}  {name:52s} expected {want!s:>10}  got {got!s:>10}")


def models():
    b = BOARDS["vamp2"]
    m = TypedModel(V=b["V"], C=b["C"], P=b["P_seeded"])
    spec = b["alternatives"]["phospho-SNAP25-Thr138"]
    return b, m, edit_relations(m, add_X=spec["add_X"])


print("\nMODEL-LEVEL LAYER (SNAP25 phospho-Thr138, S8.2)")
print("-" * 78)
b, m, alt = models()
seed = b["seed"]

# --- endpoint models ---------------------------------------------------------
check("factual fusion state, seeded orders", target_orders(m, seed), 252)
check("modified fusion state is infeasible", target_orders(alt, seed), 0)
check("cause is named, not a bare zero",
      alt.blocking_reason(frozenset(m.V)) is not None, True)

# --- the paper's worked occupancy -------------------------------------------
r = occupancy_mixture(m, alt, seed, 0.6)
check("p=0.6  Pr(competent)", round(r["P_competent"], 10), 0.4)
check("p=0.6  E[orders]", round(r["E_orders"], 6), 100.8)
check("p=0.6  Var[orders]", round(r["Var_orders"], 2), 15240.96)
check("p=0.6  E[orders | competent]", round(r["E_orders_given_competent"], 6), 252.0)

# --- the point the slider is meant to make ----------------------------------
conds = {round(occupancy_mixture(m, alt, seed, p)["E_orders_given_competent"], 6)
         for p in (0.1, 0.25, 0.5, 0.75, 0.9)}
check("conditional count invariant in p", conds, {252.0})
check("E[orders] strictly decreasing in p",
      all(occupancy_mixture(m, alt, seed, a)["E_orders"]
          > occupancy_mixture(m, alt, seed, b_)["E_orders"]
          for a, b_ in ((0.1, 0.2), (0.2, 0.5), (0.5, 0.9))), True)

# --- boundary reduction ------------------------------------------------------
check("p=0 reduces to the factual model", occupancy_mixture(m, alt, seed, 0.0)["E_orders"], 252.0)
check("p=0 has zero variance", occupancy_mixture(m, alt, seed, 0.0)["Var_orders"], 0.0)
check("p=1 is the modified model", occupancy_mixture(m, alt, seed, 1.0)["E_orders"], 0.0)
check("p=1 Pr(competent) is zero", occupancy_mixture(m, alt, seed, 1.0)["P_competent"], 0.0)

# --- edit_relations does not disturb the factual model ----------------------
check("edit_relations leaves factual V", alt.V, m.V)
check("edit_relations leaves factual P", alt.P, m.P)
check("factual model still counts 252 after edit", seeded_F(m, seed)(frozenset()), 252)
check("no-op edit is identity on counts", target_orders(edit_relations(m), seed), 252)

# --- occupancy uncertainty: p ~ Beta(alpha, beta)  (S8.2) --------------------
from sj_prob import beta_occupancy_mixture, joint_mixture

rb = beta_occupancy_mixture(m, alt, seed, 3, 2)
check("Beta(3,2)  E[p]", round(rb["E_occupancy"], 6), 0.6)
check("Beta(3,2)  E[orders]", round(rb["E_orders"], 6), 100.8)
check("Beta(3,2)  within-component variance", round(rb["Var_within"], 1), 12700.8)
check("Beta(3,2)  between-component variance", round(rb["Var_between"], 1), 2540.2)
check("Beta(3,2)  95% interval, lower", round(rb["cred_lo"], 1), 17.0)
check("Beta(3,2)  95% interval, upper", round(rb["cred_hi"], 1), 203.1)

# total variance depends on the prior MEAN only, not its spread: the two
# components telescope to L0^2 E[p](1-E[p]). Same mean, very different spread:
tight = beta_occupancy_mixture(m, alt, seed, 60, 40)
check("total variance invariant to prior spread",
      round(tight["Var_total"], 4), round(rb["Var_total"], 4))
check("but the split shifts toward heterogeneity",
      tight["Var_between"] < rb["Var_between"], True)
check("total equals the fixed-p variance at E[p]",
      round(rb["Var_total"], 2), round(252 ** 2 * 0.6 * 0.4, 2))

# --- correlated modifications: declared joint distribution (S8.2) ------------
# (xA,xB) = (0,0),(0,1),(1,0),(1,1) with weights .40,.15,.25,.20; A removes the
# SNAP25 core edge, B does not affect this target. Competent weight = 0.55.
rj = joint_mixture([(0.40, m), (0.15, m), (0.25, alt), (0.20, alt)], seed)
check("correlated pair  P(competent)", round(rj["P_competent"], 6), 0.55)
check("correlated pair  E[orders]", round(rj["E_orders"], 6), 138.6)
check("joint_mixture reduces to the two-model case",
      round(joint_mixture([(0.4, m), (0.6, alt)], seed)["E_orders"], 6),
      round(occupancy_mixture(m, alt, seed, 0.6)["E_orders"], 6))

# --- browser data pipeline ---------------------------------------------------
_, data = sj_visual.build("vamp2")
a = data["alternatives"]["phospho-SNAP25-Thr138"]
check("build() exports factual count", data["L_factual"], 252)
check("build() exports alternative count", a["L"], 0)
check("build() exports the named cause", bool(a["cause"]), True)
check("build() exports declared occupancy", a["occupancy"], 0.6)
check("boards without alternatives still build",
      sj_visual.build("eif3")[1].get("alternatives", {}), {})

print("-" * 78)
print(f"  {len(PASS)}/{len(PASS) + len(FAIL)} checks passed")
if FAIL:
    print("  FAILED: " + ", ".join(FAIL))
    raise SystemExit(1)
