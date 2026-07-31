"""Tests for the probabilistic modification layer (Supplementary S5.4)."""
import math
import pytest
import scijigsaw.cases as C
from scijigsaw.contrib.probabilistic import Modification, order_distribution


def _block_snap25(requires, excludes):
    for k in list(requires):
        requires[k] = [u for u in requires[k] if u != "SNAP25"]
    requires["SNAP25"] = ["__blocked__"]


def _fusion_competent(asm):
    return "__blocked__" not in asm.requires.get("SNAP25", [])


def test_deterministic_limits_recovered():
    # p=0 -> unmodified fusion count; p=1 -> fusion infeasible (0)
    m0 = Modification("x", 0.0, _block_snap25)
    m1 = Modification("x", 1.0, _block_snap25)
    r0 = order_distribution(C.VAMP2, [m0], target=_fusion_competent)
    r1 = order_distribution(C.VAMP2, [m1], target=_fusion_competent)
    assert r0.expected_orders == 252
    assert r1.expected_orders == 0


def test_expected_orders_linear_in_occupancy():
    for p in (0.0, 0.3, 0.6, 0.9, 1.0):
        m = Modification("x", p, _block_snap25)
        r = order_distribution(C.VAMP2, [m], target=_fusion_competent)
        assert math.isclose(r.expected_orders, (1 - p) * 252, rel_tol=1e-9)
        assert math.isclose(r.probability_nonzero(), 1 - p, rel_tol=1e-9)


def test_weights_sum_to_one():
    m = Modification("x", 0.6, _block_snap25)
    r = order_distribution(C.VAMP2, [m], target=_fusion_competent)
    assert math.isclose(sum(i["weight"] for i in r.instances), 1.0, rel_tol=1e-9)


def test_two_independent_mods_enumerate_four_instances():
    m1 = Modification("a", 0.5, _block_snap25)
    m2 = Modification("b", 0.5, lambda r, e: None)  # no-op effect
    r = order_distribution(C.VAMP2, [m1, m2], target=_fusion_competent)
    assert len(r.instances) == 4
    assert r.exact is True


def test_occupancy_out_of_range_rejected():
    with pytest.raises(ValueError):
        Modification("bad", 1.5, _block_snap25)


def test_sampling_branch_approximates_exact():
    # force sampling with max_exact=0; mean should be near the exact 100.8
    m = Modification("x", 0.6, _block_snap25)
    r = order_distribution(C.VAMP2, [m], target=_fusion_competent,
                           max_exact=0, n_samples=8000, rng_seed=1)
    assert r.exact is False
    assert abs(r.expected_orders - 100.8) < 12  # within sampling error


def test_summary_quantities_match_closed_form():
    # at p=0.6 the four quantities have exact closed forms
    m = Modification("x", 0.6, _block_snap25)
    r = order_distribution(C.VAMP2, [m], target=_fusion_competent)
    s = r.summary(total_permutations=5040)
    assert math.isclose(s["P_competent"], 0.4, rel_tol=1e-9)
    assert math.isclose(s["E_orders"], 100.8, rel_tol=1e-9)
    assert math.isclose(s["E_fraction"], 100.8 / 5040, rel_tol=1e-9)
    # Var(L) = 252^2 * p(1-p)
    assert math.isclose(s["Var_orders"], 252**2 * 0.6 * 0.4, rel_tol=1e-9)
    # E[L | competent] = 252 whenever competent
    assert math.isclose(s["E_orders_given_competent"], 252.0, rel_tol=1e-9)


def test_variance_zero_at_deterministic_limits():
    for p in (0.0, 1.0):
        m = Modification("x", p, _block_snap25)
        r = order_distribution(C.VAMP2, [m], target=_fusion_competent)
        assert math.isclose(r.variance(), 0.0, abs_tol=1e-9)


def test_joint_distribution_matches_marginal_weight():
    from scijigsaw.contrib.probabilistic import order_distribution_joint
    mA = Modification("A", 0.3, _block_snap25)
    mB = Modification("B", 0.45, lambda r, e: None)
    joint = {(0, 0): 0.40, (0, 1): 0.15, (1, 0): 0.25, (1, 1): 0.20}
    r = order_distribution_joint(C.VAMP2, [mA, mB], joint, target=_fusion_competent)
    # fusion-competent only when A absent: weight 0.40+0.15 = 0.55
    assert math.isclose(r.probability_nonzero(), 0.55, rel_tol=1e-9)
    assert math.isclose(r.expected_orders, 0.55 * 252, rel_tol=1e-9)


def test_joint_weights_must_sum_to_one():
    from scijigsaw.contrib.probabilistic import order_distribution_joint
    mA = Modification("A", 0.3, _block_snap25)
    with pytest.raises(ValueError):
        order_distribution_joint(C.VAMP2, [mA], {(0,): 0.5, (1,): 0.2},
                                 target=_fusion_competent)


def test_beta_uncertainty_analytic_exact():
    from scijigsaw.contrib.probabilistic import occupancy_uncertainty_analytic
    # Beta(3,2) mean 0.6 -> E[L] = 252*0.4 = 100.8 exactly (E linear in p)
    r = occupancy_uncertainty_analytic(252, 3, 2)
    assert math.isclose(r["E_orders"], 100.8, rel_tol=1e-9)
    assert math.isclose(r["E_P_competent"], 0.4, rel_tol=1e-9)
    lo, hi = r["L_ci95_equal_tailed"]
    assert 15 < lo < 20 and 200 < hi < 206          # ~[17.0, 203.1]
    # law of total variance: within + between = total
    assert math.isclose(r["Var_within_E_p_Var"] + r["Var_between_Var_p_E"],
                        r["Var_total"], rel_tol=1e-9)


def test_beta_uncertainty_montecarlo_reports_mcse():
    from scijigsaw.contrib.probabilistic import occupancy_uncertainty
    res = occupancy_uncertainty(C.VAMP2, lambda p: Modification("x", p, _block_snap25),
                                alpha=3, beta=2, target=_fusion_competent,
                                n_draws=3000, rng_seed=1)
    assert res["method"] == "Monte Carlo"
    assert "E_orders_mcse" in res and res["E_orders_mcse"] > 0
    assert abs(res["E_orders_mean"] - 100.8) < 15   # MC estimate near analytic


def test_counterfactual_mixture_matches_closed_form():
    from scijigsaw.contrib.probabilistic import order_distribution_counterfactual
    from scijigsaw.assembly import Assembly
    H1 = Assembly(requires={"B": ["A"], "C": ["A"], "D": ["B", "C"]}, seed="s")
    H0 = Assembly(requires={"C": ["A"], "D": ["B", "C"]}, seed="s")
    L1, L0 = H1.n_orders_permitted(), H0.n_orders_permitted()
    for q in (0.0, 0.5, 1.0):
        r = order_distribution_counterfactual({"H1": H1, "H0": H0},
                                              {"H1": q, "H0": 1 - q})
        assert math.isclose(r.expected_orders, q * L1 + (1 - q) * L0, rel_tol=1e-9)
    # variance is zero at the deterministic ends, positive in between
    r_mid = order_distribution_counterfactual({"H1": H1, "H0": H0},
                                              {"H1": 0.5, "H0": 0.5})
    assert r_mid.variance() > 0


def test_counterfactual_prior_must_sum_to_one():
    from scijigsaw.contrib.probabilistic import order_distribution_counterfactual
    from scijigsaw.assembly import Assembly
    H = Assembly(requires={"B": ["A"]}, seed="s")
    with pytest.raises(ValueError):
        order_distribution_counterfactual({"H": H}, {"H": 0.7})


def test_disjoint_factorization_identity():
    from scijigsaw.contrib.probabilistic import linext_factorizes_disjoint
    # C(5,2)*1*1 = 10
    assert linext_factorizes_disjoint(2, 1, 3, 1) == 10
    # C(4,2)*2*1 = 12
    assert linext_factorizes_disjoint(2, 2, 2, 1) == 12
