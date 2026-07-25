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
