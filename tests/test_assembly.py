"""The counts in the paper are pinned here. If they change, the paper is wrong."""
import math

import pytest

from scijigsaw import Assembly, INFLAMMASOME, VAMP2, count_30S


def test_vamp2_matches_the_paper():
    s = VAMP2.summary()
    assert s["n"] == 7
    assert s["depth"] == 3
    assert s["total"] == 5040
    assert s["permitted"] == 252           # 95.0 % eliminated (corrected by 1KIL)
    assert round(s["reduction"]) == 20


def test_inflammasome_matches_the_paper():
    s = INFLAMMASOME.summary()
    assert s["n"] == 10
    assert s["depth"] == 9
    assert s["total"] == 3628800
    assert s["permitted"] == 2             # geometry is nearly deterministic
    assert round(s["reduction"]) == 1814400


def test_30S_matches_the_paper():
    s = count_30S()
    assert s["n"] == 20
    assert s["depth"] == 3
    assert 5.5 < s["reduction"] < 6.5      # twenty units, pruned only ~6x


def test_depth_not_size_governs_pruning():
    """The paper's central claim, as an assertion."""
    infl, ribo = INFLAMMASOME.summary(), count_30S()
    assert infl["n"] < ribo["n"]                    # fewer units...
    assert infl["depth"] > ribo["depth"]            # ...but deeper
    assert infl["reduction"] > 1e5 * ribo["reduction"]   # and prunes vastly more


def test_chain_admits_exactly_one_order():
    chain = Assembly({f"u{i}": ({f"u{i-1}"} if i else set()) for i in range(6)})
    assert chain.n_orders_permitted() == 1
    assert chain.depth() == 6


def test_star_admits_every_order():
    star = Assembly({f"u{i}": set() for i in range(6)})
    assert star.n_orders_permitted() == math.factorial(6)
    assert star.reduction() == 1.0
    assert star.depth() == 1


def test_exclusion_removes_orders():
    free = Assembly({"a": set(), "b": set(), "c": set()})
    excl = Assembly({"a": set(), "b": set(), "c": set()},
                    excludes=[({"a"}, {"b"})])
    assert excl.n_orders_permitted() < free.n_orders_permitted()


def test_bridge_needs_its_seam_closed():
    """A bridge (two parents) cannot be seated before both are present."""
    A = Assembly({"a": set(), "b": set(), "bridge": {"a", "b"}})
    assert A.bridges() == ["bridge"]
    assert A.n_orders_permitted() == 2      # ab|bridge, ba|bridge -- never first


def test_cycles_are_rejected():
    with pytest.raises(ValueError, match="cyclic"):
        Assembly({"a": {"b"}, "b": {"a"}})


def test_refuses_intractable_sizes():
    """Exact counting is O(2^n n); the tool refuses rather than pretending."""
    big = Assembly({f"u{i}": set() for i in range(26)})
    with pytest.raises(ValueError, match="impractical"):
        big.n_orders_permitted()


def test_unsatisfiable_constraints_raise():
    A = Assembly({"a": set(), "b": set()}, excludes=[({"a"}, {"b"})])
    # a and b mutually exclusive and both required -> no complete order
    assert A.n_orders_permitted() == 0
    with pytest.raises(ValueError, match="unsatisfiable"):
        A.reduction()


def test_30S_result_is_not_an_artefact_of_a_weak_encoding():
    """The paper's headline contrast must survive the STRONGER encoding.

    The conservative 'any-of' tier encoding gives 6x. The published map names
    specific parents and gives 45,545x -- a 7,500-fold difference. If the contrast
    with the inflammasome depended on the conservative choice, there would be no
    result. It does not.
    """
    from scijigsaw import INFLAMMASOME, RIBOSOME_30S_SPECIFIC, count_30S

    loose = count_30S()["reduction"]
    strict = RIBOSOME_30S_SPECIFIC.summary()
    infl = INFLAMMASOME.summary()

    assert 5.5 < loose < 6.5                       # the conservative encoding
    assert 4e4 < strict["reduction"] < 5e4         # the published dependencies
    assert strict["depth"] == 3                    # still SHALLOW either way
    assert strict["n"] == 20

    # half the subunits, more than an order of magnitude more pruning
    assert infl["n"] < strict["n"]
    assert infl["reduction"] > 10 * strict["reduction"]


def test_the_structural_correction_is_recorded():
    """The 1KIL correction is a reported RESULT and must not silently regress.

    Our literature encoding had Complexin bridging VAMP2 and SNAP25. The deposited
    complexin-SNARE complex shows complexin contacting synaptobrevin (13 residues) and
    syntaxin (10), and NOT SNAP25. Permitted orders fell from 336 to 252.
    """
    from scijigsaw import VAMP2, all_results

    assert VAMP2.requires["Complexin"] == {"VAMP2", "Syntaxin-1A"}
    assert "VAMP2" in VAMP2.requires["Syntaxin-1A"], (
        "1KIL shows a direct VAMP2-syntaxin contact (34 residues); the SNARE core is "
        "a four-helix bundle, not a chain")
    assert VAMP2.summary()["permitted"] == 252

    r = all_results(60, 0)["structural_correction"]
    assert r["permitted_before"] == 336 and r["permitted_after"] == 252
    assert r["complexin_contacts_snap25_residues"] == 0
