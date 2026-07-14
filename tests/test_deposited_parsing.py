"""The extractor must work on real deposited structures.

Every test here would have FAILED at v1.0.0, because a crystallographic B-factor was
read as a confidence score and every experimental structure was therefore refused.
"""
import os
import shutil

import pytest

from scijigsaw.extract import extract, interface_residues, load_chains, looks_predicted

HERE = os.path.dirname(__file__)
EXP = os.path.join(HERE, "..", "examples", "structures_parsing")


def _run(tmp_path, pdb, chains="A,B", cutoff=5.0):
    d = tmp_path / pdb
    d.mkdir()
    shutil.copy(os.path.join(EXP, f"{pdb}.pdb"), d / f"{pdb}__partner.pdb")
    df, _ = extract(str(d), contact_cutoff=cutoff, chains=chains)
    return df.iloc[0]


def test_xray_complex_is_not_refused(tmp_path):
    """2XHE: a 567+220 residue X-ray complex. At v1.0.0 this was REFUSED."""
    r = _run(tmp_path, "2XHE")
    assert r.source == "experimental"
    assert not r.refused
    assert 40 <= r.n_interface_res <= 80, "implausible interface size"


def test_cryoem_complex_with_nonstandard_chain(tmp_path):
    """7DDO: the partner is chain C, not B."""
    r = _run(tmp_path, "7DDO", chains="A,C")
    assert not r.refused
    assert r.n_interface_res > 10


def test_amyloid_fibril_buries_its_whole_length(tmp_path):
    """2BEG: stacked beta-strands 4.8 A apart. EVERY residue must contact the
    neighbouring layer -- a structural prediction the geometry has to satisfy."""
    r = _run(tmp_path, "2BEG")
    assert r.n_interface_res == 26, "a fibril layer must be fully buried"


def test_experimental_structures_are_recognised_as_such():
    for pdb in ("2XHE", "7DDO", "2BEG", "1LCD"):
        assert looks_predicted(os.path.join(EXP, f"{pdb}.pdb")) is False


def test_interface_grows_monotonically_with_cutoff(tmp_path):
    """A 5 A interface must be a SUBSET of an 8 A one. If not, the contact search
    is wrong."""
    hub, par = load_chains(os.path.join(EXP, "2XHE.pdb"), "A", "B")
    s4, _ = interface_residues(hub, par, 4.0)
    s5, _ = interface_residues(hub, par, 5.0)
    s8, _ = interface_residues(hub, par, 8.0)
    assert s4 <= s5 <= s8
    assert len(s8) > 2 * len(s4), "the threshold must materially change the result"


def test_wrong_chain_fails_loudly(tmp_path):
    """Silence here would be worse than an error."""
    with pytest.raises(SystemExit):
        _run(tmp_path, "7DDO", chains="A,B")
