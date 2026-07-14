"""Failures that only real structures expose.

The synthetic complexes in examples/structures/ verify that the extractor
recovers interfaces it was given. They cannot expose the things every real file
has and no toy file does. Each test below FAILED before the fix.
"""
import os
import shutil

import pytest

from scijigsaw.extract import MODIFIED_AA, extract, looks_predicted

HERE = os.path.dirname(__file__)
REAL = os.path.join(HERE, "..", "examples", "structures_format_cases")


def _one(tmp_path, name):
    d = tmp_path / name
    d.mkdir()
    shutil.copy(os.path.join(REAL, f"{name}__partner.pdb"), d)
    df, _ = extract(str(d))
    return df.iloc[0]


def test_experimental_bfactor_is_not_a_confidence_score(tmp_path):
    """THE CRITICAL ONE.

    A crystallographic B-factor is a temperature factor (~28 A^2 for a good
    structure), not pLDDT. Reading it as confidence refuses every well-ordered
    experimental structure ever deposited. Before the fix this row was
    refused=True with 'confidence' 29.9.
    """
    r = _one(tmp_path, "EXPT")
    assert r.source == "experimental"
    assert not r.refused, "an experimental structure must never be refused for low B"
    assert r.n_interface_res == 21


def test_predicted_models_still_use_plddt(tmp_path):
    """The fix must not disable confidence filtering for actual predictions."""
    r = _one(tmp_path, "SEMET")
    assert r.source == "predicted"
    assert r.confidence == pytest.approx(90.0, abs=1.0)


def test_modified_residues_are_not_silently_dropped(tmp_path):
    """MSE (selenomethionine) is a HETATM but IS an amino acid. Six MSE residues
    sit in this interface; before the fix all six vanished (10 found, not 16)."""
    assert "MSE" in MODIFIED_AA
    r = _one(tmp_path, "SEMET")
    assert r.n_interface_res == 16, "MSE residues must count toward the interface"


def test_insertion_codes_are_distinct_residues(tmp_path):
    """Kabat-numbered antibodies carry 52, 52A, 52B. Keying on the residue number
    alone conflates them. Before the fix: 6 residues; after: 9."""
    r = _one(tmp_path, "ICODE")
    assert r.n_interface_res == 9


def test_predicted_vs_experimental_detection():
    assert looks_predicted(os.path.join(REAL, "SEMET__partner.pdb")) is True
    assert looks_predicted(os.path.join(REAL, "EXPT__partner.pdb")) is False
