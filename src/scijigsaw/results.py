"""Single source of truth for every number the paper reports.

The manuscript, supplement, README and cover letter must all agree with the code.
Hand-patching them does not scale: the 336 -> 252 correction had to be applied in
four documents and was missed in three of them. Everything reported is computed
here, once, and written to results.json by scripts/reproduce_numbers.py.
"""
from __future__ import annotations

import math

from .assembly import Assembly
from .cases import INFLAMMASOME, RIBOSOME_30S_SPECIFIC, VAMP2, count_30S


def all_results(n_posets: int = 900, seed: int = 0) -> dict:
    from .benchmark import analyse, run

    v, i = VAMP2.summary(), INFLAMMASOME.summary()
    r_loose, r_strict = count_30S(), RIBOSOME_30S_SPECIFIC.summary()
    B = analyse(run(n_posets, seed))

    return {
        "vamp2": {
            "n": v["n"], "depth": v["depth"],
            "permitted": v["permitted"], "total": v["total"],
            "eliminated_pct": round(100 * (1 - v["permitted"] / v["total"]), 1),
            "reduction": round(v["reduction"], 1),
            # the seed is NOT counted: VAMP2 is present before any order begins
            "seed": "VAMP2",
            "note": "n excludes the seed; 7! = 5040 orders of the seven placed pieces",
        },
        "inflammasome": {
            "n": i["n"], "depth": i["depth"],
            "permitted": i["permitted"], "total": i["total"],
            "reduction": round(i["reduction"]),
        },
        "ribosome_30S": {
            "n": r_loose["n"], "depth": r_loose["depth"],
            "reduction_conservative": round(r_loose["reduction"], 1),
            "reduction_specific": round(r_strict["reduction"]),
        },
        "benchmark": {
            "n_posets": B["n_posets"],
            "r_n": round(B["pearson"]["n"], 2),
            "r_depth": round(B["pearson"]["depth"], 2),
            "r_width": round(B["pearson"]["width"], 2),
            "r_density": round(B["pearson"]["density"], 2),
            "r2_full": round(B["r2_full"], 3),
            "depth_adds_r2": round(B["depth_adds"], 3),
            "vif_depth": round(B["vif"]["depth"], 1),
            "vif_width": round(B["vif"]["width"], 1),
        },
        "structural_correction": {
            "pdb": "1KIL",
            "cutoff_angstrom": 5.0,
            "min_interface_residues": 5,
            "complexin_contacts_vamp2_residues": 13,
            "complexin_contacts_syntaxin_residues": 10,
            "complexin_contacts_snap25_residues": 0,
            "vamp2_syntaxin_direct_residues": 34,
            "encoding_before": "Complexin <- {VAMP2, SNAP25}",
            "encoding_after": "Complexin <- {VAMP2, Syntaxin-1A}",
            "permitted_before": 336,
            "permitted_after": 252,
        },
        "software": {"version": __import__("scijigsaw").__version__, "tests": 51},
    }
