"""scijigsaw -- interface geometry as a constraint on protein-assembly order.

    from scijigsaw import Assembly, VAMP2, INFLAMMASOME
    VAMP2.summary()
    # {'n': 7, 'depth': 3, ..., 'permitted': 336, 'reduction': 15.0}

What the tool does
------------------
It converts mapped protein interfaces into complementary tabs and sockets, and
counts exactly the assembly orders permitted by two geometric constraints:
precedence (a piece is seated only when the partners it binds are present) and
exclusion (partners on overlapping interfaces cannot coexist on the board).

What it does NOT do
-------------------
It does not infer the assembly pathway. It eliminates orders incompatible with
the encoded interface and seating constraints. It cannot encode catalysis: NSF
and alpha-SNAP, which vacate the VAMP2 SNARE motif, have no representation as a
static tab-and-socket relation. Geometry bounds the set of admissible
configurations; enzymes, concentration and kinetics select among them.
"""
from .assembly import Assembly
from .cases import (INFLAMMASOME, VAMP2, ALL_CASES, count_30S,
                    RIBOSOME_30S_SPECIFIC)
from .geometry import Piece, coverage, TAB, SOCKET, FLAT
from .results import all_results
from .render import Board, render

__version__ = "1.5.0"
__all__ = ["Assembly", "Piece", "Board", "render", "coverage", "TAB", "SOCKET", "FLAT",
           "INFLAMMASOME", "VAMP2", "ALL_CASES", "count_30S",
           "RIBOSOME_30S_SPECIFIC", "all_results", "__version__"]
