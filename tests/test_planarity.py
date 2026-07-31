"""One flat tile layer suffices only while the graphs stay planar.

The manuscript (Supplementary Section S8) reports planarity of the CONNECTOR
graph the renderer draws -- not of the precedence poset used for counting, which
is a different graph. These tests call the same graph constructors the renderer
and the reproducibility script use, so there is a single construction path.
Node and edge counts are pinned as well: a planarity-only assertion would still
pass if an edge were accidentally dropped.
"""
import networkx as nx
import pandas as pd
import pytest

import scijigsaw.cases as C
from scijigsaw.render import Board

EXPECTED = {"VAMP2_tiles": (10, 12), "NLRP3": (10, 10),
            "30S": (20, 16), "30S_seeded": (21, 22)}


@pytest.fixture(scope="module")
def board():
    return Board(pd.read_csv("examples/vamp2/proteins.csv"),
                 pd.read_csv("examples/vamp2/interactions.csv"))


def test_vamp2_connector_graph_is_planar(board):
    G = board.connector_graph()
    assert (G.number_of_nodes(), G.number_of_edges()) == EXPECTED["VAMP2_tiles"]
    assert nx.check_planarity(G)[0]


def test_every_feasible_state_is_planar(board):
    states = board.feasible_states()
    assert states, "no feasible states derived"
    for state in states:
        g = board.connector_graph(state)
        assert nx.check_planarity(g)[0], f"state not planar: {sorted(state)}"


@pytest.mark.parametrize("label,case,seeded", [
    ("NLRP3", "INFLAMMASOME", False),
    ("30S", "RIBOSOME_30S_SPECIFIC", False),
    ("30S_seeded", "RIBOSOME_30S_SPECIFIC", True),
])
def test_dependency_graphs_are_planar(label, case, seeded):
    g = getattr(C, case).dependency_graph(include_seed=seeded)
    assert (g.number_of_nodes(), g.number_of_edges()) == EXPECTED[label]
    assert nx.check_planarity(g)[0]
