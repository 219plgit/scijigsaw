"""Random-poset benchmark: what governs the pruning?

Three encoded assemblies illustrate a contrast but cannot establish it. Here we
generate random constraint posets, varying n, depth, width and edge density
independently, count the linear extensions exactly, and ask whether depth still
predicts the reduction once n and density are controlled.

Finding: n is a weak predictor (r ~ 0.47), depth a strong one (r ~ 0.94), and at
fixed n depth separates the reduction by ~6 orders of magnitude. But depth, width
and density are COLLINEAR (a deep poset must be narrow: depth x width ~ n), so
they are alternative measures of one property -- the SEQUENTIALITY of the
assembly. Depth is not a privileged variable independent of the others.
"""
from __future__ import annotations

import math
import random
from typing import List

import numpy as np

from .assembly import Assembly


def random_poset(n: int, target_depth: int, edge_p: float, rng: random.Random):
    """Layered random DAG. Each unit gets a layer (bounding the depth), then
    draws one mandatory parent from a strictly earlier layer (guaranteeing the
    depth) plus further parents with probability edge_p."""
    layers = [0] + [rng.randint(0 if target_depth <= 1 else 1, max(target_depth - 1, 0))
                    for _ in range(n - 1)]
    layers.sort()
    requires = {}
    names = [f"u{i}" for i in range(n)]
    for i in range(n):
        cand = [j for j in range(i) if layers[j] < layers[i]]
        par = []
        if layers[i] > 0 and cand:
            par.append(rng.choice(cand))
            for j in cand:
                if j not in par and rng.random() < edge_p:
                    par.append(j)
        requires[names[i]] = {names[j] for j in par}
    return Assembly(requires)


def run(n_posets: int = 900, seed: int = 0,
        n_range=(9, 15), edge_p_range=(0.0, 0.6)) -> np.ndarray:
    """Returns columns: n, depth, width, density, bridge_frac, log10(reduction)."""
    rng = random.Random(seed)
    rows: List[list] = []
    for _ in range(n_posets):
        n = rng.randint(*n_range)
        td = rng.randint(1, n)
        ep = rng.uniform(*edge_p_range)
        A = random_poset(n, td, ep, rng)
        perm = A.n_orders_permitted()
        if perm == 0:
            continue
        edges = sum(len(r) for r in A.requires.values())
        density = edges / (n * (n - 1) / 2)
        rows.append([n, A.depth(), A.width(), density,
                     len(A.bridges()) / n,
                     math.log10(math.factorial(n) / perm)])
    return np.array(rows)


def analyse(A: np.ndarray) -> dict:
    n_, d_, w_, dens_, br_, y = A.T
    out = {"n_posets": len(A),
           "pearson": {k: float(np.corrcoef(x, y)[0, 1]) for k, x in
                       [("n", n_), ("depth", d_), ("width", w_),
                        ("density", dens_), ("bridge_frac", br_)]}}
    X = np.column_stack([np.ones(len(A)), n_, d_, w_, dens_])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    r2 = 1 - (resid ** 2).sum() / ((y - y.mean()) ** 2).sum()
    X0 = np.column_stack([np.ones(len(A)), n_, w_, dens_])
    b0, *_ = np.linalg.lstsq(X0, y, rcond=None)
    r2_0 = 1 - ((y - X0 @ b0) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    # variance inflation
    vif = {}
    Z = np.column_stack([n_, d_, w_, dens_])
    for i, nm in enumerate(["n", "depth", "width", "density"]):
        others = np.column_stack([np.ones(len(A))] +
                                 [Z[:, j] for j in range(4) if j != i])
        b, *_ = np.linalg.lstsq(others, Z[:, i], rcond=None)
        rr = 1 - ((Z[:, i] - others @ b) ** 2).sum() / \
            ((Z[:, i] - Z[:, i].mean()) ** 2).sum()
        vif[nm] = float(1 / max(1 - rr, 1e-9))
    out.update(r2_full=float(r2), r2_without_depth=float(r2_0),
               depth_adds=float(r2 - r2_0), vif=vif)
    return out
