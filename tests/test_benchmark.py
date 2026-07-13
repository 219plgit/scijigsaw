from scijigsaw.benchmark import run, analyse


def test_benchmark_reproduces_the_papers_finding():
    """Small but deterministic: depth must beat n as a predictor."""
    A = run(n_posets=150, seed=0)
    r = analyse(A)
    assert r["n_posets"] > 100
    # depth is the stronger marginal predictor of log10(reduction)
    assert r["pearson"]["depth"] > 0.85
    assert r["pearson"]["n"] < r["pearson"]["depth"]
    # and the collinearity we report is real, not hidden
    assert r["vif"]["depth"] > 2.0


def test_reduction_is_never_below_one():
    A = run(n_posets=60, seed=1)
    assert (A[:, 5] >= 0).all()     # log10(reduction) >= 0, i.e. reduction >= 1
