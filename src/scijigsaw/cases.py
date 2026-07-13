"""The three encoded assemblies of the paper.

Each is a CONSTRAINT ENCODING derived from the literature, not from structure.
The encoding is the modelling step, and it is where a reader should look first
if they disagree with a result.
"""
from .assembly import Assembly

# ---------------------------------------------------------------------------
# NLRP3 inflammasome -- a strictly nucleated cascade.
# Lu et al., Cell 156:1193 (2014); Xiao et al., Annu Rev Immunol 41:301 (2023).
# Poset depth 9 in 10 units: geometry is very nearly deterministic.
# ---------------------------------------------------------------------------
INFLAMMASOME = Assembly({
    "NLRP3":         set(),
    "NEK7":          {"NLRP3"},
    "NLRP3_disk":    {"NLRP3", "NEK7"},
    "NLRP3_PYDfil":  {"NLRP3_disk"},
    "ASC_PYDfil":    {"NLRP3_PYDfil"},
    "ASC_CARDfil":   {"ASC_PYDfil"},
    "Casp1_CARDfil": {"ASC_CARDfil"},
    "Casp1_cat":     {"Casp1_CARDfil"},
    "GSDMD":         {"Casp1_cat"},
    "IL1B":          {"Casp1_cat"},
})

# ---------------------------------------------------------------------------
# VAMP2 board -- bridge-rich but shallow.
# Complexin and Syt-1 are bridges (two sockets, spanning a seam).
# AP180/CALM overlap SNAP25 on the SNARE motif -> excluded while it is occupied.
# ---------------------------------------------------------------------------
VAMP2 = Assembly(
    requires={
        "SNAP25":        {"VAMP2"},                  # SNARE motif
        "Syntaxin-1A":   {"SNAP25"},
        "Munc18-1":      {"Syntaxin-1A"},
        "Complexin":     {"VAMP2", "SNAP25"},        # bridge
        "Syt-1":         {"SNAP25", "Syntaxin-1A"},  # bridge
        "Synaptophysin": {"VAMP2"},                  # vesicle face
        "SNCA":          {"VAMP2"},                  # N-terminus, aa 1-28
    },
    excludes=[({"SNAP25"}, {"AP180", "CALM"})],
    seed="VAMP2",
)

# ---------------------------------------------------------------------------
# 30S ribosomal subunit -- the Nomura assembly map.
# Mizushima & Nomura, Nature 226:1214 (1970); Held et al., JBC 249:3103 (1974).
#
# CONSERVATIVE ENCODING: a secondary requires SOME primary, a tertiary SOME
# secondary. The published map names specific dependencies, which would prune
# further, so our reduction is a LOWER BOUND. The poset is only three tiers deep
# and a secondary may follow ANY primary, so it is nearer a star than a chain:
# twenty units, pruned barely sixfold. This is the counterexample that bounds
# the method.
# ---------------------------------------------------------------------------
PRIMARY_30S = ["S4", "S7", "S8", "S15", "S17", "S20"]
SECONDARY_30S = ["S5", "S6", "S9", "S12", "S13", "S16", "S18", "S19"]
TERTIARY_30S = ["S2", "S3", "S10", "S11", "S14", "S21"]


def count_30S():
    """Exact count for the tiered 'any-of' hierarchy.

    The Assembly class encodes precedence on NAMED partners; the Nomura map is
    an 'any-of' rule, so we count it directly by DP over tier occupancies.
    """
    from functools import lru_cache
    import math
    P, Q, T = len(PRIMARY_30S), len(SECONDARY_30S), len(TERTIARY_30S)

    @lru_cache(maxsize=None)
    def f(p, q, t):
        if (p, q, t) == (P, Q, T):
            return 1
        tot = 0
        if p < P:
            tot += (P - p) * f(p + 1, q, t)
        if q < Q and p >= 1:                 # a secondary needs SOME primary
            tot += (Q - q) * f(p, q + 1, t)
        if t < T and q >= 1:                 # a tertiary needs SOME secondary
            tot += (T - t) * f(p, q, t + 1)
        return tot

    n = P + Q + T
    permitted = f(0, 0, 0)
    total = math.factorial(n)
    return dict(n=n, depth=3, total=total, permitted=permitted,
                reduction=total / permitted)


ALL_CASES = {"inflammasome": INFLAMMASOME, "vamp2": VAMP2}
