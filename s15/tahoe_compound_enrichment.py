#!/usr/bin/env python3
"""
Drug-class enrichment at the CORRECT unit of independence.

The condition-level Fisher test (OR 43.2, p = 1.9e-11) treats each
(cell line, drug, dose, plate) as independent. It is not: the same compounds
recur across cell lines and doses, and DESeq2 controls are pooled within
plates. That p-value is descriptive, not inferential.

Three tests are run instead, in increasing conservatism:

  T1  COMPOUND-LEVEL, FULL CONTRACTION. One call per compound: does it index
      a model with ALL THREE positions closed? The weaker endpoint "any
      closing call" is power-limited -- ~50% of compounds produce one, so with
      a four-compound class the smallest attainable p is ~0.07 regardless of
      the data. Both are reported; the full-contraction endpoint is primary.

  T2  COMPOUND-LABEL PERMUTATION. Compound labels shuffled while the
      condition structure is held fixed; the observed proteasome share of
      closings is compared to the permutation null. Makes no independence
      assumption across conditions at all.

  T3  PLATE-CLUSTERED LOGISTIC. closed ~ proteasome, with standard errors
      clustered by plate.

Also reports the aggregation sensitivity the review asks for: does the
direction survive going from condition level to compound level?

    python3 tahoe_compound_enrichment.py --index delta_index.csv --tau 2.0
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats

# Atlas-annotated proteasome inhibitors (moa-fine == "Proteasome inhibitor")
ANNOTATED = {"Bortezomib", "Ixazomib", "Ixazomib citrate"}
# Added on external mechanistic literature, NOT atlas-annotated:
# auranofin inhibits the proteasome-associated DUBs UCHL5 and USP14.
EXTERNAL = {"Auranofin"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", default="delta_index.csv")
    ap.add_argument("--tau", type=float, default=2.0)
    ap.add_argument("--nperm", type=int, default=20000)
    ap.add_argument("--exclude-auranofin", action="store_true",
                    help="restrict the class to atlas-annotated inhibitors")
    a = ap.parse_args()

    cls = ANNOTATED if a.exclude_auranofin else (ANNOTATED | EXTERNAL)
    print(f"proteasome class ({len(cls)}): {sorted(cls)}")
    if not a.exclude_auranofin:
        print(f"  atlas-annotated: {sorted(ANNOTATED)}")
        print(f"  external literature: {sorted(EXTERNAL)}\n")

    r = pd.read_csv(a.index)
    r = r[r.tau == a.tau].copy()
    if r.empty:
        raise SystemExit(f"no rows at tau={a.tau}; available: "
                         f"{sorted(pd.read_csv(a.index).tau.unique())}")
    r["prot"] = r.drug.isin(cls)

    # ---------------- descriptive, condition level (NOT inferential) --------
    c = r[r.direction == "closed"]
    o = r[r.direction == "opened"]
    print(f"condition level (descriptive only)")
    print(f"  closed {len(c):4d}  proteasome {c.prot.sum():3d} "
          f"({c.prot.mean():.1%})")
    print(f"  opened {len(o):4d}  proteasome {o.prot.sum():3d} "
          f"({o.prot.mean():.1%})")

    # ---------------- T1: compound-level, FULL CONTRACTION ------------------
    print("\n" + "=" * 64)
    print("T1  COMPOUND-LEVEL, FULL CONTRACTION  (primary test)")
    print("=" * 64)
    print("Endpoint: does the compound index a model with ALL THREE positions")
    print("closed, i.e. full contraction to the constitutive-only model?")
    print("NOT 'any closing call' -- roughly half of all compounds produce one,")
    print("so that endpoint has a power ceiling near p = 0.07 with a class of")
    print("four compounds and cannot reach significance however clean the data.\n")

    g = (r.groupby("drug")
           .agg(full=("n", lambda s: (s[r.loc[s.index, "direction"] == "closed"]
                                      >= 3).any() if len(s) else False),
                anyclose=("direction", lambda s: (s == "closed").any()))
           .reset_index())
    g["prot"] = g.drug.isin(cls)
    n_cmpd = len(g)

    for lbl, col in (("full contraction (all 3 closed)", "full"),
                     ("any closing call", "anyclose")):
        tb = [[int((g.prot & g[col]).sum()), int((g.prot & ~g[col]).sum())],
              [int((~g.prot & g[col]).sum()), int((~g.prot & ~g[col]).sum())]]
        fe = stats.fisher_exact(tb, alternative="greater")
        base = tb[1][0] / max(tb[1][0] + tb[1][1], 1)
        star = "  <-- PRIMARY" if col == "full" else "  (power-limited)"
        print(f"  {lbl:34s} proteasome {tb[0][0]}/{tb[0][0]+tb[0][1]}   "
              f"other {tb[1][0]}/{tb[1][0]+tb[1][1]} ({base:.1%})")
        print(f"  {'':34s} OR = {fe[0]:.1f}   p = {fe[1]:.3g}{star}\n")

    print(f"  compounds with >=1 robust call: {n_cmpd}")

    # ---------------- T2: compound-label permutation ------------------------
    print("\n" + "-" * 64)
    print(f"T2  COMPOUND-LABEL PERMUTATION  ({a.nperm:,} permutations)")
    print("-" * 64)
    drugs = r.drug.unique()
    k = len(set(drugs) & cls)
    obs = c.prot.mean()
    rng = np.random.default_rng(0)
    dmap = {d: i for i, d in enumerate(drugs)}
    cidx = c.drug.map(dmap).values
    null = np.empty(a.nperm)
    for b in range(a.nperm):
        lab = np.zeros(len(drugs), bool)
        lab[rng.choice(len(drugs), k, replace=False)] = True
        null[b] = lab[cidx].mean()
    p = (1 + (null >= obs).sum()) / (a.nperm + 1)
    print(f"observed proteasome share of closings : {obs:.3%}")
    print(f"permutation null  mean {null.mean():.3%}  "
          f"95th {np.quantile(null,.95):.3%}  max {null.max():.3%}")
    print(f"permutation p = {p:.3g}"
          + ("   (at resolution limit)" if p <= 2/(a.nperm+1) else ""))

    # ---------------- T3: plate-clustered logistic --------------------------
    print("\n" + "-" * 64)
    print("T3  PLATE-CLUSTERED LOGISTIC  closed ~ proteasome")
    print("-" * 64)
    try:
        import statsmodels.api as sm
        d = r.copy()
        d["y"] = (d.direction == "closed").astype(int)
        X = sm.add_constant(d[["prot"]].astype(int))
        m = sm.Logit(d.y, X).fit(disp=0,
                                 cov_type="cluster",
                                 cov_kwds={"groups": d["plate"]})
        print(m.summary2().tables[1].to_string())
        print(f"\n  clusters (plates): {d['plate'].nunique()}")
        if d["plate"].nunique() < 20:
            print("  NOTE: few clusters -- cluster-robust SEs are anti-conservative "
                  "here.\n        Treat T1 and T2 as primary.")
    except ImportError:
        print("  pip install statsmodels")

    # ---------------- aggregation sensitivity -------------------------------
    print("\n" + "-" * 64)
    print("AGGREGATION SENSITIVITY (review item 3)")
    print("-" * 64)
    for lvl, keys in (("condition", None),
                      ("compound x cell line", ["drug", "Cell_Name_Vevo"]),
                      ("compound", ["drug"])):
        if keys is None:
            cc, oo = len(c), len(o)
            pc, po = c.prot.sum(), o.prot.sum()
        else:
            u = r.drop_duplicates(keys + ["direction"])
            cc = (u.direction == "closed").sum()
            oo = (u.direction == "opened").sum()
            pc = ((u.direction == "closed") & u.prot).sum()
            po = ((u.direction == "opened") & u.prot).sum()
        print(f"  {lvl:22s} closed {pc:3d}/{cc:4d} ({pc/max(cc,1):5.1%})   "
              f"opened {po:3d}/{oo:4d} ({po/max(oo,1):5.1%})")
    print("\n  direction must survive every level for the claim to hold.")

    os._exit(0)


if __name__ == "__main__":
    main()
