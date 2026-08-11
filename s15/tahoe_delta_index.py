#!/usr/bin/env python3
"""
Condition indexing from Tahoe: V_z -> M_z -> Omega_z.

STATISTIC
---------
For each beta-ring position j with constitutive gene c and inducible gene i:

    Delta_j = LFC(i) - LFC(c)
            = log2 [ (i/c)_treated / (i/c)_control ]

Each LFC is a WITHIN-GENE treated-vs-control contrast, so gene-specific capture
efficiency, transcript length and expression scale largely cancel before the
two are compared. This is why Delta is used and NOT a cross-gene baseMean
ratio: DESeq2 normalised counts are not calibrated across genes, so an absolute
threshold on baseMean(i)/baseMean(c) cannot be read as a protein-occupancy
ratio.

INFERENCE, not thresholding alone
---------------------------------
The DE table carries lfcSE, so Delta has a standard error and a test:

    SE(Delta) = sqrt(SE_i^2 + SE_c^2)      z = Delta / SE(Delta)

(approximate: the two fits share a dispersion estimate, so this is mildly
conservative or anticonservative depending on correlation; it is reported as
approximate and used alongside, not instead of, the threshold sweep.)
Benjamini-Hochberg across all positions and conditions.

WHAT DELTA DOES AND DOES NOT ESTABLISH
--------------------------------------
Delta measures the CHANGE in the inducible/constitutive balance, not the
balance itself. A shift from 0.01 to 0.02 is a large Delta with no meaningful
availability change; a position sitting at parity with Delta = 0 has both
variants available and no change. Delta therefore cannot by itself define
membership of the candidate set.

Accordingly this script does not declare occupancy. It emits, for each
condition and each threshold tau, a CANDIDATE component set V_z under the
declared criterion, and leaves the resulting family of models M_1..M_k to the
existing model-uncertainty machinery. The reportable result is ROBUSTNESS:
which changes in Omega hold across every reasonable tau.

    python3 tahoe_delta_index.py --raw immuno_raw.csv --out delta_index.csv
"""
import argparse
import itertools
import os

import numpy as np
import pandas as pd
from scipy import stats

POSITIONS = {                     # position: (constitutive, inducible)
    "b1": ("PSMB6", "PSMB9"),
    "b2": ("PSMB7", "PSMB10"),
    "b5": ("PSMB5", "PSMB8"),
}
TAUS = [0.5, 1.0, 1.5, 2.0]       # pre-declared switch thresholds, log2 units
FDR_MAX = 0.05
MIN_CELLS = 50
KEYS = ["Cell_Name_Vevo", "drug", "concentration", "plate"]


def load(path):
    d = pd.read_csv(path)
    need = {"gene_name", "log2FoldChange", "Cell_Name_Vevo", "drug",
            "concentration", "plate"}
    missing = need - set(d.columns)
    if missing:
        raise SystemExit(f"missing columns: {sorted(missing)}")
    if "lfcSE" not in d.columns:
        print("NOTE: lfcSE absent from the raw file -- z-tests on Delta cannot "
              "be computed.\n      Re-fetch including lfcSE to enable "
              "inference; thresholds still work.\n")
    if {"n_cells_trt", "n_cells_ctrl"} <= set(d.columns):
        d = d[(d.n_cells_trt >= MIN_CELLS) & (d.n_cells_ctrl >= MIN_CELLS)]
    return d


def deltas(d):
    have_se = "lfcSE" in d.columns
    vals = ["log2FoldChange"] + (["lfcSE"] if have_se else [])
    w = d.pivot_table(index=KEYS, columns="gene_name", values=vals)
    rows = []
    for pos, (c, i) in POSITIONS.items():
        try:
            lc, li = w[("log2FoldChange", c)], w[("log2FoldChange", i)]
        except KeyError:
            print(f"  position {pos}: {c} or {i} absent -- skipped")
            continue
        r = pd.DataFrame({"position": pos, "delta": li - lc}, index=w.index)
        if have_se:
            se = np.sqrt(w[("lfcSE", c)] ** 2 + w[("lfcSE", i)] ** 2)
            r["se"] = se
            r["z"] = r.delta / se
            r["p"] = 2 * stats.norm.sf(np.abs(r.z))
        rows.append(r.reset_index())
    out = pd.concat(rows, ignore_index=True).dropna(subset=["delta"])
    if "p" in out.columns:
        m = out.p.notna()
        q = np.full(len(out), np.nan)
        pv = out.loc[m, "p"].values
        order = np.argsort(pv)
        adj = pv[order] * len(pv) / (np.arange(len(pv)) + 1)
        adj = np.minimum.accumulate(adj[::-1])[::-1]
        tmp = np.empty(len(pv))
        tmp[order] = np.clip(adj, 0, 1)
        q[np.where(m)[0]] = tmp
        out["fdr"] = q
    return out


def candidate_sets(dl, tau, use_fdr=True):
    """V_z under the declared criterion, in BOTH directions.

    A ring position is OPEN when both variants are candidates, CLOSED when only
    one is. A perturbation can move it either way:

      OPEN  (delta >= +tau): the balance shifts toward the inducible variant,
            which enters the candidate set. Omega grows.
      CLOSE (delta <= -tau): the balance shifts toward the constitutive
            variant, and the inducible variant leaves. Omega shrinks.

    The closing direction is the one the data actually shows (bortezomib
    shifts all three positions strongly constitutive), and the original
    gain-only criterion was blind to it.

    IMPORTANT ASYMMETRY IN WHAT THESE LICENSE. Declaring a position OPENED
    needs only the shift. Declaring one CLOSED presupposes it was open at
    baseline -- which requires a baseline LEVEL, and DESeq2 baseMean is not
    calibrated across genes. Closing calls are therefore reported relative to
    a declared baseline-open assumption, and both assumptions are carried
    through as alternative models rather than resolved here.
    """
    if use_fdr and "fdr" in dl.columns:
        dl = dl[dl.fdr <= FDR_MAX]
    opened = (dl[dl.delta >= tau].groupby(KEYS)["position"]
              .apply(frozenset).to_dict())
    closed = (dl[dl.delta <= -tau].groupby(KEYS)["position"]
              .apply(frozenset).to_dict())
    return opened, closed


def omega(n_open):
    """Compositions = 2^(positions with both variants candidate); orders = 3!.

    PLACEHOLDER for orders-per-composition -- replace with the human contact
    graph once extracted.
    """
    return 2 ** n_open, (2 ** n_open) * 6


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default="immuno_raw.csv")
    ap.add_argument("--out", default="delta_index.csv")
    a = ap.parse_args()

    d = load(a.raw)
    dl = deltas(d)
    print(f"{len(dl):,} position-conditions from "
          f"{dl.groupby(KEYS).ngroups:,} conditions\n")

    print("Delta distribution by position (log2 units):")
    print(dl.groupby("position").delta
          .describe()[["25%", "50%", "75%", "max"]].round(3).to_string())
    if "fdr" in dl.columns:
        print("\nFDR-supported |Delta| by position:")
        for pos, g in dl.groupby("position"):
            s = g[g.fdr <= FDR_MAX]
            print(f"  {pos}: {len(s):5d}/{len(g):5d} significant "
                  f"({len(s)/len(g):.1%}), median delta {s.delta.median():+.3f}")

    print("\n" + "=" * 62)
    print("THRESHOLD SWEEP -- conditions gaining >=1 inducible candidate")
    print("=" * 62)
    print(f"{'tau':>5}{'opened':>9}{'closed':>9}{'max |Omega|':>13}"
          f"{'min |Omega| (closing)':>24}")
    recs = []
    total = dl.groupby(KEYS).ngroups
    for tau in TAUS:
        op, cl = candidate_sets(dl, tau)
        mx = max([omega(len(v))[1] for v in op.values()] or [6])
        # closing starts from all three open (Omega=48) and removes positions
        mn = min([omega(3 - len(v))[1] for v in cl.values()] or [48])
        print(f"{tau:>5}{len(op):>9}{len(cl):>9}{mx:>13}{mn:>24}")
        for k, v in op.items():
            recs.append(dict(zip(KEYS, k), tau=tau, direction="opened",
                             positions="|".join(sorted(v)), n=len(v),
                             omega=omega(len(v))[1]))
        for k, v in cl.items():
            recs.append(dict(zip(KEYS, k), tau=tau, direction="closed",
                             positions="|".join(sorted(v)), n=len(v),
                             omega=omega(3 - len(v))[1]))
    print(f"\n(total conditions scored: {total:,}; "
          f"baseline Omega with no inducible candidate = 6)")

    rec = pd.DataFrame([r for r in recs if r])
    if rec.empty:
        print("\n*** NO condition crosses any threshold. ***\n"
              "Under this declared criterion the candidate set never changes, "
              "so Omega is invariant across the perturbation grid, in EITHER "
              "direction. That is a reportable null: small-molecule "
              "perturbation in this panel does not index a different "
              "immunoproteasome model.")
    else:
        rec.to_csv(a.out, index=False)
        stable = (rec.groupby(KEYS + ["direction", "positions"])
                  .tau.nunique() == len(TAUS))
        print(f"\nROBUST across all {len(TAUS)} thresholds: {int(stable.sum())} "
              f"(condition, gained-set) pairs")
        print("These are the only ones whose Omega change should be reported "
              "as a result.")
        for d_ in ("opened", "closed"):
            top = rec[(rec.tau == max(TAUS)) & (rec.direction == d_)]
            if len(top):
                print(f"\nStrongest {d_} (tau={max(TAUS)}):")
                print(top.sort_values("n", ascending=False).head(10)
                      [KEYS + ["positions", "omega"]].to_string(index=False))
        print(f"\nwrote {a.out}")

    dl.to_csv(a.out.replace(".csv", "_deltas.csv"), index=False)
    print(f"wrote {a.out.replace('.csv', '_deltas.csv')}")
    os._exit(0)


if __name__ == "__main__":
    main()
