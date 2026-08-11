#!/usr/bin/env python3
"""
Real S15 proteasome result: screen Tahoe-100M for drug-dose-cell-line
conditions with DIFFERENTIAL response among beta1/beta5/beta6, then push the
survivors through the Scientific Jigsaw encounter-order weighting.

SUBUNIT -> GENE  (the numbering inverts; verified against HGNC/UniProt)
    beta1  = PSMB6      caspase-like
    beta5  = PSMB5      chymotrypsin-like, bortezomib target
    beta6  = PSMB1      non-catalytic
Only beta5 matches its own number. State this mapping in the manuscript.

WHY "DIFFERENTIAL" IS THE WHOLE POINT
------------------------------------
The encounter-order prior is Plackett-Luce, which is scale-invariant:
w -> c*w leaves every probability, the entropy and every support EXACTLY
unchanged. The NRF1 bounce-back is a COORDINATED upregulation, so it is the
exact null. Only a change in the RATIO among beta1:beta5:beta6 moves anything.

This screen therefore looks for conditions where the three genes move APART.

THE NULL IS A PUBLISHABLE OUTCOME
---------------------------------
If no condition shows differential movement beyond the background, that is a
result: the framework's predicted null, confirmed on 100M real cells. The
script reports the full distribution, not just the top hits, precisely so that
outcome is reportable rather than looking like a failed search.

GUARDING AGAINST SELECTION ON NOISE
-----------------------------------
Ranking tens of thousands of conditions and taking the extreme is a
garden-of-forking-paths. Three defences, all pre-specified:
  1. restrict to proteasome-inhibitor MOA before ranking;
  2. compute the same statistic on a matched background of non-proteasome
     drugs, and report proteasome hits against that background;
  3. require FDR support on the genes driving the spread.

Usage
-----
    python proteasome_differential.py --inspect
    python proteasome_differential.py --moa-pattern proteasome --out screen.csv
"""
import argparse
import csv
import itertools
import math
import sys
from collections import defaultdict

REPO = "tahoebio/Tahoe-100M"

SUBUNIT = {"PSMB6": "b1", "PSMB5": "b5", "PSMB1": "b6"}
GENES = list(SUBUNIT)

FDR_MAX = 0.05        # required on the genes driving the spread
MIN_CELLS = 50        # per (cell_line, drug, dose) group
MIN_SPREAD = 0.5      # log2 units between max and min subunit response


# ----------------------------------------------------------------- Jigsaw ---
def plackett_luce(order, w):
    p, rem = 1.0, dict(w)
    for x in order:
        p *= rem[x] / sum(rem.values())
        del rem[x]
    return p


def jigsaw(w):
    """Return (N_eff, supports) for the six admissible 20S orders of S14.2.

    The ADMISSIBLE SET IS ALWAYS SIX. Weighting never changes it. Any report of
    these numbers must say so, or it reproduces exactly the conflation the
    framework exists to prevent.
    """
    ps = {o: plackett_luce(o, w) for o in itertools.permutations(["b1", "b5", "b6"])}
    H = -sum(p * math.log(p) for p in ps.values() if p > 0)
    sup = lambda f: sum(p for o, p in ps.items() if f(o))
    return math.exp(H), {
        "13S+b1":     sup(lambda o: o[0] == "b1"),
        "13S+b5":     sup(lambda o: o[0] == "b5"),
        "13S+b1+b5":  sup(lambda o: set(o[:2]) == {"b1", "b5"}),
        "13S+b5+b6":  sup(lambda o: set(o[:2]) == {"b5", "b6"}),
    }


# ------------------------------------------------------------------- data ---
def inspect():
    """Print the schema so the column names below can be corrected."""
    from datasets import load_dataset
    for cfg in ("pseudobulk_differential_expression", "drug_metadata",
                "gene_metadata", "obs_metadata"):
        try:
            ds = load_dataset(REPO, name=cfg, split="train", streaming=True)
            row = next(iter(ds))
            print(f"\n=== {cfg} ===")
            for k, v in row.items():
                s = str(v)
                print(f"  {k:<32} {s[:70]}")
        except Exception as exc:                      # noqa: BLE001
            print(f"\n=== {cfg} === FAILED: {exc}")


def load_de(cols):
    from datasets import load_dataset
    ds = load_dataset(REPO, name="pseudobulk_differential_expression",
                      split="train", streaming=True)
    g, lfc, fdr = cols["gene"], cols["lfc"], cols["fdr"]
    drug, line, dose = cols["drug"], cols["cell_line"], cols.get("dose")
    ncell = cols.get("n_cells")

    want = set(GENES)
    rows = defaultdict(dict)
    for i, r in enumerate(ds):
        if i and i % 2_000_000 == 0:
            print(f"  ...{i:,} DE rows scanned", file=sys.stderr)
        if r.get(g) not in want:
            continue
        key = (r.get(line), r.get(drug), r.get(dose) if dose else "NA")
        rows[key][r[g]] = {
            "lfc": float(r[lfc]),
            "fdr": float(r[fdr]) if r.get(fdr) is not None else float("nan"),
            "n_cells": int(r.get(ncell, 0)) if ncell else None,
        }
    return rows


def screen(rows, moa_by_drug, moa_pattern):
    out = []
    for (line, drug, dose), rec in rows.items():
        if len(rec) < 3:
            continue                        # need all three subunits measured
        nc = [v["n_cells"] for v in rec.values() if v["n_cells"] is not None]
        if nc and min(nc) < MIN_CELLS:
            continue
        lfcs = {SUBUNIT[gene]: v["lfc"] for gene, v in rec.items()}
        hi = max(lfcs, key=lfcs.get)
        lo = min(lfcs, key=lfcs.get)
        spread = lfcs[hi] - lfcs[lo]
        # FDR required on the two genes DRIVING the spread
        inv = {v: k for k, v in SUBUNIT.items()}
        fdr_ok = all(rec[inv[s]]["fdr"] <= FDR_MAX for s in (hi, lo))
        moa = (moa_by_drug.get(drug) or "").lower()
        out.append({
            "cell_line": line, "drug": drug, "dose": dose, "moa": moa,
            "is_proteasome": moa_pattern.lower() in moa,
            "lfc_b1": lfcs["b1"], "lfc_b5": lfcs["b5"], "lfc_b6": lfcs["b6"],
            "spread_log2": spread, "driver_hi": hi, "driver_lo": lo,
            "fdr_ok": fdr_ok, "min_cells": min(nc) if nc else None,
        })
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inspect", action="store_true")
    ap.add_argument("--moa-pattern", default="proteasome")
    ap.add_argument("--out", default="screen.csv")
    ap.add_argument("--cols", nargs="*", default=[],
                    help="override schema, e.g. gene=gene_symbol lfc=log2FoldChange")
    args = ap.parse_args()

    if args.inspect:
        inspect()
        return

    cols = {"gene": "gene_symbol", "lfc": "log2FoldChange", "fdr": "padj",
            "drug": "drug", "cell_line": "cell_line", "dose": "dose",
            "n_cells": "n_cells"}
    for kv in args.cols:
        k, v = kv.split("=", 1)
        cols[k] = v

    from datasets import load_dataset
    dm = load_dataset(REPO, name="drug_metadata", split="train")
    moa_col = next((c for c in ("moa-fine", "moa_fine", "moa")
                    if c in dm.column_names), None)
    moa_by_drug = ({r["drug"]: r.get(moa_col) for r in dm} if moa_col else {})
    if not moa_col:
        print("WARNING: no MOA column found; --moa-pattern filter disabled",
              file=sys.stderr)

    rows = screen(load_de(cols), moa_by_drug, args.moa_pattern)
    if not rows:
        sys.exit("No conditions with all three subunits measured. Check --cols "
                 "against --inspect output.")

    prot = [r for r in rows if r["is_proteasome"]]
    bg = [r for r in rows if not r["is_proteasome"]]

    def pct(v, q):
        v = sorted(v)
        return v[min(len(v) - 1, int(q * len(v)))] if v else float("nan")

    bg_spread = [r["spread_log2"] for r in bg]
    print(f"\nconditions scored          : {len(rows):,}")
    print(f"  proteasome-inhibitor MOA : {len(prot):,}")
    print(f"  background               : {len(bg):,}")
    if bg_spread:
        print(f"background spread  median {pct(bg_spread,0.50):.3f}  "
              f"95th {pct(bg_spread,0.95):.3f}  99th {pct(bg_spread,0.99):.3f}")

    hits = sorted((r for r in prot if r["fdr_ok"]
                   and r["spread_log2"] >= MIN_SPREAD),
                  key=lambda r: -r["spread_log2"])
    print(f"proteasome hits passing FDR<={FDR_MAX} and spread>={MIN_SPREAD}: "
          f"{len(hits)}")

    if not hits:
        print("\n*** NULL RESULT -- this is reportable. ***\n"
              "No proteasome-inhibitor condition moves the three subunits "
              "apart beyond threshold. That is exactly what scale-invariance "
              "predicts for a coordinated bounce-back response, now shown on "
              "real data rather than asserted. Report the background "
              "distribution above alongside it.")

    for r in hits[:10]:
        w = {"b1": 2 ** r["lfc_b1"], "b5": 2 ** r["lfc_b5"], "b6": 2 ** r["lfc_b6"]}
        neff, sup = jigsaw(w)
        r["N_eff"] = round(neff, 3)
        for k, v in sup.items():
            r[f"support_{k}"] = round(v, 3)
        print(f"  {r['cell_line']:<12} {r['drug']:<22} dose={r['dose']:<8} "
              f"spread={r['spread_log2']:.2f}  N_eff={neff:.2f}")

    fields = sorted({k for r in rows for k in r})
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {len(rows)} scored conditions -> {args.out}")
    print("\nREMINDER: the admissible set is SIX orders in every condition. "
          "Weighting changes p(omega), never Omega.")


if __name__ == "__main__":
    main()
