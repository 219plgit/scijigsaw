#!/usr/bin/env python3
"""
Full Tahoe-100M sweep + confound tests for the PSMB1/beta6 lag finding.

FINDING UNDER TEST (from a 150-shard pilot)
-------------------------------------------
In 7/7 conditions showing genuine proteasome bounce-back induction (2 cell
lines, 2 drugs, 3 doses), PSMB1 (beta6) was the LEAST-induced of the three
subunits. Support for beta6 joining first fell 0.333 -> 0.165-0.222.

The bounce-back is therefore coordinated in DIRECTION but not in MAGNITUDE.
Under Plackett-Luce a perfectly uniform response is exactly invisible
(scale invariance), so only the non-uniform residual is mechanistically
consequential. That is the claim to be tested here at full scale.

THREE PRE-SPECIFIED CONFOUND TESTS
----------------------------------
C1  BASELINE CEILING (the dangerous one). PSMB1 had the highest baseMean in
    the pilot. Abundant transcripts can show compressed log2FC through
    regression to the mean and saturation. If rank(baseMean) predicts
    rank(log2FC), "least induced" is an artifact, not biology.
      -> Spearman corr(baseMean, log2FC) within condition, and a test of
         whether PSMB1 lags AFTER conditioning on baseline rank.

C2  CELL-LINE DOMINANCE. 5 of 7 pilot hits were COLO 205. Two lines is not
    replication.
      -> per-cell-line sign test; report how many INDEPENDENT lines show it.

C3  DRUG SPECIFICITY. If PSMB1 lags under non-proteasome drugs too, it is a
    general property of the gene, not of the bounce-back.
      -> identical statistic on all non-proteasome drugs as background.

A finding that survives C1-C3 is reportable. One that fails C1 is not.

    pip install duckdb pandas scipy
    python tahoe_full_screen.py --shards 1026 --out-prefix full
"""
import argparse
import itertools
import math
import os
import time

REPO = "tahoebio/Tahoe-100M"
NSHARD = 1026
PROTEASOME_DRUGS = ("Bortezomib", "Ixazomib", "Ixazomib citrate")
SUB = {"PSMB6": "b1", "PSMB5": "b5", "PSMB1": "b6"}
MIN_CELLS = 50
FDR_MAX = 0.05


# ----------------------------------------------------------------- Jigsaw ---
def jigsaw(w):
    """(N_eff, {subunit: support for joining first}) over the six 20S orders.

    The admissible set is SIX in every condition. Weighting changes p(omega),
    never Omega.
    """
    ps = {}
    for o in itertools.permutations(["b1", "b5", "b6"]):
        p, r = 1.0, dict(w)
        for x in o:
            p *= r[x] / sum(r.values())
            del r[x]
        ps[o] = p
    H = -sum(p * math.log(p) for p in ps.values() if p > 0)
    first = {s: sum(p for o, p in ps.items() if o[0] == s)
             for s in ("b1", "b5", "b6")}
    return math.exp(H), first


# ------------------------------------------------------------------ fetch ---
def fetch(n, out_raw):
    import duckdb
    import pandas as pd
    from huggingface_hub import hf_hub_download
    os.environ["HF_HUB_DISABLE_PROGRESS_BARS"] = "1"

    con, frames, t0 = duckdb.connect(), [], time.time()
    genes = "','".join(SUB)
    for i in range(n):
        fn = (f"metadata/pseudobulk_differential_expression/"
              f"train-{i:05d}-of-{NSHARD:05d}.parquet")
        for attempt in range(4):
            try:
                p = hf_hub_download(REPO, filename=fn, repo_type="dataset")
                break
            except Exception:                            # noqa: BLE001
                if attempt == 3:
                    raise
                time.sleep(2 ** attempt * 5)
        # ALL drugs: the non-proteasome rows are the C3 background
        frames.append(con.execute(f"""
            SELECT gene_name, baseMean, log2FoldChange, padj,
                   n_cells_trt, n_cells_ctrl, drug, concentration,
                   Cell_Name_Vevo, plate
            FROM read_parquet('{p}')
            WHERE gene_name IN ('{genes}') AND log2FoldChange IS NOT NULL
        """).df())
        try:
            os.remove(os.path.realpath(p))
        except OSError:
            pass
        if (i + 1) % 50 == 0:
            print(f"  {i+1}/{n}  {sum(len(f) for f in frames):,} rows  "
                  f"{time.time()-t0:.0f}s", flush=True)

    df = pd.concat(frames, ignore_index=True)
    df.to_csv(out_raw, index=False)
    print(f"\n{len(df):,} rows -> {out_raw}  ({time.time()-t0:.0f}s)", flush=True)
    return df


# ----------------------------------------------------------------- screen ---
def score(df):
    import pandas as pd
    df = df.copy()
    df["su"] = df.gene_name.map(SUB)
    rows = []
    for k, s in df.groupby(["Cell_Name_Vevo", "drug", "concentration", "plate"]):
        if set(s.su) != {"b1", "b5", "b6"}:
            continue
        if s.n_cells_trt.min() < MIN_CELLS or s.n_cells_ctrl.min() < MIN_CELLS:
            continue
        lfc = dict(zip(s.su, s.log2FoldChange))
        bm = dict(zip(s.su, s.baseMean))
        q = dict(zip(s.su, s.padj))
        hi = max(lfc, key=lfc.get)
        lo = min(lfc, key=lfc.get)
        neff, first = jigsaw({x: 2 ** lfc[x] for x in lfc})
        rows.append({
            "line": k[0], "drug": k[1], "conc": k[2], "plate": k[3],
            "is_proteasome": k[1] in PROTEASOME_DRUGS,
            **{f"lfc_{x}": lfc[x] for x in lfc},
            **{f"base_{x}": bm[x] for x in bm},
            **{f"padj_{x}": q[x] for x in q},
            "spread": lfc[hi] - lfc[lo],
            "lowest": lo, "highest": hi,
            "induced": min(lfc.values()) > 0,      # true bounce-back
            "fdr_ok": bool(pd.notna(q[hi]) and pd.notna(q[lo])
                           and q[hi] <= FDR_MAX and q[lo] <= FDR_MAX),
            "N_eff": neff,
            **{f"first_{x}": first[x] for x in first},
            # rank of each subunit's BASELINE, for confound C1
            "base_rank_lowest": sorted(bm, key=bm.get).index(lo) + 1,
        })
    return pd.DataFrame(rows)


# ------------------------------------------------------------- confounds ---
def confounds(r):
    import pandas as pd
    from scipy import stats

    print("\n" + "=" * 68)
    print("PRIMARY: is beta6/PSMB1 the least-induced subunit under bounce-back?")
    print("=" * 68)
    prot = r[(r.is_proteasome) & (r.induced)]
    if prot.empty:
        print("no induced proteasome conditions -- nothing to test")
        return
    c = prot.lowest.value_counts()
    n, k = len(prot), int(c.get("b6", 0))
    p = stats.binomtest(k, n, 1 / 3, alternative="greater").pvalue
    print(f"induced proteasome conditions: {n}")
    print(f"  lowest = b6/PSMB1 : {k} ({k/n:.1%})   binomial p = {p:.3g}")
    for s in ("b1", "b5"):
        print(f"  lowest = {s:<10}: {int(c.get(s,0))}")
    print(f"\nN_eff  median {prot.N_eff.median():.3f}  min {prot.N_eff.min():.3f}"
          f"   (6.000 = no shift)")
    print(f"support(b6 first)  median {prot.first_b6.median():.3f}  (0.333 = none)")

    print("\n" + "-" * 68)
    print("C1  BASELINE CEILING  -- the dangerous confound")
    print("-" * 68)
    rho = []
    for _, s in prot.iterrows():
        bm = [s.base_b1, s.base_b5, s.base_b6]
        lf = [s.lfc_b1, s.lfc_b5, s.lfc_b6]
        if len(set(bm)) == 3:
            rho.append(stats.spearmanr(bm, lf).correlation)
    if rho:
        print(f"within-condition Spearman(baseMean, log2FC): median "
              f"{pd.Series(rho).median():+.3f}  (strong negative => artifact)")
    bmed = prot[["base_b1", "base_b5", "base_b6"]].median()
    print(f"median baseMean   b1={bmed.base_b1:.0f}  b5={bmed.base_b5:.0f}  "
          f"b6={bmed.base_b6:.0f}")
    hi_base = (prot[["base_b1", "base_b5", "base_b6"]].idxmax(axis=1)
               .str.replace("base_", "", regex=False))
    print(f"b6 has HIGHEST baseline in {(hi_base=='b6').mean():.1%} of conditions")
    both = ((prot.lowest == "b6") & (hi_base == "b6")).mean()
    print(f"b6 both highest-baseline AND least-induced: {both:.1%}")
    sub = prot[hi_base != "b6"]
    if len(sub):
        k2 = int((sub.lowest == "b6").sum())
        p2 = stats.binomtest(k2, len(sub), 1/3, alternative="greater").pvalue
        print(f"\n  *** conditions where b6 is NOT the most abundant: {len(sub)}")
        print(f"      b6 still least-induced in {k2} ({k2/len(sub):.1%}), "
              f"p = {p2:.3g}")
        print("      ^ THIS is the confound-free test. If it holds, C1 is passed.")
    else:
        print("\n  b6 is the most abundant in every condition -- C1 CANNOT be "
              "separated from the finding. Report the finding as confounded, "
              "or test with an abundance-matched gene set.")

    print("\n" + "-" * 68)
    print("C2  CELL-LINE DOMINANCE")
    print("-" * 68)
    per = prot.groupby("line").agg(n=("lowest", "size"),
                                   b6=("lowest", lambda x: (x == "b6").sum()))
    per["frac"] = per.b6 / per.n
    print(f"independent cell lines with >=1 induced condition: {len(per)}")
    print(f"lines where b6 lowest in ALL conditions: "
          f"{(per.frac == 1).sum()}/{len(per)}")
    print(per.sort_values('n', ascending=False).head(15).to_string())

    print("\n" + "-" * 68)
    print("C3  DRUG SPECIFICITY  -- background from non-proteasome drugs")
    print("-" * 68)
    bg = r[(~r.is_proteasome) & (r.induced)]
    if len(bg):
        kb = int((bg.lowest == "b6").sum())
        pb = stats.binomtest(kb, len(bg), 1/3, alternative="greater").pvalue
        print(f"induced NON-proteasome conditions: {len(bg)}")
        print(f"  lowest = b6 : {kb} ({kb/len(bg):.1%})  p = {pb:.3g}")
        print(f"  proteasome  : {k/n:.1%}   background: {kb/len(bg):.1%}")
        tbl = [[k, n - k], [kb, len(bg) - kb]]
        print(f"  Fisher exact (proteasome vs background): "
              f"p = {stats.fisher_exact(tbl, alternative='greater')[1]:.3g}")
        print("  ^ if background is also ~100%, the lag is a property of the "
              "GENE, not of proteasome inhibition.")
    else:
        print("no induced background conditions")

    print("\n" + "=" * 68)
    print("VERDICT: reportable only if C1's confound-free test holds AND C2 "
          "shows multiple independent lines AND C3 separates from background.")
    print("=" * 68)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shards", type=int, default=NSHARD)
    ap.add_argument("--out-prefix", default="full")
    ap.add_argument("--reuse", help="skip download, read this raw CSV")
    a = ap.parse_args()

    import pandas as pd
    df = pd.read_csv(a.reuse) if a.reuse else fetch(a.shards,
                                                    f"{a.out_prefix}_raw.csv")
    r = score(df)
    r.sort_values("spread", ascending=False).to_csv(
        f"{a.out_prefix}_screen.csv", index=False)
    print(f"scored {len(r):,} conditions "
          f"({r.is_proteasome.sum()} proteasome) -> {a.out_prefix}_screen.csv")
    confounds(r)
    os._exit(0)


if __name__ == "__main__":
    main()
