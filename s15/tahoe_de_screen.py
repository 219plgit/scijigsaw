#!/usr/bin/env python3
"""
Extract beta1/beta5/beta6 rows from the Tahoe-100M pseudobulk DE table and
screen for DIFFERENTIAL response, then run survivors through the Scientific
Jigsaw encounter-order weighting.

Written against the REAL schema (verified 2026):
    gene_name, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj,
    plate, n_cells_trt, n_cells_ctrl, Cell_ID_Cellosaur, Cell_ID_DepMap,
    drug, concentration, concentration_unit, Cell_Name_Vevo
1,026 shards x ~3.99M rows = ~4.1 billion rows. We need 3 genes, so DuckDB
does projection + predicate pushdown over HTTP range requests.

SUBUNIT -> GENE  (the numbering inverts)
    beta1 = PSMB6   beta5 = PSMB5   beta6 = PSMB1
Only beta5 matches its own number.

WEIGHTS ARE RELATIVE TO EACH CELL LINE'S OWN CONTROL
----------------------------------------------------
w_x = 2^log2FoldChange_x, so control is (1,1,1) by construction. This measures
the drug-induced CHANGE in relative availability, not absolute stoichiometry.
That is deliberate: DESeq2 normalised counts are not comparable across genes
(no length/efficiency correction), so absolute cross-gene ratios would be
biased. Ratios of fold-changes are not.

Plackett-Luce is scale-invariant, so a COORDINATED response gives exactly the
control answer. The NRF1 bounce-back is coordinated. A null here is therefore
the predicted result, and is reportable.

    pip install duckdb pandas

    python tahoe_de_screen.py --shards 1 --out probe.csv      # measure cost
    python tahoe_de_screen.py --shards 0 --out screen.csv     # 0 = all 1026
"""
import argparse
import itertools
import math
import os
import sys
import time

BASE = ("https://huggingface.co/datasets/tahoebio/Tahoe-100M/resolve/main/"
        "metadata/pseudobulk_differential_expression")
NSHARD = 1026

SUBUNIT = {"PSMB6": "b1", "PSMB5": "b5", "PSMB1": "b6"}

FDR_MAX = 0.05
MIN_CELLS = 50
MIN_SPREAD = 0.5      # log2 units between the most- and least-changed subunit


def shard_urls(n):
    n = n or NSHARD
    return [f"{BASE}/train-{i:05d}-of-{NSHARD:05d}.parquet" for i in range(n)]


def fetch(n, threads, batch=8, retries=6):
    """Query shards in batches, authenticated, with exponential backoff.

    HF rate-limits hard (HTTP 429). Unauthenticated requests trip it almost
    immediately when DuckDB issues parallel range requests, so a token is
    effectively required. Set HF_TOKEN in the environment.
    """
    import duckdb
    import pandas as pd

    token = os.environ.get("HF_TOKEN") or os.environ.get(
        "HUGGING_FACE_HUB_TOKEN") or ""
    if not token:
        print("WARNING: no HF_TOKEN set -- expect HTTP 429. Create a free read "
              "token at huggingface.co/settings/tokens, then:\n"
              "    export HF_TOKEN=hf_...\n", flush=True)

    con = duckdb.connect()
    con.execute("INSTALL httpfs; LOAD httpfs;")
    con.execute(f"SET threads TO {threads};")
    if token:
        con.execute(f"""CREATE OR REPLACE SECRET hf (
                            TYPE HTTP,
                            EXTRA_HTTP_HEADERS MAP {{
                                'Authorization': 'Bearer {token}'
                            }}
                        );""")
        print("using HF_TOKEN for authentication", flush=True)

    urls = shard_urls(n)
    genes = "','".join(SUBUNIT)
    print(f"querying {len(urls)} shard(s), batch={batch}, "
          f"threads={threads} ...", flush=True)

    frames, t0 = [], time.time()
    for i in range(0, len(urls), batch):
        chunk = urls[i:i + batch]
        for attempt in range(retries):
            try:
                frames.append(con.execute(f"""
                    SELECT gene_name, baseMean, log2FoldChange, padj,
                           n_cells_trt, n_cells_ctrl,
                           drug, concentration, concentration_unit,
                           Cell_Name_Vevo, Cell_ID_Cellosaur, plate
                    FROM read_parquet({chunk!r})
                    WHERE gene_name IN ('{genes}')
                      AND log2FoldChange IS NOT NULL
                """).df())
                break
            except Exception as exc:                    # noqa: BLE001
                if "429" not in str(exc) or attempt == retries - 1:
                    raise
                wait = 2 ** attempt * 5
                print(f"  429 rate-limited; retry {attempt+1}/{retries} "
                      f"in {wait}s", flush=True)
                time.sleep(wait)
        done = min(i + batch, len(urls))
        el = time.time() - t0
        print(f"  {done}/{len(urls)} shards  {sum(len(f) for f in frames):,} rows"
              f"  {el:.0f}s", flush=True)

    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    dt = time.time() - t0
    print(f"\n  {len(df):,} rows in {dt:.1f}s", flush=True)
    if n and n < NSHARD:
        print(f"  EXTRAPOLATION: all {NSHARD} shards ~ "
              f"{dt*NSHARD/n/60:.0f} min", flush=True)
    return df


# ----------------------------------------------------------------- Jigsaw ---
def jigsaw(w):
    """(N_eff, supports) over the SIX admissible 20S orders of S14.2.

    The admissible set is six in every condition. Weighting changes p(omega),
    never Omega. Any table of these numbers must say so.
    """
    def pl(order):
        p, rem = 1.0, dict(w)
        for x in order:
            p *= rem[x] / sum(rem.values())
            del rem[x]
        return p
    ps = {o: pl(o) for o in itertools.permutations(["b1", "b5", "b6"])}
    H = -sum(p * math.log(p) for p in ps.values() if p > 0)
    s = lambda f: sum(p for o, p in ps.items() if f(o))
    return math.exp(H), {
        "13S+b1":    s(lambda o: o[0] == "b1"),
        "13S+b5":    s(lambda o: o[0] == "b5"),
        "13S+b1+b5": s(lambda o: set(o[:2]) == {"b1", "b5"}),
        "13S+b5+b6": s(lambda o: set(o[:2]) == {"b5", "b6"}),
    }


def screen(df):
    import pandas as pd
    df = df.copy()
    df["subunit"] = df["gene_name"].map(SUBUNIT)
    keys = ["Cell_Name_Vevo", "drug", "concentration", "concentration_unit", "plate"]

    out = []
    for key, g in df.groupby(keys, dropna=False):
        if set(g["subunit"]) != {"b1", "b5", "b6"}:
            continue                       # need all three measured
        if g["n_cells_trt"].min() < MIN_CELLS or g["n_cells_ctrl"].min() < MIN_CELLS:
            continue
        lfc = dict(zip(g["subunit"], g["log2FoldChange"]))
        fdr = dict(zip(g["subunit"], g["padj"]))
        hi, lo = max(lfc, key=lfc.get), min(lfc, key=lfc.get)
        spread = lfc[hi] - lfc[lo]
        rec = dict(zip(keys, key))
        rec.update({
            "lfc_b1": lfc["b1"], "lfc_b5": lfc["b5"], "lfc_b6": lfc["b6"],
            "padj_b1": fdr["b1"], "padj_b5": fdr["b5"], "padj_b6": fdr["b6"],
            "spread_log2": spread, "driver_hi": hi, "driver_lo": lo,
            # FDR required on the two subunits DRIVING the spread, not all three
            "fdr_ok": bool(pd.notna(fdr[hi]) and pd.notna(fdr[lo])
                           and fdr[hi] <= FDR_MAX and fdr[lo] <= FDR_MAX),
            "n_cells_trt": int(g["n_cells_trt"].min()),
        })
        w = {s_: 2 ** lfc[s_] for s_ in ("b1", "b5", "b6")}
        neff, sup = jigsaw(w)
        rec["N_eff"] = neff
        rec.update({f"support_{k}": v for k, v in sup.items()})
        out.append(rec)
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shards", type=int, default=1,
                    help="0 = all 1026; start with 1 to measure cost")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--batch", type=int, default=8,
                    help="shards per query; lower if you hit 429")
    ap.add_argument("--out", default="screen.csv")
    ap.add_argument("--drug-pattern", default="",
                    help="case-insensitive substring, e.g. bortezomib")
    a = ap.parse_args()

    df = fetch(a.shards, a.threads, a.batch)
    if df.empty:
        sys.exit("no rows -- check gene names against gene_metadata")

    res = screen(df)
    if res.empty:
        sys.exit("no conditions with all three subunits and adequate cell counts")

    res = res.sort_values("spread_log2", ascending=False)
    res.to_csv(a.out, index=False)

    q = res["spread_log2"].quantile([.5, .95, .99])
    print(f"\nconditions scored: {len(res):,}")
    print(f"spread_log2  median {q[.5]:.3f}   95th {q[.95]:.3f}   "
          f"99th {q[.99]:.3f}")
    print(f"N_eff        min {res['N_eff'].min():.3f}  "
          f"median {res['N_eff'].median():.3f}   (6.000 = no shift)")

    hits = res[(res.fdr_ok) & (res.spread_log2 >= MIN_SPREAD)]
    if a.drug_pattern:
        hits = hits[hits.drug.str.lower().str.contains(a.drug_pattern.lower())]
    print(f"\nhits (FDR<={FDR_MAX}, spread>={MIN_SPREAD}"
          f"{', drug~'+a.drug_pattern if a.drug_pattern else ''}): {len(hits)}")

    if hits.empty:
        print("\n*** NULL -- reportable. ***\n"
              "No condition moves the three subunits apart beyond threshold. "
              "That is exactly what Plackett-Luce scale-invariance predicts for "
              "a coordinated bounce-back, now shown on real data. Report the "
              "distribution above, not just the absence of hits.")
    else:
        cols = ["Cell_Name_Vevo", "drug", "concentration", "spread_log2",
                "N_eff", "support_13S+b5"]
        print(hits[cols].head(15).to_string(index=False))

    print(f"\nwrote {len(res)} -> {a.out}")
    print("REMINDER: the admissible set is SIX orders in every row. "
          "Weighting changes p(omega), never Omega.")
    os._exit(0)


if __name__ == "__main__":
    main()
