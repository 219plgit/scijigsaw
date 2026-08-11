#!/usr/bin/env python3
"""
Probe Tahoe-100M WITHOUT the `datasets` streaming machinery.

`load_dataset(..., streaming=True)` resolves thousands of shard references
before yielding a row, and leaves non-daemon fsspec threads alive afterwards.
Both symptoms look like a hang. This script avoids it entirely: it lists the
repo, downloads only the small metadata Parquet files, and reads them with
pandas. Every step prints before it starts, so you always know where you are.

    pip install huggingface_hub pandas pyarrow

    python tahoe_probe.py --list
    python tahoe_probe.py --schema
    python tahoe_probe.py --drugs bortezomib carfilzomib MG-132 proteasome
"""
import argparse
import os
import sys

REPO = "tahoebio/Tahoe-100M"


def _hub():
    try:
        from huggingface_hub import HfApi, hf_hub_download
    except ImportError:
        sys.exit("pip install huggingface_hub pandas pyarrow")
    return HfApi(), hf_hub_download


def list_files(prefix=None):
    api, _ = _hub()
    print(f"[1/1] listing files in {REPO} ...", flush=True)
    files = api.list_repo_files(REPO, repo_type="dataset")
    meta = [f for f in files if not f.startswith("data/")]
    print(f"\n{len(files)} files total; {len(meta)} outside data/\n", flush=True)
    for f in sorted(meta):
        if prefix and not f.startswith(prefix):
            continue
        print("  " + f, flush=True)
    ndata = len(files) - len(meta)
    print(f"\n  data/  ({ndata} expression shards, not listed)", flush=True)
    return files


def grab(path):
    _, dl = _hub()
    print(f"  downloading {path} ...", flush=True)
    return dl(REPO, filename=path, repo_type="dataset")


def schema():
    import pandas as pd
    api, _ = _hub()
    files = api.list_repo_files(REPO, repo_type="dataset")

    targets = [f for f in files
               if f.endswith(".parquet") and not f.startswith("data/")]
    # one representative shard per logical table
    seen, picks = set(), []
    for f in sorted(targets):
        table = os.path.dirname(f) or os.path.basename(f)
        if table in seen:
            continue
        seen.add(table)
        picks.append(f)

    for f in picks:
        print(f"\n=== {f} ===", flush=True)
        try:
            df = pd.read_parquet(grab(f))
            print(f"  shape: {df.shape}", flush=True)
            print(f"  columns: {list(df.columns)}", flush=True)
            with pd.option_context("display.width", 200,
                                   "display.max_columns", 50):
                print(df.head(3).to_string()[:2000], flush=True)
        except Exception as exc:                       # noqa: BLE001
            print(f"  FAILED: {exc}", flush=True)
    print("\n=== schema probe complete ===", flush=True)


def drugs(patterns):
    import pandas as pd
    api, _ = _hub()
    files = api.list_repo_files(REPO, repo_type="dataset")
    cand = [f for f in files if "drug_metadata" in f and f.endswith(".parquet")]
    if not cand:
        sys.exit("no drug_metadata parquet found; run --list")
    df = pd.read_parquet(grab(cand[0]))
    print(f"\ncolumns: {list(df.columns)}\n", flush=True)

    pats = [p.lower() for p in patterns]
    blob = df.astype(str).agg(" ".join, axis=1).str.lower()
    hit = df[blob.apply(lambda b: any(p in b for p in pats))]
    if hit.empty:
        print("no matches -- widen the patterns", flush=True)
        return
    with pd.option_context("display.width", 250, "display.max_columns", 50,
                           "display.max_colwidth", 60):
        print(hit.to_string(), flush=True)
    print(f"\n{len(hit)} matching drugs. Copy the EXACT strings from the drug "
          f"column -- names carry formulation suffixes.", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--prefix")
    ap.add_argument("--schema", action="store_true")
    ap.add_argument("--drugs", nargs="+", metavar="PATTERN")
    a = ap.parse_args()

    if a.list:
        list_files(a.prefix)
    elif a.schema:
        schema()
    elif a.drugs:
        drugs(a.drugs)
    else:
        ap.print_help()
    os._exit(0)


if __name__ == "__main__":
    main()
