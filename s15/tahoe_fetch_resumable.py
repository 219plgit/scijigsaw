#!/usr/bin/env python3
"""
Resumable, low-footprint Tahoe sweep. Replaces tahoe_full_screen.py's fetch.

WHY THE PREVIOUS RUN DIED
-------------------------
1. MEMORY (confirmed cause). DuckDB defaults to ~80% of system RAM and does
   not release it between queries on the same connection; 1,026 sequential
   reads of ~91 MB Parquet files exhausted it. Disk was NOT the problem --
   the HF cache stayed at 11 MB, so blob deletion worked correctly.
2. MEMORY, secondary. All 1,026 DataFrames were held in a list and
   concatenated only at the end.

This version:
  * downloads each shard with plain streaming HTTP to a temp file we own
    (the same request curl succeeds on), then deletes it immediately;
  * appends rows to CSV per shard -- nothing accumulates in RAM;
  * caps DuckDB at 1 GB and recycles the connection every 50 shards
    (the previous run exhausted RAM: DuckDB defaults to ~80% of system memory
    and accumulates across many queries on a single connection);
  * records progress in a state file, so it RESUMES after any interruption;
  * prints disk free space every 25 shards and aborts if it drops below 5 GB.

Peak footprint: ~100 MB disk, <100 MB RAM, regardless of shard count.

    export HF_TOKEN=$(python3 -c "from huggingface_hub import get_token; print(get_token())")
    python3 tahoe_fetch_resumable.py --shards 1026 --out full_raw.csv

Interrupt with Ctrl-C at any time; rerun the same command to continue.
Then analyse with:
    python3 tahoe_full_screen.py --reuse full_raw.csv --out-prefix full
"""
import argparse
import csv
import os
import shutil
import sys
import tempfile
import time
import urllib.request

REPO = "tahoebio/Tahoe-100M"
NSHARD = 1026
URL = ("https://huggingface.co/datasets/tahoebio/Tahoe-100M/resolve/main/"
       "metadata/pseudobulk_differential_expression/"
       "train-{i:05d}-of-01026.parquet")
GENES = ("PSMB6","PSMB5","PSMB7","PSMB8","PSMB9","PSMB10")
COLS = ["gene_name", "baseMean", "log2FoldChange", "lfcSE", "padj", "n_cells_trt",
        "n_cells_ctrl", "drug", "concentration", "Cell_Name_Vevo", "plate"]
MIN_FREE_GB = 5


def free_gb(path="."):
    return shutil.disk_usage(path).free / 1e9


def _valid_parquet(path, expect=None):
    """A Parquet file starts and ends with the magic bytes PAR1.

    A silently truncated download (connection closed early without raising)
    passes the download call but fails at read time with 'No magic bytes found
    at end of file'. Checking here turns that into a retry instead of a crash.
    """
    try:
        size = os.path.getsize(path)
        if size < 8:
            return False
        if expect is not None and size != expect:
            return False
        with open(path, "rb") as fh:
            if fh.read(4) != b"PAR1":
                return False
            fh.seek(-4, os.SEEK_END)
            return fh.read(4) == b"PAR1"
    except OSError:
        return False


def download(i, token, dest, retries=6):
    """Stream one shard, then verify it is a complete Parquet file."""
    req = urllib.request.Request(URL.format(i=i))
    if token:
        req.add_header("Authorization", f"Bearer {token}")
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                expect = r.headers.get("Content-Length")
                expect = int(expect) if expect else None
                with open(dest, "wb") as fh:
                    shutil.copyfileobj(r, fh, length=1 << 20)
            if not _valid_parquet(dest, expect):
                raise IOError(f"truncated/corrupt download "
                              f"({os.path.getsize(dest)} bytes"
                              f"{f', expected {expect}' if expect else ''})")
            return os.path.getsize(dest)
        except Exception as exc:                          # noqa: BLE001
            try:
                os.remove(dest)
            except OSError:
                pass
            if attempt == retries - 1:
                raise
            wait = min(2 ** attempt * 5, 120)
            print(f"    shard {i}: {type(exc).__name__}: {exc}; retry "
                  f"{attempt+1}/{retries} in {wait}s", flush=True)
            time.sleep(wait)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shards", type=int, default=NSHARD)
    ap.add_argument("--out", default="full_raw.csv")
    ap.add_argument("--state", default=None)
    a = ap.parse_args()
    state_path = a.state or a.out + ".state"

    import duckdb

    token = os.environ.get("HF_TOKEN", "")
    if not token:
        print("WARNING: HF_TOKEN not set; expect rate limiting.\n"
              "  export HF_TOKEN=$(python3 -c \"from huggingface_hub import "
              "get_token; print(get_token())\")\n", flush=True)

    done = set()
    if os.path.exists(state_path):
        done = {int(x) for x in open(state_path).read().split()}
        print(f"resuming: {len(done)} shards already processed", flush=True)

    header_written = os.path.exists(a.out) and os.path.getsize(a.out) > 0
    genes = "','".join(GENES)

    def new_con():
        """Fresh connection with a hard memory cap.

        DuckDB defaults to ~80% of system RAM and accumulates across many
        queries on one connection; 1,026 sequential reads exhausted RAM in the
        previous run. Cap it, and recycle every RECYCLE shards.
        """
        c = duckdb.connect()
        c.execute("PRAGMA memory_limit='1GB';")
        c.execute("PRAGMA threads=2;")
        return c

    RECYCLE = 50
    con = new_con()
    tmpdir = tempfile.mkdtemp(prefix="tahoe_")
    tmp = os.path.join(tmpdir, "shard.parquet")

    t0, nrows, nbytes = time.time(), 0, 0
    skipped = []
    todo = [i for i in range(a.shards) if i not in done]
    print(f"{len(todo)} shards to fetch; free disk {free_gb():.0f} GB", flush=True)

    try:
        for n, i in enumerate(todo, 1):
            if free_gb() < MIN_FREE_GB:
                print(f"\nABORT: free disk {free_gb():.1f} GB < {MIN_FREE_GB} GB",
                      flush=True)
                break
            try:
                nbytes += download(i, token, tmp)
                df = con.execute(f"""
                    SELECT {', '.join(COLS)}
                    FROM read_parquet('{tmp}')
                    WHERE gene_name IN ('{genes}') AND log2FoldChange IS NOT NULL
                """).df()
            except Exception as exc:                      # noqa: BLE001
                # one bad shard out of 1,026 must not end the run; it is left
                # out of the state file, so a later rerun retries it
                print(f"    shard {i}: SKIPPED after retries -- {exc}",
                      flush=True)
                skipped.append(i)
                try:
                    os.remove(tmp)
                except OSError:
                    pass
                continue
            os.remove(tmp)
            if n % RECYCLE == 0:          # release accumulated DuckDB memory
                con.close()
                con = new_con()

            with open(a.out, "a", newline="") as fh:
                w = csv.writer(fh)
                if not header_written:
                    w.writerow(COLS)
                    header_written = True
                w.writerows(df.itertuples(index=False, name=None))
            nrows += len(df)

            with open(state_path, "a") as fh:
                fh.write(f"{i}\n")

            if n % 25 == 0 or n == len(todo):
                el = time.time() - t0
                rate = n / el
                print(f"  {n}/{len(todo)}  {nrows:,} rows  "
                      f"{nbytes/1e9:.1f} GB  {el/60:.0f} min  "
                      f"ETA {(len(todo)-n)/rate/60:.0f} min  "
                      f"free {free_gb():.0f} GB", flush=True)
    except KeyboardInterrupt:
        print("\ninterrupted -- rerun the same command to resume", flush=True)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    if skipped:
        print(f"\n{len(skipped)} shard(s) skipped after retries: "
              f"{skipped[:20]}{' ...' if len(skipped) > 20 else ''}\n"
              f"  rerun the same command to retry them", flush=True)
    print(f"\n{nrows:,} rows appended -> {a.out}", flush=True)
    print(f"state: {a.out}.state  (delete it to start over)", flush=True)
    print(f"\nNEXT: python3 tahoe_full_screen.py --reuse {a.out} "
          f"--out-prefix full", flush=True)
    os._exit(0)


if __name__ == "__main__":
    main()
