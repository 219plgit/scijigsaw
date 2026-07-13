#!/usr/bin/env python3
"""OPTIONAL. Runtime and memory of the exact count (Supplement S7).

Slow by construction: it measures an O(2^n n) algorithm and traces allocations.
It is NOT part of `reproduce_numbers.py`, which must finish quickly.

    python scripts/benchmark_runtime.py --max-n 20 --memory
"""
import argparse, random, time
from scijigsaw.benchmark import random_poset

ap = argparse.ArgumentParser()
ap.add_argument("--max-n", type=int, default=20)
ap.add_argument("--memory", action="store_true",
                help="also trace peak memory (adds large overhead)")
a = ap.parse_args()

rng = random.Random(1)
print(f"{'n':>4}{'time (s)':>11}" + (f"{'peak MB':>10}" if a.memory else ""))
for n in range(10, a.max_n + 1, 2):
    P = random_poset(n, max(2, n // 3), 0.25, rng)
    if a.memory:
        import tracemalloc
        tracemalloc.start()
    t0 = time.perf_counter()
    P.n_orders_permitted()
    dt = time.perf_counter() - t0
    if a.memory:
        _, peak = tracemalloc.get_traced_memory(); tracemalloc.stop()
        print(f"{n:>4}{dt:>11.2f}{peak/1e6:>10.1f}")
    else:
        print(f"{n:>4}{dt:>11.2f}")
print("\nDoubling per unit, as O(2^n n) requires; n = 25 is the practical ceiling.")
