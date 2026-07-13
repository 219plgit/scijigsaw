# scijigsaw

**Interface geometry as a constraint on the assembly order of protein complexes.**

`scijigsaw` converts mapped protein interfaces into complementary tabs and sockets,
and counts *exactly* the assembly orders permitted by two geometric constraints:

- **precedence** — a piece may be seated only when the partners it binds are present,
  so a *bridge* piece spanning a seam (Complexin over VAMP2+SNAP25) cannot be seated
  until that seam is closed;
- **exclusion** — partners whose interfaces overlap cannot occupy the board together.

Permitted orders are the **linear extensions of the resulting constraint poset**,
counted by subset dynamic programming in `O(2^n · n)`.

## What it does not do

It **does not infer the assembly pathway.** It eliminates orders incompatible with the
*encoded interface and seating constraints*. Three limits are load-bearing, and stated
in the code as well as the paper:

1. **The seating model is not biology.** Precedence demands *all* partners be present
   before a piece is seated. A multivalent protein may bind one partner, remain
   partially engaged, and only then encounter a second. The count is a **lower bound**
   on the biologically accessible set.
2. **It cannot encode catalysis.** NSF/α-SNAP vacate the VAMP2 SNARE motif by an
   ATP-dependent process with no static tab-and-socket representation. Geometry bounds
   the admissible configurations; enzymes and kinetics select among them.
3. **Contact coverage (`n/N`) is not affinity.** It proxies the *extent* of interface
   engaged. A 5/5 electrostatically mismatched interface may bind more weakly than a
   2/5 hydrophobic one. `N` is a display parameter (default 5).

## Install

```bash
pip install -e ".[dev]"      # or: conda env create -f environment.yml
pytest -q                    # 18 tests
```

## Use

```bash
# TIER 2 -- render a board from an interaction table
scijigsaw-render examples/vamp2/proteins.csv examples/vamp2/interactions.csv \
    --out vamp2_board.svg

# TIER 1 -- derive that table from structures, then render it
scijigsaw-extract examples/structures --contact-cutoff 5.0 --out interactions.csv
scijigsaw-render  examples/vamp2/proteins.csv interactions.csv --out board.svg

# exact counting and the benchmark
scijigsaw-count all
scijigsaw-bench
```

```python
from scijigsaw import Assembly, VAMP2, INFLAMMASOME

VAMP2.summary()
# {'n': 7, 'depth': 3, 'total': 5040, 'permitted': 336, 'reduction': 15.0}

INFLAMMASOME.summary()
# {'n': 10, 'depth': 9, 'total': 3628800, 'permitted': 2, ...}
```

## The result

**Sequentiality, not size, governs how much geometry can eliminate.**

| assembly | units | depth | permitted / n! | reduction |
|---|---:|---:|---|---:|
| NLRP3 inflammasome | 10 | 9 | 2 / 3,628,800 | 1,814,400× |
| VAMP2 board | 7 | 3 | 336 / 5,040 | 15× |
| 30S ribosomal subunit | 20 | 3 | 4.2×10¹⁷ / 2.4×10¹⁸ | **6×** |

The 30S has **three times** the subunits of the inflammasome and prunes **six-fold**.
Across 900 random posets, `n` is a weak predictor of the reduction (r = 0.49) and depth
a strong one (r = 0.93); at fixed `n`, depth separates it by six orders of magnitude.
Depth, width and density are collinear, so they are alternative measures of one
property — *sequentiality* — not independent variables.

## Reproduce the paper

```bash
python scripts/reproduce_numbers.py    # ~25 s -- every number; exits non-zero on mismatch
python scripts/reproduce_figures.py    # ~40 s -- Figures 1-4 and S1, PDF + PNG
python scripts/benchmark_runtime.py    # OPTIONAL, slow: the O(2^n) runtime table (S7)
```

Every number in the manuscript is asserted in `reproduce_numbers.py` and pinned in
`tests/test_assembly.py`. **If the code and the paper disagree, the script fails and CI
goes red.** The runtime benchmark is deliberately separate: it measures an exponential
algorithm and traces allocations, so it is slow by construction and must not sit in the
default reproduction path.

Figure 2 is rendered **from the CSVs by the renderer**, not drawn by hand.

## Thresholds are choices

The extractor exposes `--contact-cutoff`, `--overlap` and `--min-confidence`. They are
not defaults to be ignored. In the bundled example the α-synuclein interface is **absent
at 5 Å and present at 8 Å** — the rendered claim depends on an analytical threshold, not
on the biology switching on. Report the value you used, and sweep it.

## Status of validation

The bundled structures in `examples/structures/` are **synthetic**, built with known
interfaces to verify that the extractor recovers what it is given. That is a correctness
check, **not** evidence of accuracy on real structures. Benchmarking against curated
interface annotations (precision/recall, site-clustering agreement, cutoff sensitivity)
across real complexes **remains to be done**, and is stated as such in the paper.

## Layout

```
src/scijigsaw/
  geometry.py    tabs, sockets, subsite ladders, contact coverage
  render.py      TIER 2: interaction table -> laid-out, labelled, saved board
  extract.py     TIER 1: structures -> interfaces, sites, overlaps
  assembly.py    the constraint poset; exact linear-extension counting
  cases.py       the three encoded assemblies
  benchmark.py   random-poset generator + regression/collinearity analysis
  cli.py         scijigsaw-render / -extract / -count / -bench
tests/           24 tests; the paper's numbers are pinned here
scripts/         reproduce_numbers.py, reproduce_figures.py, benchmark_runtime.py
examples/        VAMP2 input tables; synthetic test structures
```

[![tests](https://github.com/219plgit/scijigsaw/actions/workflows/test.yml/badge.svg)](https://github.com/219plgit/scijigsaw/actions)
[![licence](https://img.shields.io/badge/licence-MIT-blue.svg)](https://github.com/219plgit/scijigsaw/blob/main/LICENSE)

## Authors

Pietro Liò and Maria Teresa Liò
*Department of Computer Science and Technology, University of Cambridge, Cambridge, UK*

## Citation

See `CITATION.cff`. Please cite both the software (archived DOI) and the paper.

## Licence

MIT.
