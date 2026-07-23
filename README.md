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
pytest -q                    # 51 tests
```

## Use

```bash
# TIER 2 -- render a board from an interaction table
scijigsaw-render examples/vamp2/proteins.csv examples/vamp2/interactions.csv \
    --out vamp2_board.svg

# TIER 1 -- derive that table from structures, then render it
scijigsaw-extract examples/structures --contact-cutoff 5.0 --out interactions.csv
scijigsaw-render  examples/vamp2/proteins.csv interactions.csv --out board.svg

# PRINT & BUILD -- an easy-to-cut set of tiles you assemble by hand
scijigsaw-tiles examples/vamp2/proteins.csv examples/vamp2/interactions.csv \
    --out vamp2_kit.pdf

# exact counting and the benchmark
scijigsaw-count all
scijigsaw-bench
```

## Print a physical kit

`scijigsaw-render` draws the board already assembled; `scijigsaw-tiles` does the
opposite. It lays every protein out as a **separate tile** on A4 pages, so you
can print (at 100%), cut along the solid outlines, and build the complex by hand.
Each interface is a complementary **tab (protrudes) and socket (indents) with its
own keyed shape** — round / wedge / square, in four sizes — so a tab fits *only*
its true partner. A piece that will not fit signals an interface it does not
match: the kit is self-correcting.

Two sets come out of the same cut geometry, so a solved model and a class set are
the same physical pieces:

- **student** — each tile carries only the protein name. Learners build the
  complex by matching connector shapes and reasoning about the biology.
- **teacher** — the answer key: adds the connector number on every tab/socket,
  the interface contact coverage *n/N*, and the precedence (bridge, *seat last*)
  and exclusion (dashed, *either/or*) cues, plus an assembly-key page.

```bash
scijigsaw-tiles examples/vamp2/proteins.csv examples/vamp2/interactions.csv \
    --out vamp2_kit.pdf --variant both      # writes _student.pdf and _teacher.pdf
```

Output is vector and true-to-scale (`.pdf` multi-page A4, or `.svg` single sheet).

```python
from scijigsaw import Assembly, VAMP2, INFLAMMASOME

VAMP2.summary()
# {'n': 7, 'depth': 3, 'total': 5040, 'permitted': 252, 'reduction': 20.0}

INFLAMMASOME.summary()
# {'n': 10, 'depth': 9, 'total': 3628800, 'permitted': 2, ...}
```

## The result

**Sequentiality, not size, governs how much geometry can eliminate.**

| assembly | units | depth | permitted / n! | reduction |
|---|---:|---:|---|---:|
| NLRP3 inflammasome | 10 | 9 | 2 / 3,628,800 | 1,814,400× |
| VAMP2 board | 7 | 3 | 252 / 5,040 | 20× |
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

## Reading real structures

Three things that every real structure has, and no synthetic test file does. Each broke
the extractor; each is now tested in `tests/test_structure_formats.py`.

- **A B-factor is not a confidence score.** In an experimental structure it is a
  temperature factor (~28 Å² for a well-ordered crystal). Read as pLDDT, that scores
  below any sane confidence floor, and the tool would **refuse every PDB entry ever
  deposited**. Predicted models are now detected from the header (`--confidence auto`),
  and experimental structures are never refused for low B. Override with
  `--confidence plddt|none`.
- **Modified residues are HETATM but are still amino acids.** Selenomethionine (MSE)
  and phospho-Ser/Thr/Tyr sit in interfaces all the time. Skipping all HETATM silently
  deletes them: six MSE residues vanished from a sixteen-residue interface.
- **Insertion codes are distinct residues.** Kabat-numbered antibodies carry 52, 52A,
  52B. Keying on the residue number alone conflates them.

## Status of validation

**The boards in the paper are hand-encodings of published interface assignments, not
structure-derived results.** The extractor has never been run on the deposited structures
of the assemblies it depicts.

`scripts/validate_board_from_structures.py` is the test that would change that: it applies
the extractor to the SNARE structures themselves (1SFC, 3C98, 1KIL, 5CCG) and asks whether
the derived interfaces reproduce the board — whether VAMP2's SNARE motif really does contact
SNAP25 and syntaxin, whether Munc18-1 really does contact syntaxin and *not* VAMP2, whether
complexin really does touch two chains at once and is therefore a bridge. **That test has not
been run.** Run it before believing the board.

Two edges cannot be settled this way at all: AP180/CALM competing for VAMP2's SNARE motif has
no deposited complex known to us, and α-synuclein is intrinsically disordered, so its VAMP2
N-terminal interface is not crystallographically resolved. Those remain literature-encoded.

## Parser robustness (a lesser thing)

The bundled structures are **synthetic**, built with known interfaces to verify that the
extractor recovers what it is given, plus the three realistic stress cases above. That is
a correctness check, **not** evidence of accuracy on real structures. Benchmarking against
curated interface annotations (precision/recall, site-clustering agreement, cutoff
sensitivity) across experimental complexes **remains to be done**, and the paper says so.

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
tests/           51 tests; the paper's numbers are pinned here
scripts/         reproduce_numbers.py, reproduce_figures.py, benchmark_runtime.py
examples/        VAMP2 input tables; synthetic test structures
```

[![tests](https://github.com/219plgit/scijigsaw/actions/workflows/test.yml/badge.svg)](https://github.com/219plgit/scijigsaw/actions)
[![licence](https://img.shields.io/badge/licence-MIT-blue.svg)](https://github.com/219plgit/scijigsaw/blob/main/LICENSE)

## Authors

Pietro Liò and Maria Teresa Liò
*Department of Computer Science and Technology, University of Cambridge, Cambridge, UK*

## Citation

Cite the tagged release and the paper:

> Liò, P. and Liò, M.T. (2026) *scijigsaw: interface geometry as a constraint on
> protein-assembly order*, v1.6.0.
> https://github.com/219plgit/scijigsaw/releases/tag/v1.6.0

See `CITATION.cff`. A Zenodo DOI can be minted later by enabling the Zenodo–GitHub
hook and cutting a new release; nothing in the code needs to change.

## Licence

MIT.
