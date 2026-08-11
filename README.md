# scijigsaw

**Interface geometry as a constraint on the assembly order of protein complexes.**

`scijigsaw` converts declared protein interfaces into complementary tabs and sockets and
counts *exactly* the assembly histories the evidence permits. It does not predict a
pathway: it enumerates the space of pathways that a declared model does not exclude, and
names the relation responsible when something is excluded.

## Two enumerators

**Seed-anchored orders.** Components join one growing assembly containing a declared
seed. Permitted orders are the linear extensions of the constraint poset, counted by
subset dynamic programming in `O(2^n · n)`.

**Assembly trees.** Independently formed subcomplexes merge. A subcomplex is a connected
subset of the contact graph; an event is a binary merger. This is the general case, and
it reduces exactly to the seed-anchored count when every relation is a prerequisite and
every merger adds a single component to the seed-containing assembly.

The second exists because the first cannot represent a subcomplex that forms away from
the seed — a limitation found, not assumed, when the yeast eIF3 `g:i` dimer proved
unbuildable.

## Typed relations

Relations are declared with an explicit type, and never inferred from graph shape.

- **contact** — undirected; two components can form an interface. No claim about order.
- **prerequisite** — directed; one component must precede another.
- **exclusion** — symmetric; two components compete for a surface.
- **observed subcomplex** — an annotation, not a constraint: an experimentally supported
  subset, queried against the enumerated space rather than converted into a temporal claim.

Each relation carries its evidence provenance, used to raise consistency warnings and
never to reassign a type. A relation derived from a structure but declared as a temporal
prerequisite is flagged for review.

A bare `requires` retains its original meaning — every listed partner is a prerequisite —
so encodings and counts published under the earlier model are unchanged.

```python
from scijigsaw import Assembly, VAMP2, INFLAMMASOME

VAMP2.summary()
# {'n': 7, 'depth': 3, 'total': 5040, 'permitted': 252, 'reduction': 20.0}

VAMP2.n_trees()          # 252 under legacy semantics
VAMP2.support(["SNAP25", "Syntaxin-1A"])
# {'n_containing': 0, 'n_total': 252, 'support': 0.0, 'reading': 'impossible'}
```

`support()` reports the exact fraction of admissible histories containing a subset,
distinguishing **impossible** (0), **optional** (0 < s < 1) and **necessary** (1).

## What it does not do

It **does not infer the assembly pathway.** It eliminates histories incompatible with the
*declared* interface and seating constraints. Four limits are load-bearing, and stated in
the code as well as the paper:

1. **The seating model is not biology.** Precedence demands *all* partners be present
   before a piece is seated. A multivalent protein may bind one partner, remain partially
   engaged, and only then encounter a second. Counts are exact for the declared model but
   are not formal bounds in either direction: kinetics may select a subset of the admitted
   routes, and an over-restrictive encoding may omit accessible ones.
2. **It cannot encode catalysis.** NSF/α-SNAP vacate the VAMP2 SNARE motif by an
   ATP-dependent process with no static tab-and-socket representation. Geometry bounds the
   admissible configurations; enzymes and kinetics select among them.
3. **Contact coverage (`n/N`) is not affinity.** It proxies the *extent* of interface
   engaged. A 5/5 electrostatically mismatched interface may bind more weakly than a 2/5
   hydrophobic one. `N` is a display parameter (default 5).
4. **Tree enumeration needs contacts, not just precedence.** A classical assembly map
   records which components depend on which for incorporation, not which touch. Boards
   encoded that way have a disconnected contact graph and admit no merger history; the
   30S subunit is such a case, and the software reports the responsible units rather than
   a bare zero.

## Install

```bash
pip install -e ".[dev]"      # or: conda env create -f environment.yml
pytest -q                    # 71 tests
```

## Reproduce the paper

```bash
python analysis/reproduce_all.py       # every reported number, each with its expected value
python scripts/reproduce_numbers.py    # ~25 s -- exits non-zero on mismatch
python scripts/reproduce_figures.py    # ~40 s -- figures, PDF + PNG
python scripts/benchmark_runtime.py    # OPTIONAL, slow: the O(2^n) runtime table
```

Every number in the manuscript is asserted in these scripts and pinned in
`tests/test_assembly.py`. **If the code and the paper disagree, the script fails and CI
goes red.** The runtime benchmark is deliberately separate: it measures an exponential
algorithm and traces allocations, so it must not sit in the default reproduction path.

Supporting analyses under `analysis/`:

| script | what it does |
|---|---|
| `reproduce_all.py` | regenerates every reported number |
| `validate_support.py` | brute-force validation of the support recurrence (6,005 subset queries) |
| `extract_snare_contacts.py` | interface extraction from 5W5C and 3C98 |
| `af_ensemble.py` | AlphaFold ensemble analysis: interface recurrence, robustness, support |

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

## The result

**Sequentiality, not size, governs how much geometry can eliminate.**

| assembly | units | depth | permitted / n! | reduction |
|---|---:|---:|---|---:|
| NLRP3 inflammasome | 10 | 9 | 2 / 3,628,800 | 1,814,400× |
| VAMP2 board | 7 | 3 | 252 / 5,040 | 20× |
| 30S, conservative tiers | 20 | 3 | 4.2×10¹⁷ / 2.4×10¹⁸ | **6×** |
| 30S, specific dependencies | 20 | 3 | 5.3×10¹³ / 2.4×10¹⁸ | 45,545× |

The 30S has **three times** the subunits of the inflammasome and, under the conservative
tiered encoding, prunes **six-fold**. Across 900 random posets the logarithmic reduction
`log10(n!/L)` correlates weakly with `n` (r = 0.49) and strongly with depth (r = 0.93);
at fixed `n`, depth separates it by six orders of magnitude. Depth, width and density are
collinear — alternative measures of one property, *sequentiality*, not independent
variables.

## Validation against real structures

The extractor was compared with **PDBePISA** across seven heteromeric complexes from
distinct structural and functional families — PDB **1KIL** (SNARE/complexin), **1K8K**
(Arp2/3), **1BXR** (carbamoyl-phosphate synthetase), **1WBJ** (tryptophan synthase),
**2NS1** (AmtB–GlnK), **1I7Q** (anthranilate synthase) and **1D09** (aspartate
transcarbamoylase). Taking PISA's identity-symmetry interfaces (≥ 200 Å²) as reference,
agreement was exact: **precision, recall, F1 and Cohen's κ all 1.00 over 73 chain pairs**
(36 interfaces, 37 non-interfaces). A cutoff-sensitivity sweep is reported for 1KIL and
1K8K.

The agreement is not a threshold coincidence. In several of these entries a crystal-packing
contact is *larger* than a genuine interface (1990 vs 493 Å² in 2NS1; 1624 vs 1363 Å² in
1WBJ), so buried area alone cannot separate the classes. The extractor is correct on them
because it works within the deposited assembly and never encounters the symmetry mates
that generate such contacts. This is a consistency benchmark against a reference method,
not an independent biological ground truth.

**The board was corrected by structure.** Applying the extractor to the deposited
complexin–SNARE structure (1KIL) revised the literature-derived encoding: complexin's
central helix contacts synaptobrevin (13 residues) and syntaxin (10) and **not** SNAP-25
at the declared threshold. Correcting it tightened the admissible set from **336 to 252**
orders. Individual declared relations were likewise sourced by extraction — synaptotagmin
contacts from **5W5C**, Munc18-1 from **3C98**, proteasome β-subunit contacts from
**9RM0** and **9RLA**.

Two edges cannot be settled this way. AP180/CALM competing for VAMP2's SNARE motif has no
deposited complex known to us, and α-synuclein is intrinsically disordered, so its VAMP2
N-terminal interface is not crystallographically resolved. Those remain literature-encoded
and are labelled as such.

## Held-out biological tests

Three, on two independent systems, exercising two distinct failure modes.

**yeast eIF3** — the seed-anchored model recovers two of the three subcomplexes seen by
native mass spectrometry; the third contains no seed and cannot be represented at all.
Typed relations plus merger trees recover all three (6 orders → 14 trees).

**neuronal SNARE** — an external hold-out. The SNAP-25:syntaxin acceptor complex, observed
in solution and formed without synaptobrevin, has support **0/252** under the seed-anchored
board and **284/1,242** under merger trees. The converse error is also reported: the model
admits a synaptobrevin:syntaxin pair that native MS did not capture.

**20S proteasome** — a temporal hold-out and a *retrospective* diagnosis. The long-assumed
sequential β-subunit order admits one history and declares two subsequently resolved cryo-EM
intermediates impossible; the structure-supported encoding admits all of them. Here the
representation was adequate and the fault was a temporal assumption.

## Thresholds are choices

The extractor exposes `--contact-cutoff`, `--overlap` and `--min-confidence`. They are not
defaults to be ignored. In the bundled example the α-synuclein interface is **absent at
5 Å and present at 8 Å** — the rendered claim depends on an analytical threshold, not on
the biology switching on. Report the value you used, and sweep it.

## Reading real structures

Three things that every real structure has, and no synthetic test file does. Each broke
the extractor; each is now tested in `tests/test_structure_formats.py`.

- **A B-factor is not a confidence score.** In an experimental structure it is a
  temperature factor (~28 Å² for a well-ordered crystal). Read as pLDDT, that scores below
  any sane confidence floor, and the tool would **refuse every PDB entry ever deposited**.
  Predicted models are detected from the header (`--confidence auto`); experimental
  structures are never refused for low B. Override with `--confidence plddt|none`.
- **Modified residues are HETATM but are still amino acids.** Selenomethionine (MSE) and
  phospho-Ser/Thr/Tyr sit in interfaces all the time. Skipping all HETATM silently deletes
  them: six MSE residues vanished from a sixteen-residue interface.
- **Insertion codes are distinct residues.** Kabat-numbered antibodies carry 52, 52A, 52B.
  Keying on the residue number alone conflates them.

A fourth, learned the hard way: **mmCIF `label_asym_id` and `auth_asym_id` are different
labels.** Mapping one onto the other silently dropped three chains from a proteasome
extraction and produced a confident wrong answer. `analysis/` scripts derive chain
identity from entity records and verify with a second parser.

## Print a physical kit

`scijigsaw-render` draws the board already assembled; `scijigsaw-tiles` does the opposite.
It lays every protein out as a **separate tile** on A4 pages, so you can print (at 100%),
cut along the solid outlines, and build the complex by hand.

Each interface is a complementary **tab (protrudes) and socket (indents) with its own keyed
shape** — round / wedge / square, in four sizes — so a tab fits *only* its true partner. A
piece that will not fit signals an interface it does not match: the kit is self-correcting.

Two sets come out of the same cut geometry:

- **student** — each tile carries only the protein name. Learners build the complex by
  matching connector shapes and reasoning about the biology.
- **teacher** — the answer key: adds the connector number on every tab/socket, the
  interface contact coverage *n/N*, and the precedence (bridge, *seat last*) and exclusion
  (dashed, *either/or*) cues, plus an assembly-key page.

```bash
scijigsaw-tiles examples/vamp2/proteins.csv examples/vamp2/interactions.csv \
    --out vamp2_kit.pdf --variant both      # writes _student.pdf and _teacher.pdf
```

Output is vector and true-to-scale (`.pdf` multi-page A4, or `.svg` single sheet). The
tiles are generated by the same renderer, from the same interaction table, as the reported
counts — so revising the evidence revises the artefact. Educational effectiveness has not
been evaluated and no student results are reported.

## Optional extensions

The core never imports `contrib`, so the deterministic tool runs unchanged whether or not
the extensions are used.

```python
from scijigsaw.contrib.probabilistic import Modification, order_distribution
```

## Layout

```
src/scijigsaw/
  geometry.py    tabs, sockets, subsite ladders, contact coverage
  render.py      TIER 2: interaction table -> laid-out, labelled, saved board
  extract.py     TIER 1: structures -> interfaces, sites, overlaps
  assembly.py    typed relations; both enumerators; exact subcomplex support
  cases.py       the encoded assemblies
  benchmark.py   random-poset generator + regression/collinearity analysis
  cli.py         scijigsaw-render / -extract / -count / -bench / -tiles
tests/           71 tests; the paper's numbers are pinned here
scripts/         reproduce_numbers.py, reproduce_figures.py, benchmark_runtime.py
analysis/        reproduce_all.py, validate_support.py, extract_snare_contacts.py,
                 af_ensemble.py
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
> protein-assembly order*, v2.1.0.
> https://github.com/219plgit/scijigsaw/releases/tag/v2.1.0

See `CITATION.cff`.

## Licence

MIT.
