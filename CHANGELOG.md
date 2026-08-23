# Changelog
## v3.0.0

Release accompanying Liò & Liò, *Scientific Jigsaw: evidence-typed
probabilistic reasoning and hypothesis navigation for protein-complex
assembly* (PLOS, submitted). Major version: adds the evidence-weighted
inference layer and the interactive hypothesis-navigation interfaces on
top of the unchanged exact enumerators.

### Evidence-weighted exact inference
- Strictly positive mechanism weights `W(ω; E) > 0` over the admissible
  space only: soft evidence redistributes probability but can never
  create or delete a mechanism. At unit potentials every weighted count
  reproduces the uniform enumerator exactly (`Z = F`, `Z_B = G_B`).
- Weighted continuation recurrence `Q(S)` (seed-anchored) and weighted
  merger recurrence: exact partition functions, branch probabilities
  `Pr(p next | S, E)`, weighted subcomplex support, Shannon entropy
  `H₂` and effective mechanism number `N_eff = 2^H₂`.
- Interface-level structural potentials
  `φ_m(λ) = exp(λ(a·r̄_m + b·log(1+n_m)))` with λ-sweep reporting.
  Includes the additivity lemma check: edge-additive potentials are
  provably uninformative (all trees score identically), and any
  interface potential is flat on an acyclic contact graph.
- Forward-branch calibration verified against brute-force enumeration
  to 1e-12; tree-enumeration uniqueness certified by explicit
  construction and deduplication.
- Abundance-proportional encounter priors (GTEx v8/v10 demonstration on
  the fixed 252-history VAMP2 space) with per-tissue N_eff reporting.

### Probabilistic assembly tree and navigation API
- The probabilistic assembly tree as a navigational view of the exact
  (uniform or weighted) distribution: branch probabilities, hypothesis
  mass Pr(H | E) for declared mechanism families, conditioning on
  observed intermediates, branch-level entropy for locating
  high-uncertainty decision points.
- Navigation API exposing TEST / WHY / COMPARE / UPDATE and the `walk`
  operation over both representations, with conditioning trails and
  provenance preserved on every query.

### User interfaces
- `sj_session.py` — interactive command-line session; biological-language
  commands (`ask exact A B`, `ask before A B`, `see before A B`,
  `do drop-p A B`, `compare`, `walk`, `back`); supports seed-anchored
  and merger-tree modes without re-encoding.
- `sj_visual.py` — self-contained browser view; precomputes the factual
  space and one counterfactual per removable relation into a single
  HTML file; click-to-condition, exact-subcomplex testing and relation
  deletion with killed/resurrected routes and entropy change, no server.
  (Seed-anchored view; merger-tree visualisation is a named extension.)

### Cryo-EM observation-set analysis
- Joint analysis of resolved compositions entered as queries: pairwise
  mutual-exclusivity detection, minimum number of coexisting mechanisms
  by exhaustive covering, and leave-one-out attribution of the branching
  conclusion to specific depositions. Reproduces the two-route result on
  the six deposited yeast proteasome precursor structures.

### Counterfactual comparison
- Declared alternative models M′ with exact recomputation; differences
  decomposed into killed / resurrected / reweighted mechanisms with
  entropy change in bits; per-contact in-silico interface deletion
  (load-bearing vs redundant contacts).

### Encodings, data and reproduction
- Complete biological encodings shipped as boards: VAMP2/SNARE,
  eIF3 (declared + AF3 ensembles), proteasome late-stage, NLRP3, 30S,
  mTORC1 FRB alternative occupancy.
- AlphaFold 3 eIF3 model ensemble (10 predictions) and extracted
  contact sets included for reproduction.
- Machine-readable outputs and analysis scripts reproducing every
  number and figure in the manuscript (`reproduce_all.py`).

### Metadata
- CITATION.cff and pyproject.toml updated for the PLOS manuscript;
  Python 3.11–3.13 tested; stale v2.0.0 release URLs fixed. Zenodo
  archiving (and DOI insertion) is deferred to manuscript acceptance;
  .zenodo.json is kept current so the GitHub–Zenodo webhook archives the
  acceptance-time release correctly.

### Compatibility
- The deterministic core is unchanged; all v2.1.0 counts are
  reproduced exactly. The weighted layer is additive and reduces to the
  uniform enumerator at unit potentials.

## v2.1.0

- Typed relations: `contacts` and `prerequisites` declared separately, with
  provenance recorded per relation and used only to raise consistency warnings.
  A bare `requires` retains its previous meaning, so all published counts are
  unchanged.
- Assembly-tree enumeration (`n_trees`): counts binary merger histories over
  connected subsets, recovering the seed-anchored count exactly when every
  relation is a prerequisite and mergers are singleton additions to the seed.
  Verified against exhaustive enumeration on 2,156 connected graphs and 1,157
  cases carrying prerequisites.
- Exact experimental-subcomplex support (`support`, `support_table`):
  distinguishes impossible, optional and necessary intermediates. Verified
  against explicit tree enumeration on 6,005 subset queries.
- Diagnostics: `check_contact_graph` and `diagnose` report isolated units or a
  disconnected contact graph rather than a bare zero; cyclic prerequisites are
  rejected at construction.
- Declared prerequisites no longer create contacts: a temporal relation is not
  evidence of a physical interface.
- Interface extraction validated against PDBePISA across seven complexes
  (73 chain pairs, precision/recall/F1/kappa all 1.00).
- Legend rendering: section headings given clearance above; dashed
  alternative-occupancy row height corrected.
- New `analysis/` scripts: `reproduce_all.py`, `validate_support.py`,
  `extract_snare_contacts.py`, `af_ensemble.py`.
## v2.0.0

- Probabilistic layer (`contrib.probabilistic`) now reports P_competent,
  E[L], E[L/N], Var(L) and E[L | competent] via `MixtureResult.summary`;
  adds `order_distribution_joint` (declared correlated modifications) and
  `occupancy_uncertainty` (Beta-distributed occupancy, credible intervals).
- Counterfactual histories: `order_distribution_counterfactual` treats a
  declared prior over alternative constraint graphs as the S5.5 mixture
  (additive over exclusive histories); `linext_factorizes_disjoint`
  implements the exact disjoint-union identity
  L(P u Q) = C(n_P+n_Q, n_P) * L(P) * L(Q) for independent sub-posets.
- New analysis scripts: `analyse_disassembly.py` (dual-poset disassembly
  counts + counterfactual demo), `analyse_crowding.py` (q = L/n! as the
  probability a random arrival order is admissible), and
  `analyse_gtex_abundance.py` (abundance-weighted order distribution from
  user-provided GTEx TPM, with age-monotone overlap and donor-age
  stratification; GTEx data files are not shipped).
- `analyse_ordering_checks.py` adds the exact tier-label permutation test,
  complementary fractions, enrichment, edge concordance and Kendall tau-b.
- Renderer: removed the in-figure BRIDGE PIECES and ALTERNATIVE OCCUPANCY
  annotations (now described in figure captions) and tightened the board /
  alternative-occupancy panel spacing.

- **Two configurations.** The deterministic core and an extended (complete)
  configuration that adds optional extensions under `scijigsaw.contrib`.
- Moved the probabilistic layer to `scijigsaw.contrib.probabilistic`; the
  core no longer imports or exports it, so the deterministic tool is
  unchanged. Major version bump reflects the new capability and the
  reorganised import path.

## v1.6.1

Fixes the real-structure benchmark harness, which had never been executed.

- `benchmark/run_benchmark.py` looked for structures with the extensions
  `.cif/.mmcif/.pdf/.ent`; `.pdf` was a typo for `.pdb`, so downloaded PDB files
  were never found and every structure was skipped, yielding all-NaN results.
- It shelled out to `scijigsaw-extract` with `--min-residues`, which is not a
  flag of that command (`--min-interface-residues`), and passed a single
  multi-chain file where the extractor expects a directory of two-chain
  hub/partner complexes.
- The labels name proteins while structures contain chains, and no mapping
  existed. The harness now takes `--chain-map` (columns pdb,chain,protein) and
  provides `--list-chains`, which prints every chain with its header description
  so the mapping can be checked against the deposited file rather than assumed.
- Interface detection now calls `scijigsaw.extract.interface_residues` directly,
  the same function the extractor uses, over all chain combinations for each
  labelled protein pair, sweeping contact cutoff and minimum interface size.
- No library code changed; reported manuscript values are unaffected.
- The documented protein-table schema did not match the one the code read: the
  supplement listed `protein_id`, `functional_class` and `conservation_tier`
  while the reader required `name`, `function` and `age`, so a table written
  from the manuscript could not be loaded. Both spellings are now accepted,
  `conservation_tier`/`age` is genuinely optional as documented (it previously
  raised `KeyError`), and a table missing a required column now reports which
  column and what was found instead of a pandas error.
- Command-line entry points report missing files and malformed tables as a
  message rather than a traceback.
- Added `scijigsaw.probabilistic`: a probabilistic layer over declared
  modifications (Supplementary S5.4). Each modification carries a declared
  occupancy; `order_distribution` enumerates the 2^k instances (or samples
  for large k) and reuses the exact enumerator on each, returning the
  expected order count and P(target feasible). `cases.vamp2_phospho_mixture`
  reproduces the VAMP2 result (E=100.8 fusion-competent orders at p=0.6).
  Occupancies are declared inputs, never inferred.
- Legend ("How to read the pieces") gains an Optional-overlay section for
  declared modifications, shown as an extension point rather than a core
  channel; the legend canvas is slightly taller to fit it. render_legend
  now takes a `top` argument and legend() a `height` argument.
- Added `cases.vamp2_phospho_thr138`, a worked, evidence-declared example of a
  post-translational modification. SNAP25 Thr138 phosphorylation (PKA) is
  declared to remove the SNAP25 core interface; the fusion-competent state
  becomes infeasible and the board selects an AP180/CALM alternative-occupancy
  state. Modifications are represented from declared evidence, never inferred.
- Figure 3(A) now plots order-space reduction against poset DEPTH with component
  count shown as colour and a least-squares trend line, replacing a plot against
  component count that stacked ~130 points at each integer n and hid the depth
  signal. reproduce_numbers.py additionally reports the lazy-evaluation reach.
- Planarity is now checked on the CONNECTOR graph the renderer lays out, per
  feasible state, rather than on the precedence poset. These are different
  graphs: for the VAMP2 example the connector graph has 10 nodes and 12 edges
  (the printed tile set) against the poset's 8 and 10, because the alternative
  occupants AP180 and CALM are tiles but not enumerated units. New
  `Board.connector_graph`, `Board.exclusion_pairs`, `Board.feasible_states` and
  `Assembly.dependency_graph` give the script and the tests one construction
  path. Node and edge counts are pinned, so a dropped edge fails the test even
  though the graph would remain planar. The 30S dependency graph is reported
  both without the 16S seed (20/16) and with it (21/22).
- `scripts/reproduce_numbers.py` also reports the planarity of each encoded board,
  and the test suite asserts it. The manuscript states that one flat tile layer
  suffices for the assemblies presented; that holds only while their interaction
  graphs are planar, so the claim is now checked rather than asserted. This also
  puts the declared `networkx` dependency to use, which no module had imported.
- New `--variant backs` emits a reverse sheet for the tiles, carrying the protein
  name, its classification and, where the protein table provides them, a sequence
  accession and a structure identifier. The layout is mirrored horizontally so a
  long-edge duplex flip lands each reverse on its own tile. The three columns
  (`accession`, `pdb`, `class`) are optional additions to the protein table.
- Teacher tiles now label each connector with the partner protein it mates with.
  Previously a connector showed only a numbered, coloured badge, so the sheet said
  which connector it was but not where it went, and the pairing had to be looked
  up in the key table.
- The teacher key page printed "None/5" in the n/N column for boards without
  structure-derived coverage (that is, for every curated board, including the
  VAMP2 example). It now prints an em dash, matching the tiles themselves.
- Printable kits now carry the board name in the sheet header, so the teacher
  and student sets read "VAMP2 board \u2014 Scientific Jigsaw cut-out kit
  (teacher answer key / student class set)" rather than a generic title. The
  name is the protein with the most encoded interfaces and is used only when
  that maximum is unique; `--title` still overrides it.

## v1.6.0

Structural-meaning and enumeration release accompanying the revised manuscript.
All reported counts are unchanged.

- **Enumerator.** The subset recurrence is now evaluated lazily, with maximality
  pruning, so only states reachable under the precedence and exclusion rules (the
  order ideals of the poset) are allocated. The VAMP2 board visits 40 of 128
  states; random posets of depth n/2 at n=20 visit 174 of 1,048,576. The worst
  case is unchanged at O(n*2^n) time and O(2^n) memory -- a poset with few
  relations still reaches every subset -- so the n=25 guard remains.
- **Renderer: n/N is no longer derived from affinity.** Previously the displayed
  n/N was computed from `kd_nM`, which contradicted the documented meaning. It is
  now taken only from explicit structural columns (`coverage`, `coverage_denom`)
  and is left blank for curated relations, so the literature-curated VAMP2 board
  displays no n/N. Core-clique selection now uses functional homogeneity and the
  alternative-occupancy occupant is chosen by interface degree; neither uses
  affinity.
- **Extractor.** New `--min-interface-residues` filter (default 3), matching the
  documented minimum interface size. New output fields
  `hub_interface_residue_count`, `overlap_group_union_residue_count`,
  `structural_footprint_fraction`, `coverage` and `coverage_denom`, so the
  structural footprint fraction |I_j|/|S_C| is computed from the same residue
  sets as the overlap coefficient.
- **Tiles and legend.** "binding sites" restated as "encoded interfaces";
  connector geometry described as symbolic ("a key, not a molecular surface");
  the n/N row restated structurally; fixed a right-edge text-clipping bug.
- **Figures.** Accessory connectors now seat flush with the core, and the core
  carries the tab throughout (bridges and accessories receive sockets), in both
  the board renderer and the assembly-steps script. Assembly-step labels resized
  and nudged to remove name overlaps; the in-figure payoff line moved to the
  caption; the benchmark figure overlays the three encoded assemblies.
- **Documentation.** The overlap coefficient |I_A n I_B| / min(|I_A|,|I_B|) and
  the single-linkage partner clustering are documented to match the code.

## v1.5.0

Physical output: build the complex with scissors.

- New `scijigsaw-tiles` prints an **easy-to-cut set of tiles** (multi-page A4 PDF,
  or a single SVG) in which every protein is a separate tile. Each interface is a
  complementary **tab and socket with its own keyed shape** (round / wedge /
  square, four sizes), so a tab fits only its true partner — the kit is
  self-correcting, a *control of error* in cardboard.
- Two variants from the same cut geometry: a **student** set (protein name only —
  a puzzle) and a **teacher** answer key that adds the connector number, coverage
  n/N, and the precedence (bridge, *seat last*) and exclusion (dashed,
  *either/or*) cues that `scijigsaw-count` enumerates, with an assembly-key page.
- Topology, colour and evolutionary rings are reused from the renderer; only the
  layout differs. Tests include a page-collision check (overlapping pieces cannot
  be cut) and a check that the student set carries no solution labels.

Reproducibility and robustness release accompanying the manuscript.

- Every value reported in the paper is pinned in the test suite and regenerated
  by `scripts/reproduce_numbers.py`, which exits non-zero on any mismatch; every
  figure by `scripts/reproduce_figures.py`.
- Random-poset benchmark (`scijigsaw-bench`): 900 posets with collinearity
  diagnostics (variance-inflation factors), showing that sequentiality — not unit
  count — governs how much of the order space geometry eliminates.
- Structure-extractor robustness on real deposited files (selenomethionine,
  alternate locations, insertion codes, multi-model files), with parsing tests.
- 51 tests across Python 3.11–3.13 under continuous integration.

## v1.4.0

Applied the extractor to the deposited complexin–SNARE complex (PDB 1KIL) and corrected
the encoding it contradicted.

- `Complexin` bridges **VAMP2 and Syntaxin-1A**, not VAMP2 and SNAP25; it makes no contact
  with SNAP25. Added the direct VAMP2–Syntaxin edge (34 interface residues): the SNARE core
  is a four-helix bundle, not a chain.
- **Permitted orders 336 → 252 of 5,040** (93.3% → 95.0% eliminated).
- Core detection rewritten: a bundle is a clique, so the core is the most tightly bound
  clique, with a longest-path fallback for nucleated cascades.

## v1.3.0

- The 30S hierarchy is reported under both encodings. Conservative "any-of": 6×. The
  specific parents named by the published map: 45,545×. The inflammasome (10 subunits,
  depth 9) still prunes 40× more than the 30S (20 subunits, depth 3) — depth, not size.
- `RIBOSOME_30S_SPECIFIC` added, with a test pinning the contrast.

## v1.2.0

- The renderer now draws the paper's Figure 2 from `interactions.csv`. It previously
  produced only a hub-and-spoke board.
- Schema extended to partner–partner edges; backbone, bridges, pendants and contenders are
  derived from the table rather than declared.
- A shared site means competition only when the two partners do **not** bind each other; if
  they do, the piece bridges them.

## v1.1.0

Parser fixes for deposited structures — none of which synthetic test files can expose.

- A crystallographic B-factor is a temperature factor, not a confidence score. Reading it as
  pLDDT refused every experimental structure. Predicted models are now detected from the
  header (`--confidence auto|plddt|none`).
- Modified residues (MSE, SEP, TPO, PTR…) are HETATM but are amino acids; discarding all
  HETATM deleted them from interfaces.
- Insertion codes denote distinct residues (Kabat 52, 52A, 52B).

## v1.0.0

First public release.

- `scijigsaw-extract` — structures → interfaces, sites, overlaps. Contact cutoff, overlap
  threshold and confidence floor are user-set flags.
- `scijigsaw-render` — interaction table → jigsaw board.
- `scijigsaw-count` — exact linear-extension counting under the encoded seating model.
- `scijigsaw-bench` — 900-poset benchmark.
- `scripts/reproduce_numbers.py` regenerates every number in the paper and exits non-zero if
  the code and the manuscript disagree; `scripts/reproduce_figures.py` regenerates the figures.

## Known limitations

- The seating model requires all partners present before a piece is seated. A multivalent
  protein may bind one partner and meet a second later, so counts are a lower bound on
  biologically accessible pathways.
- Catalysis cannot be encoded: NSF/α-SNAP vacate the VAMP2 SNARE motif by an ATP-dependent
  process with no static tab-and-socket representation.
- Contact coverage (n/N) proxies buried interface extent, not affinity.
- The extractor has not been benchmarked against curated interface annotations.

## v1.5.0

- `scijigsaw.results.all_results()` — a single source of truth for every number the paper
  reports. The 336→252 correction had to be applied by hand in four documents and was missed
  in three; the manuscript, supplement, README and cover letter are now generated from one
  computed record.
- Test pinning the 1KIL structural correction, so it cannot silently regress.
