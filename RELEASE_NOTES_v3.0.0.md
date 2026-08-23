# scijigsaw v3.0.0 — evidence-weighted inference and hypothesis navigation

This is the release accompanying the PLOS manuscript *Scientific Jigsaw:
evidence-typed probabilistic reasoning and hypothesis navigation for
protein-complex assembly*. It adds a complete probabilistic layer and two
interactive interfaces on top of the exact enumerators of v2.x, which are
unchanged and reproduce all previously published counts.

## Highlights

**Weighted exact inference.** Strictly positive evidence potentials
redistribute probability among already-admissible mechanisms — soft
evidence can never create or delete a mechanism. Exact partition
functions, branch probabilities, weighted subcomplex support, entropy and
the effective mechanism number N_eff, all reducing to the uniform
enumerator at unit potentials. Includes the additivity lemma: a large
natural class of structural potentials is provably uninformative about
assembly order, and any interface potential is flat on an acyclic contact
graph.

**Probabilistic assembly tree + navigation.** TEST / WHY / COMPARE /
UPDATE and a `walk` operation over the exact mechanism space, in two
interfaces:

- `sj_session.py` — command-line session, biological-language commands,
  seed-anchored and merger-tree modes;
- `sj_visual.py` — a self-contained HTML page (no server): click a node
  to condition, select components to test an exact subcomplex, delete a
  relation to see killed/resurrected routes and the entropy change.

**Cryo-EM observation sets.** Jointly analyse structurally resolved
compositions as queries: detect mutually exclusive observations, compute
the minimum number of coexisting assembly routes, and attribute the
branching conclusion to specific depositions by leave-one-out.

**Counterfactual comparison.** Edit one relation, state or structural
hypothesis and receive the exact decomposition into killed, resurrected
and reweighted mechanisms, with the change in bits — including
per-contact in-silico interface deletion.

**Reproduction.** Complete biological encodings (VAMP2/SNARE, eIF3 with
the 10-model AlphaFold 3 ensemble, proteasome, NLRP3, 30S, mTORC1),
benchmark labels, analysis scripts and machine-readable outputs
reproducing every number and figure in the manuscript.

## Verification

- Weighted recurrences reproduce all uniform counts at unit potentials
  (Z = F, Z_B = G_B, machine precision).
- Tree-enumeration uniqueness certified by explicit construction and
  deduplication on all evaluated boards.
- Forward-branch calibration verified against brute-force enumeration to
  1e-12.

## Compatibility

Python 3.11–3.13. The deterministic core API is unchanged from v2.1.0.

## Citation

Please cite both the software (see CITATION.cff) and the accompanying
manuscript, referencing this release by its tag (`v3.0.0`). A Zenodo
archive with a DOI will be created upon acceptance of the manuscript and
linked here.
