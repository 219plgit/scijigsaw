#!/usr/bin/env python3
"""Abundance-weighted assembly-order distribution for the VAMP2 fusion state,
grounded in GTEx tissue expression (v8 and v10). Runs locally; no data leaves
your machine.

For each sample: the 7 component TPMs form an abundance vector; the admissible
orders (252 of 5040) are weighted by an abundance-proportional sequential-
encounter model (at each step the next binder is drawn proportional to the
abundance of the still-admissible unplaced components). We then summarise the
weighted distribution per sample, aggregate to mean +/- variance within each
brain tissue, and compare across tissues and across GTEx releases.

CEILING: this is a declared, evidence-grounded PRIOR over admissible orders,
using expression as a measured proxy for the encounter-order prior. It is not a
kinetic model and not an independent validation of the framework.
"""
import gzip
import itertools
import math
import os
from collections import defaultdict

# ---------------------------------------------------------------- config
HERE = os.path.dirname(os.path.abspath(__file__))
V10_GCT = os.path.join(HERE, "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz")
V8_GCT  = os.path.join(HERE, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")
V10_ANN = os.path.join(HERE, "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt")
V8_ANN  = os.path.join(HERE, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
V10_PHENO = os.path.join(HERE, "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")
V8_PHENO  = os.path.join(HERE, "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")

# set to None to use ALL brain tissues; or list exact SMTSD labels to restrict
TISSUE_FILTER = None          # e.g. ["Brain - Cortex", "Brain - Hippocampus"]
INCLUDE_NERVE = False         # also include "Nerve - Tibial" if True

# Ensembl base IDs (version suffix stripped when matching)
GENES = {
    "SNAP25":  "ENSG00000132639",
    "STX1A":   "ENSG00000106089",
    "STXBP1":  "ENSG00000136854",
    "CPLX1":   "ENSG00000168993",
    "CPLX2":   "ENSG00000145920",
    "SYT1":    "ENSG00000067715",
    "SYP":     "ENSG00000102003",
    "SNCA":    "ENSG00000145335",
}
# model component  ->  how to build its abundance from GTEx symbols
COMPONENT_TO_SYMBOLS = {
    "SNAP25":        ["SNAP25"],
    "Syntaxin-1A":   ["STX1A"],
    "Munc18-1":      ["STXBP1"],
    "Complexin":     ["CPLX1", "CPLX2"],   # (c) summed
    "Syt-1":         ["SYT1"],
    "Synaptophysin": ["SYP"],
    "SNCA":          ["SNCA"],
}
SEED = "VAMP2"
COMPONENTS = list(COMPONENT_TO_SYMBOLS.keys())

# VAMP2 precedence (must match scijigsaw cases.VAMP2). Bridges Syt-1 & Complexin
# require the core seam; encode the same prerequisites used in the paper.
REQUIRES = {
    "Syntaxin-1A": {"SNAP25"},
    "Munc18-1":    {"Syntaxin-1A"},
    "Complexin":   {"Syntaxin-1A"},
    "Syt-1":       {"SNAP25", "Syntaxin-1A"},
}

# ---------------------------------------------------------------- helpers
def strip_ver(ensembl):
    return ensembl.split(".")[0]

def load_annotations(path):
    """Return {sample_id: SMTSD_tissue}."""
    tissue = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            smtsd = header.index("SMTSD")
        except ValueError:
            smtsd = 6  # fallback
        sampid = header.index("SAMPID") if "SAMPID" in header else 0
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > smtsd:
                tissue[f[sampid]] = f[smtsd]
    return tissue

def load_ages(path):
    """Return {subject_id: AGE_bracket} from SubjectPhenotypesDS."""
    ages = {}
    if not os.path.exists(path):
        return ages
    with open(path, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            ai = header.index("AGE"); si = header.index("SUBJID")
        except ValueError:
            return ages
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > ai:
                ages[f[si]] = f[ai]
    return ages

def subject_of(sample_id):
    parts = sample_id.split("-")
    return "-".join(parts[:2]) if len(parts) >= 2 else sample_id

def extract_genes(gct_path, cache_path):
    """Stream the big .gct once; cache a small TSV of the 9 genes x samples.
    Returns (sample_ids, {symbol: [tpm,...]})."""
    if os.path.exists(cache_path):
        with open(cache_path, encoding="utf-8") as fh:
            samples = fh.readline().rstrip("\n").split("\t")[1:]
            data = {}
            for line in fh:
                p = line.rstrip("\n").split("\t")
                data[p[0]] = [float(x) for x in p[1:]]
        return samples, data
    wanted = {v: k for k, v in GENES.items()}  # base ensembl -> symbol
    samples = None
    data = {}
    op = gzip.open
    with op(gct_path, "rt") as fh:
        fh.readline(); fh.readline()             # #1.2 ; dims
        header = fh.readline().rstrip("\n").split("\t")
        samples = header[2:]
        for line in fh:
            tab1 = line.find("\t")
            base = strip_ver(line[:tab1])
            if base in wanted:
                p = line.rstrip("\n").split("\t")
                data[wanted[base]] = [float(x) for x in p[2:]]
                if len(data) == len(GENES):
                    break
    with open(cache_path, "w", encoding="utf-8") as out:
        out.write("symbol\t" + "\t".join(samples) + "\n")
        for sym, vals in data.items():
            out.write(sym + "\t" + "\t".join(f"{v:.6g}" for v in vals) + "\n")
    return samples, data

def admissible_orders():
    ns = COMPONENTS
    out = []
    for p in itertools.permutations(ns):
        placed = {SEED}; ok = True
        for u in p:
            if not REQUIRES.get(u, set()) <= placed:
                ok = False; break
            placed.add(u)
        if ok:
            out.append(p)
    return out

ORDERS = admissible_orders()

AGE_TIER = {
    "SNAP25": 0, "Syntaxin-1A": 0, "Munc18-1": 0,
    "Complexin": 1, "Syt-1": 1, "Synaptophysin": 1,
    "SNCA": 2,
}

def _age_monotone(order):
    return all(AGE_TIER[order[i]] <= AGE_TIER[order[i+1]]
               for i in range(len(order) - 1))

AGE_MONOTONE_SET = {o for o in ORDERS if _age_monotone(o)}
UNIFORM_AGE_FRAC = len(AGE_MONOTONE_SET) / len(ORDERS)

def order_weight(order, abund):
    """Abundance-proportional sequential-encounter probability of one order."""
    placed = {SEED}; remaining = set(order); w = 1.0
    for u in order:
        cand = [c for c in remaining if REQUIRES.get(c, set()) <= placed]
        s = sum(abund[c] for c in cand)
        if s <= 0:
            return 0.0
        w *= abund[u] / s
        placed.add(u); remaining.discard(u)
    return w

def weighted_distribution(abund):
    ws = [order_weight(o, abund) for o in ORDERS]
    tot = sum(ws)
    if tot <= 0:
        return None
    ps = [w / tot for w in ws]
    # summaries: top-order mass, effective number of orders (exp entropy)
    top = max(ps)
    ent = -sum(p * math.log(p) for p in ps if p > 0)
    eff = math.exp(ent)               # effective # of orders (uniform => 252)
    top_idx = ps.index(top)
    age_mass = sum(ps[i] for i, o in enumerate(ORDERS) if o in AGE_MONOTONE_SET)
    return {"top_mass": top, "eff_orders": eff, "top_order": ORDERS[top_idx], "age_mass": age_mass, "ps": ps}

def sample_abundance(sym_tpm, sidx):
    a = {}
    for comp, syms in COMPONENT_TO_SYMBOLS.items():
        a[comp] = sum(sym_tpm[s][sidx] for s in syms if s in sym_tpm)
    return a

def analyse(gct, ann, cache, label):
    print(f"\n===== {label} =====")
    if not os.path.exists(gct):
        print(f"  (missing {gct} -- skipping)"); return
    tissue = load_annotations(ann)
    samples, sym_tpm = extract_genes(gct, cache)
    # group sample indices by tissue
    by_t = defaultdict(list)
    for i, s in enumerate(samples):
        t = tissue.get(s)
        if not t:
            continue
        keep = t.startswith("Brain") or (INCLUDE_NERVE and t.startswith("Nerve"))
        if TISSUE_FILTER is not None:
            keep = t in TISSUE_FILTER
        if keep:
            by_t[t].append(i)
    if not by_t:
        print("  no matching tissues found"); return
    all_eff, all_top = [], []
    for t in sorted(by_t):
        effs, ages = [], []
        for i in by_t[t]:
            d = weighted_distribution(sample_abundance(sym_tpm, i))
            if d:
                effs.append(d["eff_orders"]); ages.append(d["age_mass"])
        if not effs:
            continue
        n = len(effs)
        me = sum(effs)/n; ve = sum((x-me)**2 for x in effs)/n
        ma = sum(ages)/n; va = sum((x-ma)**2 for x in ages)/n
        print("  %-34s n=%4d  eff %5.1f+/-%4.1f (=%.2fx252)  age_mass %.3f+/-%.3f [unif %.3f]"
              % (t, n, me, ve**0.5, me/len(ORDERS), ma, va**0.5, UNIFORM_AGE_FRAC))
        all_eff.append(me); all_top.append(ma)
    if all_eff:
        gm = sum(all_eff)/len(all_eff)
        gv = sum((x-gm)**2 for x in all_eff)/len(all_eff)
        am = sum(all_top)/len(all_top)
        av = sum((x-am)**2 for x in all_top)/len(all_top)
        print("  ---- across %d tissues: eff_orders %.1f+/-%.1f (=%.2fx252)"
              % (len(all_eff), gm, gv**0.5, gm/len(ORDERS)))
        print("       age-monotone mass %.3f+/-%.3f vs uniform %.3f  ->  %.2fx enrichment"
              % (am, av**0.5, UNIFORM_AGE_FRAC, am/UNIFORM_AGE_FRAC))

def analyse_by_age(gct, ann, pheno, cache, label):
    """Stratify pooled brain samples by donor age bracket; report concentration
    and age-monotone enrichment per bracket."""
    print(f"\n===== {label}: by donor age (pooled brain tissues) =====")
    if not (os.path.exists(gct) and os.path.exists(pheno)):
        print("  (missing gct or phenotype file -- skipping)"); return
    tissue = load_annotations(ann)
    ages = load_ages(pheno)
    samples, sym_tpm = extract_genes(gct, cache)
    by_age = defaultdict(list)
    for i, sid in enumerate(samples):
        t = tissue.get(sid)
        if not t or not t.startswith("Brain"):
            continue
        a = ages.get(subject_of(sid))
        if a:
            by_age[a].append(i)
    for a in sorted(by_age):
        effs, agem = [], []
        for i in by_age[a]:
            d = weighted_distribution(sample_abundance(sym_tpm, i))
            if d:
                effs.append(d["eff_orders"]); agem.append(d["age_mass"])
        if not effs:
            continue
        n = len(effs)
        me = sum(effs)/n; ma = sum(agem)/n
        print("  age %-6s n=%4d  eff %5.1f (=%.2fx252)  age_mono_mass %.3f (enrich %.2fx vs %.3f)"
              % (a, n, me, me/len(ORDERS), ma, ma/UNIFORM_AGE_FRAC, UNIFORM_AGE_FRAC))


def main():
    print(f"admissible orders: {len(ORDERS)} (uniform baseline)")
    analyse(V10_GCT, V10_ANN, os.path.join(HERE, "vamp2_gtex_v10.tsv"), "GTEx v10")
    analyse(V8_GCT,  V8_ANN,  os.path.join(HERE, "vamp2_gtex_v8.tsv"),  "GTEx v8")
    analyse_by_age(V10_GCT, V10_ANN, V10_PHENO, os.path.join(HERE, "vamp2_gtex_v10.tsv"), "GTEx v10")
    analyse_by_age(V8_GCT,  V8_ANN,  V8_PHENO,  os.path.join(HERE, "vamp2_gtex_v8.tsv"),  "GTEx v8")
    print("\nNote: eff_orders is the effective number of assembly orders under the")
    print("abundance-weighted distribution (exp of entropy); the uniform count is")
    print(f"{len(ORDERS)}. Lower means abundance concentrates the distribution onto")
    print("fewer orders. This is a declared abundance PRIOR (expression proxy),")
    print("not a kinetic prediction.")

if __name__ == "__main__":
    main()
