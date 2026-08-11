"""Reproduce every number reported in the manuscript.

Run from a checkout with scijigsaw installed:
    python reproduce_all.py
"""
from scijigsaw.assembly import Assembly as A
from itertools import permutations

def line(): print("-"*66)

# ---------------------------------------------------------------- boards
VAMP2_REQ = {"SNAP25":{"VAMP2"},"Syntaxin-1A":{"VAMP2","SNAP25"},
             "Munc18-1":{"Syntaxin-1A"},"Complexin":{"VAMP2","Syntaxin-1A"},
             "Syt-1":{"SNAP25","Syntaxin-1A"},"Synaptophysin":{"VAMP2"},
             "SNCA":{"VAMP2"}}
VAMP2_EXC = [({"SNAP25"},{"AP180","CALM"})]
NLRP3_REQ = {"NLRP3":set(),"NEK7":{"NLRP3"},"NLRP3_disk":{"NLRP3","NEK7"},
             "NLRP3_PYDfil":{"NLRP3_disk"},"ASC_PYDfil":{"NLRP3_PYDfil"},
             "ASC_CARDfil":{"ASC_PYDfil"},"Casp1_CARDfil":{"ASC_CARDfil"},
             "Casp1_cat":{"Casp1_CARDfil"},"GSDMD":{"Casp1_cat"},"IL1B":{"Casp1_cat"}}
EIF3_REQ  = {"Prt1":set(),"Tif34":{"Prt1"},"Tif35":{"Tif34"},"Tif32":{"Prt1"},
             "Nip1":{"Tif32"}}
R30_REQ   = {"S4":set(),"S7":set(),"S8":set(),"S15":set(),"S17":set(),"S20":set(),
             "S16":{"S4"},"S5":{"S4","S8"},"S12":{"S4"},"S6":{"S8","S15"},
             "S18":{"S6"},"S9":{"S7"},"S19":{"S7"},"S13":{"S7"},"S2":{"S5"},
             "S3":{"S5"},"S10":{"S9"},"S14":{"S19"},"S11":{"S6"},"S21":{"S13"}}

print("== seed-anchored order counts ==")
for n,(kw,exp) in {
  "VAMP2": (dict(requires=VAMP2_REQ, excludes=VAMP2_EXC, seed="VAMP2"), 252),
  "NLRP3": (dict(requires=NLRP3_REQ), 2),
  "eIF3":  (dict(requires=EIF3_REQ, seed="Prt1"), 6),
  "30S":   (dict(requires=R30_REQ, seed="16S"), 53417729664000),
}.items():
    g = A(**kw).n_orders_permitted()
    print(f"   {n:6} {g:>18,}  expect {exp:>18,}  {'OK' if g==exp else 'FAIL'}")

line(); print("== tree enumerator recovers the order count (contact-complete boards) ==")
for n,kw,exp in [("VAMP2",dict(requires=VAMP2_REQ,excludes=VAMP2_EXC,seed="VAMP2"),252),
                 ("NLRP3",dict(requires=NLRP3_REQ),2),
                 ("eIF3", dict(requires=EIF3_REQ,seed="Prt1"),6)]:
    g=A(**kw).n_trees(); print(f"   {n:6} n_trees={g:>6}  expect {exp}  {'OK' if g==exp else 'FAIL'}")

line(); print("== 30S: disconnected contact graph ==")
R=A(requires=R30_REQ, seed="16S")
print(f"   n_trees = {R.n_trees()}  (expect 0)")
print(f"   {R.check_contact_graph()}")

line(); print("== eIF3 exact support ==")
OBS=[{"Tif35","Tif34"},{"Prt1","Tif34"},{"Prt1","Tif35","Tif34"}]
for nm,kw in [("seed-anchored", dict(requires=EIF3_REQ,seed="Prt1",observed_subcomplexes=OBS)),
              ("typed",dict(contacts=[("Prt1","Tif34"),("Tif34","Tif35"),
                                      ("Prt1","Tif32"),("Tif32","Nip1")],
                            prerequisites=[("Tif34","Tif35")],seed="Prt1",
                            observed_subcomplexes=OBS))]:
    a=A(**kw); print(f"   {nm} [{a.n_trees()} trees]")
    for r in a.support_table():
        print(f"      {'+'.join(r['subset']):26} {r['n_containing']}/{r['n_total']}  {r['reading']}")

line(); print("== eIF3 ablation ==")
seed="Prt1"; ns=["Tif34","Tif35","Tif32","Nip1"]; req={"Tif34":set(),"Tif35":{"Tif34"},
                                                        "Tif32":set(),"Nip1":set()}
def ok(o):
    p={seed}
    for u in o:
        if not req[u]<=p: return False
        p.add(u)
    return True
print(f"   typed single-addition orders: {sum(1 for p in permutations(ns) if ok(p))} (expect 12)")

line(); print("== SNARE: acceptor complex and sensitivity ==")
CON=[("VAMP2","SNAP25"),("VAMP2","Syntaxin-1A"),("VAMP2","Complexin"),
     ("SNAP25","Syntaxin-1A"),("SNAP25","Syt-1"),("Syntaxin-1A","Complexin"),
     ("Syntaxin-1A","Syt-1"),("Syntaxin-1A","Munc18-1"),
     ("VAMP2","Synaptophysin"),("VAMP2","SNCA")]
SOBS=[{"SNAP25","Syntaxin-1A"},{"SNAP25","VAMP2"},
      {"VAMP2","SNAP25","Syntaxin-1A"},{"VAMP2","SNAP25","Syntaxin-1A","Complexin"}]
V=A(requires=VAMP2_REQ, excludes=VAMP2_EXC, seed="VAMP2", observed_subcomplexes=SOBS)
print(f"   seed-anchored [{V.n_orders_permitted()} orders]")
for r in V.support_table():
    print(f"      {'+'.join(r['subset']):40} {r['n_containing']}/{r['n_total']}  {r['reading']}")
CPLX=[("VAMP2","Complexin"),("Syntaxin-1A","Complexin")]
SYT=[("SNAP25","Syt-1"),("Syntaxin-1A","Syt-1")]; MUN=[("Syntaxin-1A","Munc18-1")]
for lab,pre,exp in [("all prerequisites",CPLX+SYT+MUN,1242),
                    ("complexin omitted",SYT+MUN,2272),
                    ("Syt-1 omitted",CPLX+MUN,2750),
                    ("contacts only",[],5412)]:
    T=A(contacts=CON,prerequisites=pre,excludes=VAMP2_EXC,seed="VAMP2")
    r=T.support(["SNAP25","Syntaxin-1A"])
    print(f"   {lab:20} {T.n_trees():>6,} trees (expect {exp:,})  "
          f"acceptor {r['n_containing']}/{r['n_total']} = {r['support']:.3f}")

line(); print("== proteasome ==")
POBS=[{"S","b1"},{"S","b5"},{"S","b1","b5"},{"S","b5","b6"}]
for lab,req_ in [("canonical sequential",{"S":set(),"b5":{"S"},"b6":{"b5"},"b1":{"b5","b6"}}),
                 ("structure-supported", {"S":set(),"b1":{"S"},"b5":{"S"},"b6":{"S"}})]:
    a=A(requires=req_, seed="S", observed_subcomplexes=POBS)
    row=" ".join(f"{'+'.join(sorted(r['subset']))}={r['n_containing']}/{r['n_total']}"
                 for r in a.support_table())
    print(f"   {lab:22} {a.n_orders_permitted()} orders  {row}")
line(); print("done")
