from typed_assembly import TypedAssembly as TA

BOARDS = {
 "VAMP2": (dict(requires={"SNAP25":{"VAMP2"},"Syntaxin-1A":{"VAMP2","SNAP25"},
    "Munc18-1":{"Syntaxin-1A"},"Complexin":{"VAMP2","Syntaxin-1A"},
    "Syt-1":{"SNAP25","Syntaxin-1A"},"Synaptophysin":{"VAMP2"},"SNCA":{"VAMP2"}},
    excludes=[({"SNAP25"},{"AP180","CALM"})], seed="VAMP2"), 252),
 "NLRP3": (dict(requires={"NLRP3":set(),"NEK7":{"NLRP3"},
    "NLRP3_disk":{"NLRP3","NEK7"},"NLRP3_PYDfil":{"NLRP3_disk"},
    "ASC_PYDfil":{"NLRP3_PYDfil"},"ASC_CARDfil":{"ASC_PYDfil"},
    "Casp1_CARDfil":{"ASC_CARDfil"},"Casp1_cat":{"Casp1_CARDfil"},
    "GSDMD":{"Casp1_cat"},"IL1B":{"Casp1_cat"}}), 2),
 "eIF3": (dict(requires={"Prt1":set(),"Tif34":{"Prt1"},"Tif35":{"Tif34"},
    "Tif32":{"Prt1"},"Nip1":{"Tif32"}}, seed="Prt1"), 6),
 "30S": (dict(requires={"S4":set(),"S7":set(),"S8":set(),"S15":set(),
    "S17":set(),"S20":set(),"S16":{"S4"},"S5":{"S4","S8"},"S12":{"S4"},
    "S6":{"S8","S15"},"S18":{"S6"},"S9":{"S7"},"S19":{"S7"},"S13":{"S7"},
    "S2":{"S5"},"S3":{"S5"},"S10":{"S9"},"S14":{"S19"},"S11":{"S6"},
    "S21":{"S13"}}, seed="16S"), 53417729664000),
}

ok = True
print("=== A2.2 published counts preserved ===")
for nm, (kw, exp) in BOARDS.items():
    got = TA(**kw).n_orders_permitted()
    good = got == exp; ok &= good
    print(f"  {nm:7} {got:>18,}  expect {exp:>18,}  {'OK' if good else 'FAIL'}")

print("\n=== A2.7 typed form == legacy form (VAMP2) ===")
legacy = TA(**BOARDS["VAMP2"][0])
pre = [(p, u) for u, ps in BOARDS["VAMP2"][0]["requires"].items() for p in ps]
typed = TA(prerequisites=pre, excludes=[({"SNAP25"},{"AP180","CALM"})],
           seed="VAMP2")
a, b = legacy.n_orders_permitted(), typed.n_orders_permitted()
print(f"  legacy={a}  typed={b}  {'OK' if a==b else 'FAIL'}")
ok &= (a == b)

print("\n=== A2.5 cyclic prerequisite must raise ===")
try:
    TA(requires={"A":{"B"},"B":{"A"}}, seed="S")
    print("  FAIL: no error raised"); ok = False
except ValueError as e:
    print(f"  OK: {e}")

print("\n=== A2.6 disconnected contact graph reported ===")
d = TA(**BOARDS["30S"][0])
print(f"  30S: {d.check_contact_graph()}")
print(f"  n_trees = {d.n_trees()}")

print("\n=== eIF3 typed: all three native-MS subcomplexes ===")
E = TA(contacts=[("Prt1","Tif34"),("Tif34","Tif35"),("Prt1","Tif32"),
                 ("Tif32","Nip1")],
       prerequisites=[("Tif34","Tif35")], seed="Prt1",
       observed_subcomplexes=[{"Tif35","Tif34"},{"Prt1","Tif34"},
                              {"Prt1","Tif35","Tif34"}])
print(f"  trees = {E.n_trees()}")
for e, v in E.coverage().items():
    print(f"    {'+'.join(sorted(e)):26} formable: {v}")

print(f"\n=== {'ALL PASS' if ok else 'FAILURES PRESENT'} ===")
