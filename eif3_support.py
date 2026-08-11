from scijigsaw.assembly import Assembly as A

OBS = [{'Tif35','Tif34'}, {'Prt1','Tif34'}, {'Prt1','Tif35','Tif34'}]

seed_anchored = A(requires={'Prt1': set(), 'Tif34': {'Prt1'},
                            'Tif35': {'Tif34'}, 'Tif32': {'Prt1'},
                            'Nip1': {'Tif32'}},
                  seed='Prt1', observed_subcomplexes=OBS)

typed = A(contacts=[('Prt1','Tif34'), ('Tif34','Tif35'),
                    ('Prt1','Tif32'), ('Tif32','Nip1')],
          prerequisites=[('Tif34','Tif35')], seed='Prt1',
          observed_subcomplexes=OBS)

for name, a in [('seed-anchored', seed_anchored), ('typed', typed)]:
    print(f"\n{name}: {a.n_trees()} trees")
    for r in a.support_table():
        print("   %-24s %d/%d  %.3f  %s"
              % ('+'.join(r['subset']), r['n_containing'],
                 r['n_total'], r['support'], r['reading']))
