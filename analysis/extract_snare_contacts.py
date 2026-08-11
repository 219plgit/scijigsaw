from Bio.PDB import PDBParser
import numpy as np, warnings
warnings.filterwarnings("ignore")
CUT,MINRES=5.0,3
def heavy(c):
    co=[];ri=[]
    for r in c:
        if r.id[0]!=" ": continue
        for a in r:
            if a.element!="H": co.append(a.coord); ri.append(r.id[1])
    return np.array(co),ri
def run(pdb,names,pairs):
    model=next(PDBParser(QUIET=True).get_structure(pdb,f"{pdb}.pdb").get_models())
    ch={c.id:c for c in model if c.id in names}
    at={k:heavy(v) for k,v in ch.items()}
    print(f"=== {pdb} ({len(ch)} chains) ===")
    for x,y in pairs:
        if x not in at or y not in at:
            print(f"   {names[x]}--{names[y]}: chain missing"); continue
        ca,ra=at[x]; cb,rb=at[y]; sa=set();sb=set()
        for i in range(len(ca)):
            d=np.linalg.norm(cb-ca[i],axis=1)
            if (d<=CUT).any():
                sa.add(ra[i])
                for j in np.where(d<=CUT)[0]: sb.add(rb[j])
        ok = len(sa)>=MINRES or len(sb)>=MINRES
        print(f"   {names[x]:14}--{names[y]:14} {len(sa):3}+{len(sb):3}  "
              f"{'QUALIFYING' if ok else 'below threshold'}")
N5={'A':'VAMP2','B':'Syntaxin-1A','C':'SNAP25-N','D':'SNAP25-C','E':'Complexin','F':'Syt-1'}
run('5W5C',N5,[('C','F'),('D','F'),('B','F'),('A','F'),('C','B'),('D','B')])
print()
N3={'A':'Munc18-1','B':'Syntaxin-1A'}
run('3C98',N3,[('A','B')])
