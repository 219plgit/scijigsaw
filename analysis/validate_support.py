import importlib.util, random
from itertools import combinations
spec=importlib.util.spec_from_file_location("a","assembly_support.py")
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); A=m.Assembly

def brute_all_trees(units, edges, pres=()):
    """Enumerate actual trees as nested frozensets; return list of node-sets."""
    idx={u:i for i,u in enumerate(units)}; n=len(units)
    con=[0]*n; pre=[0]*n
    for u,v in edges: con[idx[u]]|=1<<idx[v]; con[idx[v]]|=1<<idx[u]
    for u,v in pres: pre[idx[v]]|=1<<idx[u]
    def conn(S):
        if S==0: return False
        f=(S&-S).bit_length()-1; seen=1<<f; st=[f]
        while st:
            i=st.pop(); nb=con[i]&S&~seen
            while nb:
                j=(nb&-nb).bit_length()-1; seen|=1<<j; st.append(j); nb&=nb-1
        return seen==S
    def mok(Aa,Bb):
        C=Aa|Bb; linked=False; mm=Aa
        while mm:
            i=(mm&-mm).bit_length()-1
            if con[i]&Bb: linked=True
            mm&=mm-1
        if not linked: return False
        mm=C
        while mm:
            i=(mm&-mm).bit_length()-1
            if pre[i]&~C: return False
            mm&=mm-1
        return True
    from functools import lru_cache
    @lru_cache(maxsize=None)
    def trees(S):
        """return tuple of frozensets-of-nodes, one per distinct tree"""
        if S&(S-1)==0: return (frozenset([S]),)
        out=[]; seen=set(); sub=(S-1)&S
        while sub:
            Aa,Bb=sub,S^sub; k=frozenset((Aa,Bb))
            if Aa and Bb and k not in seen:
                seen.add(k)
                if conn(Aa) and conn(Bb) and mok(Aa,Bb):
                    for ta in trees(Aa):
                        for tb in trees(Bb):
                            out.append(frozenset(ta|tb|{S}))
            sub=(sub-1)&S
        return tuple(out)
    return trees((1<<n)-1), idx

random.seed(3); tested=bad=0
print("=== support() vs explicit tree enumeration ===")
for n in range(3,6):
    units=[f"u{i}" for i in range(n)]; allp=list(combinations(units,2))
    for trial in range(120):
        edges=random.sample(allp, random.randint(n-1,len(allp)))
        a0=A(contacts=edges, units=units)
        if a0.check_contact_graph(): continue
        pres=[(u,v) for u,v in random.sample(edges,min(1,len(edges)))]
        try: a=A(contacts=edges, prerequisites=pres, units=units)
        except ValueError: continue
        if a.check_contact_graph(): continue
        allt, idx = brute_all_trees(units, edges, pres)
        if not allt: continue
        for k in range(1,n+1):
            for c in combinations(units,k):
                E=sum(1<<idx[u] for u in c)
                truth=sum(1 for t in allt if E in t)
                res=a.support(list(c))
                tested+=1
                if res["n_containing"]!=truth or res["n_total"]!=len(allt):
                    bad+=1
                    if bad<4:
                        print(f"  MISMATCH n={n} E={c} got={res['n_containing']}/{res['n_total']} truth={truth}/{len(allt)}")
print(f"  {tested} subset queries tested, {bad} mismatches")
print("RESULT:", "VALIDATED" if bad==0 else "FAILURES")
