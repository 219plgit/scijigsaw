from typed_assembly import TypedAssembly as TA
from itertools import combinations
from functools import lru_cache
import random

def brute(units, edges, pres=()):
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
    def mok(A,B):
        C=A|B; linked=False; m=A
        while m:
            i=(m&-m).bit_length()-1
            if con[i]&B: linked=True
            m&=m-1
        if not linked: return False
        m=C
        while m:
            i=(m&-m).bit_length()-1
            if pre[i]&~C: return False
            m&=m-1
        return True
    @lru_cache(maxsize=None)
    def F(S):
        if S&(S-1)==0: return 1
        tot=0; seen=set(); sub=(S-1)&S
        while sub:
            A,B=sub,S^sub
            k=frozenset((A,B))
            if A and B and k not in seen:
                seen.add(k)
                if conn(A) and conn(B) and mok(A,B): tot+=F(A)*F(B)
            sub=(sub-1)&S
        return tot
    return F((1<<n)-1)

random.seed(1); tested=0; bad=0
print("=== A2.3 F(S) vs brute force: exhaustive n=2..5, sampled n=6 ===")
for n in range(2,7):
    units=[f"u{i}" for i in range(n)]
    allp=list(combinations(units,2))
    for mask in range(1,1<<len(allp)):
        edges=[allp[i] for i in range(len(allp)) if mask>>i&1]
        if len(edges)<n-1: continue
        if n==6 and random.random()>0.05: continue
        A=TA(contacts=edges, units=units)
        if A.check_contact_graph(): continue      # skip disconnected
        f,b=A.n_trees(), brute(units,edges)
        tested+=1
        if f!=b: bad+=1; print(f"  MISMATCH n={n} {edges} fast={f} brute={b}")
print(f"  connected graphs tested: {tested}   mismatches: {bad}")

print("\n=== with prerequisites (n=2..5, sampled) ===")
tested2=bad2=0
for n in range(3,6):
    units=[f"u{i}" for i in range(n)]
    allp=list(combinations(units,2))
    for trial in range(400):
        k=random.randint(n-1,len(allp))
        edges=random.sample(allp,k)
        A0=TA(contacts=edges, units=units)
        if A0.check_contact_graph(): continue
        pres=[(u,v) for u,v in random.sample(edges,min(2,len(edges)))]
        try: A=TA(contacts=edges, prerequisites=pres, units=units)
        except ValueError: continue     # cyclic prerequisite, correctly rejected
        f,b=A.n_trees(), brute(units,edges,pres)
        tested2+=1
        if f!=b: bad2+=1; print(f"  MISMATCH n={n} edges={edges} pres={pres} {f} vs {b}")
print(f"  cases tested: {tested2}   mismatches: {bad2}")

print(f"\n=== RESULT: {'F(S) VALIDATED' if bad==0 and bad2==0 else 'FAILURES'} ===")
