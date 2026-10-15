"""Independent exact finite combinatorial / tagged-count algebra for RI-001.
No FIN-specific A7 values are simulated here; g=0 rates are exact Fractions.
"""
from fractions import Fraction
from itertools import product
from math import comb

def hypergeom(N,M,k,h):
    if not all(isinstance(x,int) for x in (N,M,k,h)) or N<0 or not 0<=M<=N or not 0<=k<=N:
        raise ValueError('invalid N/M/k/h')
    if h<0 or h>k or h>M or k-h>N-M:return Fraction(0)
    return Fraction(comb(M,h)*comb(N-M,k-h),comb(N,k))

def joint_hypergeom(N,M,k1,k2,h1,h2):
    if not all(isinstance(x,int) for x in (N,M,k1,k2,h1,h2)) or N<0 or not 0<=M<=N or min(k1,k2)<0 or k1+k2>N:
        raise ValueError('invalid')
    if h1<0 or h1>k1 or h2<0 or h2>k2 or h1+h2>M or (k1-h1)+(k2-h2)>N-M:
        return Fraction(0)
    return Fraction(comb(M,h1)*comb(N-M,k1-h1)*comb(M-h1,h2)*comb(N-M-k1+h1,k2-h2),comb(N,k1)*comb(N-k1,k2))

def one_refresh_counts(n,tagged,q=12):
    """g=0 exact off-diagonal generator for (counts, tagged labels), q=12.
    Tuple result -> Fraction rate. tagged is ordered. Every copy attempts at rate one.
    """
    assert len(n)==q and sum(n)>0
    c=[tagged.count(i) for i in range(q)]
    if any(c[i]>n[i] for i in range(q)):raise ValueError('incompatible tag/count')
    outs={}
    def put(i,j,ntag):
        nn=list(n);nn[i]-=1;nn[j]+=1
        key=(tuple(nn),tuple(ntag))
        outs[key]=outs.get(key,Fraction(0))+Fraction(1,q)
    for a,i in enumerate(tagged):
        for j in range(q):
            if j!=i:
                b=list(tagged);b[a]=j
                put(i,j,b)
    for i in range(q):
        for j in range(q):
            if i!=j and n[i]-c[i]>0:
                nn=list(n);nn[i]-=1;nn[j]+=1
                key=(tuple(nn),tuple(tagged))
                outs[key]=outs.get(key,Fraction(0))+Fraction(n[i]-c[i],q)
    return outs

def count_refresh(n,q=12):
    outs={}
    for i in range(q):
        for j in range(q):
            if i!=j and n[i]:
                nn=list(n);nn[i]-=1;nn[j]+=1
                key=tuple(nn)
                outs[key]=outs.get(key,Fraction(0))+Fraction(n[i],q)
    return outs

def verify_example(n,b):
    q=len(n)
    a=one_refresh_counts(n,b,q)
    c=count_refresh(n,q)
    agg={}
    for (nn,_),r in a.items():agg[nn]=agg.get(nn,Fraction(0))+r
    assert agg==c
    assert sum(a.values())==Fraction(sum(n)*(q-1),q)
    return True

if __name__=='__main__':
    print(hypergeom(12,2,4,3),joint_hypergeom(12,2,4,3,3,0))
