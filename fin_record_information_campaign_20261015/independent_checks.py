"""Independent reconstruction checks for finite rational RI algebra (NOT for limit theorems).
Producer ri_exact is used only as an implementation under test; reference paths
use different direct combinatorics and direct labelled-copy refresh summation.
"""
from fractions import Fraction
from itertools import combinations
from math import comb
from pathlib import Path
from collections import defaultdict
import json
from ri_exact import hypergeom,joint_hypergeom,one_refresh_counts,count_refresh
checks=[]
def check(name,condition):
    checks.append({'name':name,'pass':bool(condition)})
    if not condition:raise AssertionError(name)

def direct_hypergeom(N,M,k,h):
    plus=set(range(M))
    count=0
    for group in combinations(range(N),k):
        if sum(a in plus for a in group)==h: count+=1
    return Fraction(count,comb(N,k))

def alt_joint(N,M,k1,k2,h1,h2):
    # sequential draws with independent combinatorial form, zero support guarded
    if h1<0 or h1>k1 or h2<0 or h2>k2:return Fraction(0)
    if h1>M or k1-h1>N-M:return Fraction(0)
    if h2>M-h1 or k2-h2>N-M-(k1-h1):return Fraction(0)
    return Fraction(comb(M,h1)*comb(N-M,k1-h1),comb(N,k1))*Fraction(comb(M-h1,h2)*comb(N-M-k1+h1,k2-h2),comb(N-k1,k2))

def direct_tag_jumps(n,b,q=12):
    # construct an ordered labelled population then refresh each old slot to every new label
    lab=list(b)
    for i,nn in enumerate(n):lab.extend([i]*(nn-b.count(i)))
    assert len(lab)==sum(n)
    out=defaultdict(Fraction)
    for a,old in enumerate(lab):
        for new in range(q):
            if new==old:continue
            lab2=lab.copy();lab2[a]=new
            cnt=tuple(lab2.count(i) for i in range(q))
            out[(cnt,tuple(lab2[:len(b)]))]+=Fraction(1,q)
    return dict(out)

for N in [2,3,4,8,12]:
 for M in sorted(set([0,1,N//2,N-1,N])):
  for k in range(N+1):
   for h in range(k+1):
    check('hypergeom_by_enumeration_%s_%s_%s_%s'%(N,M,k,h),hypergeom(N,M,k,h)==direct_hypergeom(N,M,k,h))
   check('hyg_norm_%s_%s_%s'%(N,M,k),sum((hypergeom(N,M,k,h) for h in range(k+1)),Fraction())==1)
check('impossible_specific',joint_hypergeom(12,2,4,3,3,0)==0)
for N in [3,4,6,10,12]:
 for M in sorted(set([0,1,N//2,N])):
  for k1 in range(N+1):
   for k2 in [0,N-k1]:
    check('joint_norm_%s_%s_%s_%s'%(N,M,k1,k2),sum((joint_hypergeom(N,M,k1,k2,h1,h2) for h1 in range(k1+1) for h2 in range(k2+1)),Fraction())==1)
    for h1 in range(k1+1):
     for h2 in range(k2+1):
      check('joint_sequential_%s_%s_%s_%s_%s_%s'%(N,M,k1,k2,h1,h2),joint_hypergeom(N,M,k1,k2,h1,h2)==alt_joint(N,M,k1,k2,h1,h2))
for n,b in [((1,)*12,(0,1)),((3,)+ (0,)*10+(2,),(0,11)),((2,2,1,1,1,1,1,1,1,1,1,1),(0,1))]:
 check('tag_g0_full_direct_%d'%len(checks),one_refresh_counts(n,b)==direct_tag_jumps(n,b))
 agg=defaultdict(Fraction)
 for (nn,_),r in direct_tag_jumps(n,b).items():agg[nn]+=r
 check('count_quotient_direct_%d'%len(checks),dict(agg)==count_refresh(n))
for N,M,k1,k2 in [(12,5,3,4),(12,6,2,3),(8,3,1,2),(20,11,4,5)]:
 p=Fraction(M,N)
 E1=Fraction(k1)*p;E2=Fraction(k2)*p
 E12=sum(Fraction(h1*h2)*joint_hypergeom(N,M,k1,k2,h1,h2) for h1 in range(k1+1) for h2 in range(k2+1))
 check('cross_cov_exact_%d'%N,E12-E1*E2==-Fraction(k1*k2)*p*(1-p)/Fraction(N-1))
path=Path(__file__).parent/'INDEPENDENT_CHECKS.json'
path.write_text(json.dumps({'checks':len(checks),'pass':sum(x['pass'] for x in checks),'fail':sum(not x['pass'] for x in checks),'scope':'independent direct finite rational reconstruction of producer finite identities; not limiting/filter proof','details':checks},indent=2)+'\n')
print('checks',len(checks),'passed',sum(x['pass'] for x in checks))
