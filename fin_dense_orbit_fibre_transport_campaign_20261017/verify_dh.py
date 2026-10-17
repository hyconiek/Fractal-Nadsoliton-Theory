#!/usr/bin/env python3
"""Bounded exact tests ONLY; does not certify asymptotics, local Gibbs, or all-Hall."""
import math,itertools,random,json
from fractions import Fraction

def action(t,rot,reflect):
 return tuple(t[((-j if reflect else j)+rot)%12] for j in range(12))

def orbit(t):
 return sorted(set(action(t,rot,ref) for rot in range(0,12,2) for ref in (False,True)))

def minimum_pair_distance(items):
 return min(sum(abs(a-b) for a,b in zip(x,y)) for i,x in enumerate(items) for y in items[i+1:])

def centered(n_even,n_odd,k):
 vals=[0]*12
 for p,K in [(0,n_even),(1,n_odd)]:
  b,r=divmod(K,6)
  for i in range(6):vals[2*i+p]=b+int(i<r)+k*(2*i-5)
 return tuple(vals)

nchecks=0
base=tuple(v for i in range(6) for v in (2*i-5,2*i-5))
assert len(orbit(base))==12;nchecks+=1
assert minimum_pair_distance(orbit(base))==36;nchecks+=1
assert sum(base[::2])==sum(base[1::2])==0;nchecks+=1
# Symbolic integer Fourier C/S and orbit distances on central/critical-count fibres.
for N in (12000,24000,120000,1200000,12000000):
 k=math.isqrt(N) if False else math.floor(N**(1/3))
 for h in (0,math.floor(N**.75/4),-math.floor(N**.75/4)):
  n=centered(N//2+h,N//2-h,k)
  assert sum(n)==N and min(n)>0;nchecks+=1
  assert sum(n[::2])-sum(n[1::2])==2*h;nchecks+=1
  O=orbit(n)
  assert len(O)==12;nchecks+=1
  assert minimum_pair_distance(O)>=36*k-12;nchecks+=1
  C=sum(n[j]*(1 if j%4==0 else -1 if j%4==2 else 0) for j in range(12))
  S=sum(n[j]*(1 if j%4==1 else -1 if j%4==3 else 0) for j in range(12))
  assert abs(C)<=6*k+3 and abs(S)<=6*k+3;nchecks+=1
  # If distinct x,x' share one legal upward target y, ||x-x'||_1 <= 4.
  if N==12000 and h==0:
   preds=[]
   # representative y from actual 12-label integer count family, 36 downward predecessors
   y=tuple(n)
   for a in range(0,12,2):
    for b in range(1,12,2):
     x=list(y);x[a]-=1;x[b]+=1
     preds.append(tuple(x))
   assert len(preds)==36 and len(set(preds))==36;nchecks+=1
   assert all(sum(abs(u-v) for u,v in zip(x,xx))<=4 for x,xx in itertools.combinations(preds,2));nchecks+=1
# exact complement identity as a generic algebraic identity for any two equal-total vectors
for v in range(1,37):
 O=[Fraction(i+1,30) for i in range(v)]
 I=list(reversed(O)); T=set(range(0,v,2)); A=set(range(v))-T
 Nbr=set(range(v))-set(range(0,v,3))
 R=set(range(v))-Nbr
 # exact only if use Nbr=neighbourhood in some graph: algebra depends only on totals
 assert sum((I[i] for i in Nbr),Fraction())-sum((O[i] for i in T),Fraction())==sum((O[i] for i in A),Fraction())-sum((I[i] for i in R),Fraction())
 nchecks+=1
# Fraction-only inequality for E upper correction, nondeparture row lower and Hall margin
assert Fraction(18,100)-Fraction(7,1000)==Fraction(173,1000);nchecks+=1
for num in range(1,31):
 E=Fraction(num,10**8)
 A=Fraction(1000,173)*E+Fraction(1,10**8)
 assert Fraction(173,1000)*A>=E
 nchecks+=1
print(json.dumps({'bounded_checks':nchecks,'failures':0,'gplus_order':12,'base_pattern_orbit_min_l1':36,'claim':'checks are finite exact algebraic and integer orbit geometry; NOT full Hall or asymptotic MT theorem'},indent=2))
