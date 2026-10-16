#!/usr/bin/env python3
"""Bounded, exact-integer checks for the new finite Fourier/quadrant lemma.
Not a proof of the full q12/D_N all-N signed-cut bound.
"""
from fractions import Fraction
from collections import Counter
import json
from pathlib import Path
out=Path(__file__).resolve().parent
j=range(12)
C=[1,0,-1,0]*3
S=[0,1,0,-1]*3
sig=[1 if k%2==0 else -1 for k in j]
up=[(a,b) for a in j if sig[a]==-1 for b in j if sig[b]==1]
down=[(a,b) for a in j if sig[a]==1 for b in j if sig[b]==-1]
checks=[]
def record(name,cond,detail=''):
 checks.append({'test':name,'pass':bool(cond),'detail':str(detail)})
record('36 upward pairs',len(up)==36)
record('36 downward pairs',len(down)==36)
for direction,pairs in [('up',up),('down',down)]:
 buckets=Counter((C[b]-C[a],S[b]-S[a]) for a,b in pairs)
 record(direction+' exactly nine per diagonal',buckets==Counter({(i,k):9 for i in (-1,1) for k in (-1,1)}),buckets)
 for a,b in pairs:
  dC,dS=C[b]-C[a],S[b]-S[a]
  record(f'{direction} ({a},{b}) exact parity/fourier',dC*dC+dS*dS==2 and sig[b]-sig[a]==(2 if direction=='up' else -2))
# Exact R step identities with integer N, C,S representative values.
for N in (120,1200,120000):
 for c,s in [(-20,8),(0,0),(9,-14),(300,250),(0,-11)]:
  for dC,dS in ((-1,-1),(-1,1),(1,-1),(1,1)):
   r_old=Fraction(2*(c*c+s*s),N)
   r_new=Fraction(2*((c+dC)**2+(s+dS)**2),N)
   expected=Fraction(4*(c*dC+s*dS+1),N)
   record(f'R finite difference {N}:{c}:{s}:{dC}:{dS}',r_new-r_old==expected)
# Check exact continuous Lipschitz bound for T(z) independently in one harmless numeric calibration;
# not the actual exact A7 model or a real D_N classification.
import math
for gamma in (0.2,0.5,0.9):
 L=4/(3*gamma*math.sqrt(2*math.e)); tmin=2*math.log(2)/gamma; tmax=2*math.log(8/3)/gamma
 for N in (1200000,120000000000):
  eta=4*math.sqrt(tmax)/math.sqrt(N)+4/N+2*L*N**(-.75)
  if eta<tmin/2 and (N*tmin/4)**.5 >= 1+L/2*N**.25:
   for c,s in [(-math.ceil((N*tmin/4)**.5),0),(math.ceil((N*tmin/4)**.5),10)]:
    m=abs(c)+abs(s)
    inward=4*(1-m)/N
    outward=4*(1+m)/N
    record(f'conditional shell safe radial step gamma={gamma}, N={N}, C={c},S={s}',inward<=-2*L*N**(-.75) and outward>=2*L*N**(-.75))
result={'scope':'Exact finite Fourier and rational radial checks, plus numeric checks of explicitly stated sufficient threshold; not actual D_N or all-N cut proof',
 'total_checks':len(checks),'passed':sum(x['pass'] for x in checks),'failed':[x for x in checks if not x['pass']],
 'quadrants':{'up':dict((str(k),v) for k,v in Counter((C[b]-C[a],S[b]-S[a]) for a,b in up).items()),'down':dict((str(k),v) for k,v in Counter((C[b]-C[a],S[b]-S[a]) for a,b in down).items())}}
(out/'verification.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({'total':result['total_checks'],'passed':result['passed'],'failed':len(result['failed'])}))
if result['failed']: raise SystemExit(1)
