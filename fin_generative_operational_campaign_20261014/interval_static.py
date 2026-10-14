#!/usr/bin/env python3
"""Independent interval construction of 12-bin N2 static laws for supplied exact decimals.
Certified arithmetic scope is conditional on mpmath.iv backend, NOT all-operation formal proof.
"""
import json, math
from decimal import Decimal, ROUND_FLOOR, ROUND_CEILING
from pathlib import Path
import mpmath as mp
iv=mp.iv;iv.dps=65;Q=12
HERE=Path(__file__).resolve().parent
w=[iv.mpf(0)]
for r in range(1,Q):
 d=min(r,Q-r)
 w.append(iv.cos(iv.mpf('0.18575')*d+iv.mpf('0.1625'))/(1+iv.power(iv.mpf(d),iv.mpf('1.8'))))
la=[]
for k in range(Q):
 la.append(sum((w[r]*(1-iv.cos(2*iv.pi*k*r/Q)) for r in range(1,Q)),iv.mpf(0)))
a=[sum((la[k]/6*iv.cos(2*iv.pi*k*d/Q) for k in (3,4,5)),iv.mpf(0))+la[6]/12*((-1)**d) for d in range(Q)]
trace=Q*a[0]
potts=[trace/11*(iv.mpf(1 if d==0 else 0)-iv.mpf(1)/12) for d in range(Q)]
flat=[sum((trace/7)/6*iv.cos(2*iv.pi*k*d/Q) for k in (3,4,5))+trace/7/12*((-1)**d) for d in range(Q)]
rows={'FIN_A7':a,'FULL_POTTS_TRACE':potts,'FLAT_P7_TRACE':flat}
def probs(v):
 z=[iv.exp(iv.mpf('1.5')*x) for x in v];s=sum(z,iv.mpf(0));return [x/s for x in z]
ps={k:probs(v) for k,v in rows.items()}

def bracket(x):
 lo=Decimal(mp.libmp.to_str(x._mpi_[0],90));hi=Decimal(mp.libmp.to_str(x._mpi_[1],90))
 unit=Decimal('0.000000000000001')
 return [str(lo.quantize(unit,rounding=ROUND_FLOOR)),str(hi.quantize(unit,rounding=ROUND_CEILING))]

def nonnegative_abs(x):
 if x.a>0 or x.b<0: return abs(x)
 return iv.mpf([0,max(-x.a,x.b)])
result={'arithmetic_scope':'exact-decimal kernel, exact pi target; interval arithmetic conditional on correctness of mpmath.iv',
        'precision_decimal_digits':iv.dps,
        'N2_g':'3','q':12,
        'N2_FIN_first_mode':bracket(a[0]),'individual_law_enclosures':{k:[bracket(x) for x in p] for k,p in ps.items()},
        'pairwise_tv_public_outward':{}}
ks=list(ps)
for i,x in enumerate(ks):
 for y in ks[i+1:]:
  dd=sum((nonnegative_abs(ps[x][j]-ps[y][j]) for j in range(Q)),iv.mpf(0))/2
  result['pairwise_tv_public_outward'][f'{x}_vs_{y}']=bracket(dd)
p=HERE/'RG-E03/STATIC_INTERVAL_CERTIFICATE.json';p.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps(result['pairwise_tv_public_outward'],indent=2))
