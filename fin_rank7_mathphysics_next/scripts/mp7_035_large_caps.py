#!/usr/bin/env python3
import json, math, pathlib
import numpy as np
from mpmath import iv, mp
mp.dps=80; iv.dps=50
R7P=pathlib.Path('/mnt/data/r7p_source/unpacked/fin_rank7_followup')
c=json.load(open(R7P/'certificates/R7P-026_equal_energy_event.json'))

def rat(s):
 a,b=s.split('/') if '/' in s else (s,'1'); return mp.mpf(a)/mp.mpf(b)
def bnd(x): return float(x.a),float(x.b)
lam={}
for k,v in c['spectral_intervals'].items(): lam[int(k)]=iv.mpf([str(rat(v[0])),str(rat(v[1]))])
# interval X7 exact trig constants via iv trig
X=[]
for j in range(12):
 row=[]
 for k in (3,4,5):
  sc=iv.sqrt(lam[k]/6); ang=iv.pi*iv.mpf(2*k*j)/12
  row.extend([sc*iv.cos(ang),sc*iv.sin(ang)])
 row.append(iv.sqrt(lam[6]/12)*((-1)**j)); X.append(row)
# exact root center intervals -> choose midpoint floats as cap center and include source root uncertainty
sbox=c['root_box'][:4]
cent=[(float(a)+float(b))/2 for a,b in sbox]
theta0=[cent[0],0,cent[1],0,cent[2],0,cent[3]]
g=iv.mpf(c['root_box'][4])
unc=max((float(b)-float(a))/2 for a,b in sbox)

def calc(r):
 # coordinate halfwidth includes root uncertainty
 rr=r+unc
 T=[iv.mpf([str(v-rr),str(v+rr)]) for v in theta0]
 h=[sum(X[j][a]*T[a] for a in range(7)) for j in range(12)]
 e=[iv.exp(z) for z in h]; Z=sum(e,iv.mpf(0)); p=[z/Z for z in e]
 mu=[sum(p[j]*X[j][a] for j in range(12)) for a in range(7)]
 H=[[None]*7 for _ in range(7)]
 for a in range(7):
  for b in range(7):
   cov=sum(p[j]*(X[j][a]-mu[a])*(X[j][b]-mu[b]) for j in range(12))
   H[a][b]=(1/g if a==b else iv.mpf(0))-cov
 mid=np.zeros((7,7)); rad=np.zeros((7,7))
 for a in range(7):
  for b in range(7):
   lo,hi=bnd(H[a][b]);mid[a,b]=(lo+hi)/2;rad[a,b]=(hi-lo)/2
 mid=(mid+mid.T)/2; rad=np.maximum(rad,rad.T)
 eig=np.linalg.eigvalsh(mid); err=max(rad.sum(axis=1)); lb=eig[0]-err
 return {'coordinate_halfwidth':r,'euclidean_ball_radius_certified':r,'mid_min_eig':float(eig[0]),'radius_norm_bound':float(err),'lambda_min_lower':float(lb)}

rs=[0.001,0.002,0.003,0.004,0.005,0.0075,0.01,0.015,0.02,0.03]
out=[calc(r) for r in rs]
print(json.dumps(out,indent=2))
