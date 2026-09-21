#!/usr/bin/env python3
import json, math, pathlib, os
from mpmath import iv, mp
mp.dps=60; iv.dps=40
WORK=pathlib.Path(os.environ.get('MP7_WORK_ROOT',pathlib.Path(__file__).resolve().parents[1]))
R7P=pathlib.Path(os.environ.get('R7P_ROOT','/mnt/data/r7p_source/unpacked/fin_rank7_followup'))
c=json.load(open(R7P/'certificates/R7P-026_equal_energy_event.json'))
w=json.load(open(WORK/'results/MP7-034_local_phase_weights.json'))

def rat(s):
 a,b=s.split('/') if '/' in s else (s,'1'); return mp.mpf(a)/mp.mpf(b)
def absup(x): return max(abs(float(x.a)),abs(float(x.b)))
lam={int(k):iv.mpf([str(rat(v[0])),str(rat(v[1]))]) for k,v in c['spectral_intervals'].items()}
X=[]
for j in range(12):
 row=[]
 for k in (3,4,5):
  sc=iv.sqrt(lam[k]/6); ang=iv.pi*iv.mpf(2*k*j)/12
  row.extend([sc*iv.cos(ang),sc*iv.sin(ang)])
 row.append(iv.sqrt(lam[6]/12)*((-1)**j)); X.append(row)
cent=[(float(a)+float(b))/2 for a,b in c['root_box'][:4]]
theta_loc=[cent[0],0,cent[1],0,cent[2],0,cent[3]]; theta_uni=[0.0]*7
unc=max((float(b)-float(a))/2 for a,b in c['root_box'][:4])
glo=float(c['root_box'][4][0]); ghi=float(c['root_box'][4][1])
# localized exact-root central lower from MP7-034 tiny validated box
m0_loc=w['root_box_full_H7_strong_convexity_lower_bound']
# uniform exact H diag; worst retained lambda/g endpoints
m0_uni=min(1/ghi-float(rat(c['spectral_intervals'][k][1]))/12 for k in ['3','4','5','6'])

def B3(theta0,r,add_unc=False):
 rr=r+(unc if add_unc else 0)
 T=[iv.mpf([str(v-rr),str(v+rr)]) for v in theta0]
 h=[sum(X[j][a]*T[a] for a in range(7)) for j in range(12)]
 e=[iv.exp(z) for z in h]; Z=sum(e,iv.mpf(0)); p=[z/Z for z in e]
 mu=[sum(p[j]*X[j][a] for j in range(12)) for a in range(7)]
 Y=[[X[j][a]-mu[a] for a in range(7)] for j in range(12)]
 sq=0.0
 for a in range(7):
  for b in range(7):
   for cc in range(7):
    t=-sum(p[j]*Y[j][a]*Y[j][b]*Y[j][cc] for j in range(12)); sq+=absup(t)**2
 return math.sqrt(sq)
out={'m0_localized':m0_loc,'m0_uniform':m0_uni,'localized':[],'uniform':[]}
for r in [0.005,0.01,0.015,0.02,0.03,0.04,0.05,0.075,0.1]:
 for name,t,m0,u in [('localized',theta_loc,m0_loc,True),('uniform',theta_uni,m0_uni,False)]:
  B=B3(t,r,u); ml=m0-B*r
  out[name].append({'radius':r,'T3_frob_upper':B,'lipschitz_lambda_min_lower':ml})
print(json.dumps(out,indent=2))
