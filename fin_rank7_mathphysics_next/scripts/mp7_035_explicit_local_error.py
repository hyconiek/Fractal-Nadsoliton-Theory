#!/usr/bin/env python3
import json, math, pathlib
import numpy as np
from mpmath import iv, mp
mp.dps=80; iv.dps=50
R7P=pathlib.Path('/mnt/data/r7p_source/unpacked/fin_rank7_followup')
c=json.load(open(R7P/'certificates/R7P-026_equal_energy_event.json'))
weights=json.load(open('/mnt/data/fin_rank7_mathphysics_next/results/MP7-034_local_phase_weights.json'))

def rat(s):
 a,b=s.split('/') if '/' in s else (s,'1'); return mp.mpf(a)/mp.mpf(b)
def bnd(x): return float(x.a),float(x.b)
def absup(x): return max(abs(float(x.a)),abs(float(x.b)))
lam={}
for k,v in c['spectral_intervals'].items(): lam[int(k)]=iv.mpf([str(rat(v[0])),str(rat(v[1]))])
X=[]
for j in range(12):
 row=[]
 for k in (3,4,5):
  sc=iv.sqrt(lam[k]/6); ang=iv.pi*iv.mpf(2*k*j)/12
  row.extend([sc*iv.cos(ang),sc*iv.sin(ang)])
 row.append(iv.sqrt(lam[6]/12)*((-1)**j)); X.append(row)
g=iv.mpf(c['root_box'][4]); gmid=(float(c['root_box'][4][0])+float(c['root_box'][4][1]))/2
sbox=c['root_box'][:4]; cent=[(float(a)+float(b))/2 for a,b in sbox]; unc=max((float(b)-float(a))/2 for a,b in sbox)
theta_loc=[cent[0],0,cent[1],0,cent[2],0,cent[3]]
theta_uni=[0.0]*7

def cap_bounds(theta0,r):
 rr=r+(unc if theta0 is theta_loc else 0)
 T=[iv.mpf([str(v-rr),str(v+rr)]) for v in theta0]
 h=[sum(X[j][a]*T[a] for a in range(7)) for j in range(12)]
 e=[iv.exp(z) for z in h]; Z=sum(e,iv.mpf(0)); p=[z/Z for z in e]
 mu=[sum(p[j]*X[j][a] for j in range(12)) for a in range(7)]
 Y=[[X[j][a]-mu[a] for a in range(7)] for j in range(12)]
 H=[[None]*7 for _ in range(7)]
 for a in range(7):
  for b in range(7):
   cov=sum(p[j]*Y[j][a]*Y[j][b] for j in range(12))
   H[a][b]=(1/g if a==b else iv.mpf(0))-cov
 mid=np.zeros((7,7)); rad=np.zeros((7,7))
 for a in range(7):
  for b in range(7):
   lo,hi=bnd(H[a][b]); mid[a,b]=(lo+hi)/2; rad[a,b]=(hi-lo)/2
 mid=(mid+mid.T)/2; rad=np.maximum(rad,rad.T)
 eig=np.linalg.eigvalsh(mid); err=max(rad.sum(axis=1)); mlb=float(eig[0]-err)
 # third derivative tensor of Phi = - third centered moment
 sq=0.0
 for a in range(7):
  for b in range(7):
   for cc in range(7):
    t=-sum(p[j]*Y[j][a]*Y[j][b]*Y[j][cc] for j in range(12))
    sq += absup(t)**2
 B3=math.sqrt(sq)
 return {'radius':r,'lambda_min_lower':mlb,'T3_frobenius_upper':B3,'mid_min_eig':float(eig[0]),'matrix_radius_norm':float(err)}

loc=cap_bounds(theta_loc,0.01)
uni=None
for rr in [0.01,0.005,0.002,0.001,0.0005]:
    z=cap_bounds(theta_uni,rr)
    if z['lambda_min_lower']>0:
        uni=z; break
if uni is None: raise RuntimeError('uniform cap positivity failed')
# central Hessian determinant bounds from MP7-034. For tail constant use safe upper sqrt(detH)=sqrt(detG/g^7)
detG_loc_hi=weights['det_G7_localized'][1]; detG_uni_hi=weights['det_G7_uniform'][1]
glo=float(c['root_box'][4][0])
sqrt_detH_loc_hi=math.sqrt(detG_loc_hi/(glo**7))
sqrt_detH_uni_hi=math.sqrt(detG_uni_hi/(glo**7))

def rel_error_bound(cap,sqrt_detH_hi,N,cscale):
 m=cap['lambda_min_lower']; B=cap['T3_frobenius_upper']; d=7
 rho=cscale*math.sqrt(math.log(N)/N)
 if rho>cap['radius']: return None
 delta=B*cscale**3*(math.log(N)**1.5)/(6*math.sqrt(N))
 # Gaussian H0 tail and actual strong-convexity tail, union bound in eigen/isotropic coordinates
 # Use m for both, conservative.
 qg=min(1.0,2*d*math.exp(-N*m*rho*rho/(2*d)))
 C=sqrt_detH_hi/(m**(d/2))
 qa=min(1e300, C*2*d*math.exp(-N*m*rho*rho/(2*d)))
 lower=max(0.0,math.exp(-delta)*(1-qg))
 upper=math.exp(delta)+qa
 return {'rho':rho,'delta':delta,'gaussian_tail_bound':qg,'actual_tail_over_gaussian_bound':qa,'ratio_lower':lower,'ratio_upper':upper,'relative_error_upper':max(1-lower,upper-1)}

def find_N(cap,sqrt_detH_hi,target):
 m=cap['lambda_min_lower']
 # scan c^2*m/14 in [0.5,4] and powers of ten/binary refine N
 best=None
 for exponent in np.linspace(0.5,4.0,36):
  cscale=math.sqrt(14*exponent/m)
  # find first N on log grid satisfying cap and target
  lo=None
  for logN in np.linspace(4,30,521):
   N=10**logN; z=rel_error_bound(cap,sqrt_detH_hi,N,cscale)
   if z and z['relative_error_upper']<=target:
    lo=N; break
  if lo is None: continue
  # binary in log space between /1.3 and current-ish broad bracket
  a=max(1.0,lo/10**0.05); b=lo
  for _ in range(60):
   mid=math.sqrt(a*b); z=rel_error_bound(cap,sqrt_detH_hi,mid,cscale)
   if z and z['relative_error_upper']<=target:b=mid
   else:a=mid
  z=rel_error_bound(cap,sqrt_detH_hi,b,cscale)
  cand=(b,cscale,exponent,z)
  if best is None or b<best[0]: best=cand
 return best

out={'task':'MP7-035-explicit-local-error','scientific_state':'PROVED_INTERVAL_ASSISTED_EXPLICIT_LOCAL_LAPLACE_ERROR','localized_cap':loc,'uniform_cap':uni,'bounds':{}}
for target in [0.25,0.10,0.05]:
 bl=find_N(loc,sqrt_detH_loc_hi,target); bu=find_N(uni,sqrt_detH_uni_hi,target)
 # both caps simultaneously: use max N, evaluate each at its own optimized c
 N=max(bl[0],bu[0]); zl=rel_error_bound(loc,sqrt_detH_loc_hi,N,bl[1]); zu=rel_error_bound(uni,sqrt_detH_uni_hi,N,bu[1])
 out['bounds'][str(target)]={'N0_both_caps':N,'localized_c':bl[1],'uniform_c':bu[1],'localized_at_N0':zl,'uniform_at_N0':zu}
out['scope']='Explicit relative error for each fixed cap at g_eq only. It does not yet include an explicit global complement-mass bound or uniform g_eq+c/N continuation error.'
print(json.dumps(out,indent=2))
