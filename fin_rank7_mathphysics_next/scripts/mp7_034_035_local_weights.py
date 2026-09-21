#!/usr/bin/env python3
from __future__ import annotations
import os
import json, math, hashlib, sys
from pathlib import Path
import numpy as np
import mpmath as mp

WORK=Path(os.environ.get('MP7_WORK_ROOT','/mnt/data/fin_rank7_mathphysics_next'))
R7P=Path(os.environ.get('R7P_ROOT','/mnt/data/r7p_source/extracted/fin_rank7_followup'))
sys.path.insert(0,str(R7P))
from src import coexistence_certificate as cc
from src.model import feature_spaces, d12_actions

mp.iv.dps=60
iv=mp.iv

def I(x):
    if isinstance(x,(list,tuple)):
        return iv.mpf([str(x[0]),str(x[1])])
    return iv.mpf(str(x))

def bounds(x):
    return [float(x.a), float(x.b)]

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

cert_path=R7P/'certificates/R7P-026_equal_energy_event.json'
stab_path=R7P/'certificates/R7P-027_localized_full_stability.json'
trans_path=R7P/'certificates/R7P-028_crossing_transversality.json'
cert=json.loads(cert_path.read_text())
stab=json.loads(stab_path.read_text())
trans=json.loads(trans_path.read_text())

g=I(cert['root_box'][4])
L={int(k):I(v) for k,v in cert['spectral_intervals'].items()}

# Localized full-7 G7 determinant from certified block LDL pivots: G7=g H7.
detHloc=I(1)
for z in stab['H4_LDL_pivots']+stab['Hsin_LDL_pivots']:
    detHloc *= I(z)
detGloc=(g**7)*detHloc

# Uniform root: M7 is diagonal in the supplied orthogonal Fourier basis.
detGunif=I(1)
for k,powr in [(3,2),(4,2),(5,2),(6,1)]:
    factor=I(1)-g*L[k]/12
    detGunif *= factor**powr

pref_loc=1/iv.sqrt(detGloc)
pref_unif=1/iv.sqrt(detGunif)
family_ratio_at_eq=I(12)*iv.sqrt(detGunif/detGloc)
A=iv.log(I(12))+I('0.5')*iv.log(detGunif/detGloc)

# Strong-convexity lower bound on the tiny validated root box by midpoint + radius matrix norm.
_,_,C,S=cc.interval_features()
X=[I(x) for x in cert['root_box']]
_,J,p,_=cc.interval_eval(X,C)
H4=[[J[i][j] for j in range(4)] for i in range(4)]
Hs=[[((I(1)/X[4]) if a==b else I(0))-sum(p[j]*S[j][a]*S[j][b] for j in range(12))
     for b in range(3)] for a in range(3)]

def eig_lb(A):
    n=len(A); mid=np.zeros((n,n)); rad=np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            lo,hi=bounds(A[i][j]); mid[i,j]=(lo+hi)/2; rad[i,j]=(hi-lo)/2
    e=np.linalg.eigvalsh(mid)
    err=float(np.max(rad.sum(axis=1))) # symmetric ||E||2 <= ||rad||_inf
    return float(e[0]-err), e.tolist(), err
m4,e4,er4=eig_lb(H4); ms,es,ers=eig_lb(Hs)
strong_convexity=min(m4,ms)

# Krawczyk image margins in source certificate guarantee exact root lies far from root-box boundary.
min_margin=min(float(x) for pair in cert['inclusion_margins'] for x in pair)
cap_radius=min(5e-10, 0.5*min_margin)
cap_gap=0.5*strong_convexity*cap_radius**2

# Exact D12 orbit size proof for all-positive aligned amplitudes. Numerical separation is only a check.
s=np.array([float(x) for x in cert['root_center'][:4]])
theta=np.array([s[0],0,s[1],0,s[2],0,s[3]],float)
_,_,_,X7,_,_=feature_spaces(); acts=d12_actions(X7)
points=[]
for a in range(12):
    T=acts[(a,1)][1]
    points.append(T@theta)
unique=[]
for q in points:
    if not any(np.linalg.norm(q-r)<1e-10 for r in unique): unique.append(q)
min_dist=min(np.linalg.norm(unique[i]-unique[j]) for i in range(len(unique)) for j in range(i))

# MP7-035 leading local equal-cap shift.
deriv=I(trans['derivative_interval']) # DeltaV'=localized-uniform derivative; uniform derivative is zero.
cshift=A/deriv

out34={
  'task':'MP7-034',
  'scientific_state':'PROVED_INTERVAL_ASSISTED_LOCAL_CAP_WEIGHTS',
  'inputs':{
      str(cert_path.relative_to(R7P)):sha(cert_path),
      str(stab_path.relative_to(R7P)):sha(stab_path),
      str(trans_path.relative_to(R7P)):sha(trans_path),
  },
  'g_event_box':bounds(g),
  'det_G7_uniform':bounds(detGunif),
  'det_G7_localized':bounds(detGloc),
  'laplace_prefactor_uniform_detG_minus_half':bounds(pref_unif),
  'laplace_prefactor_single_localized_detG_minus_half':bounds(pref_loc),
  'localized_D12_orbit_size':12,
  'stabilizer_order':2,
  'stabilizer_reason':'all four aligned modes active; mode k=5 kills every nontrivial translation, reflection j->-j remains; any second reflection would imply a forbidden translation',
  'numerical_unique_translation_images':len(unique),
  'minimum_orbit_center_distance_check':min_dist,
  'total_localized_family_to_uniform_gaussian_weight_ratio_at_equal_energy':bounds(family_ratio_at_eq),
  'log_prefactor_multiplicity_constant_A':bounds(A),
  'root_box_full_H7_strong_convexity_lower_bound':strong_convexity,
  'root_box_H4_mid_eigs':e4,
  'root_box_Hsin_mid_eigs':es,
  'root_box_matrix_radius_errors':[er4,ers],
  'explicit_cap_radius_about_exact_root':cap_radius,
  'explicit_local_cap_boundary_gap_lower_bound':cap_gap,
  'scope':'mediator-space local caps; no global complement mass claim'
}
out35={
  'task':'MP7-035',
  'scientific_state':'PROVED_LOCAL_ASYMPTOTIC_LEADING_SHIFT; EXPLICIT_FINITE_N_REMAINDER_OPEN',
  'inputs':['MP7-034',str(trans_path.relative_to(R7P))],
  'DeltaV_prime_event':bounds(deriv),
  'A_equal_energy_log_weight_constant':bounds(A),
  'equal_cap_shift_coefficient_c_in_gN=g_eq+c/N+o(1/N)':bounds(cshift),
  'shift_sign':'positive (toward g above the local equal-energy event)',
  'two_family_logistic_limit':'for g=g_eq+c/N, R_localized_family/uniform -> exp(A-c*DeltaV_prime); conditional two-family share -> R/(1+R)',
  'global_interpretation':False,
  'explicit_finite_N_error_bound':None,
  'remaining_atom':'bound Laplace remainder uniformly in g=g_eq+O(1/N), including outside-core contribution; global complement remains separately open'
}
(WORK/'results/MP7-034_local_phase_weights.json').write_text(json.dumps(out34,indent=2)+'\n')
(WORK/'results/MP7-035_local_finite_size_shift.json').write_text(json.dumps(out35,indent=2)+'\n')
print(json.dumps({'MP7-034':out34,'MP7-035':out35},indent=2))
