"""MP7-041: exact symmetry-locked amplitude-to-phase continuation near full-phase root 14."""
from __future__ import annotations
import os
import json,math,sys
from pathlib import Path
import numpy as np, mpmath as mp
ROOT=Path(os.environ.get('MP7_WORK_ROOT',Path(__file__).resolve().parents[1]))
R7P=Path(os.environ.get('R7P_ROOT','/mnt/data/r7p_source/unpacked/fin_rank7_followup'))
sys.path.insert(0,str(R7P))
from src.model import feature_spaces
from src.coexistence_certificate import _I,_bound,interval_ldlt
from src.phase_cumulants import full_phase_value_grad_hess
mp.iv.dps=60; iv=mp.iv
N=12
FIX=np.array([0.1131879146,0.1698528641,0.2269339093,-0.3380663037],float)
PH=np.array([math.pi/2,2*math.pi/3,5*math.pi/6],float) # odd translation a=1 of aligned +z6 field
REL=0.001

def phase_hessian_iv(A,P):
    rs=A[:3]; z6=A[3]; ks=[3,4,5]
    h=[]; hp=[]; hpp=[]
    for jj in range(N):
        hv=z6*_I((-1)**jj)/iv.sqrt(12); row=[]; row2=[]
        for r,k,p in zip(rs,ks,P):
            ang=_I(2*math.pi*k*jj/N)+p
            hv += r*iv.cos(ang)/iv.sqrt(3)
            row.append(-r*iv.sin(ang)/iv.sqrt(3)); row2.append(-r*iv.cos(ang)/iv.sqrt(3))
        h.append(hv);hp.append(row);hpp.append(row2)
    ew=[iv.exp(x) for x in h];Z=sum(ew,_I(0));prob=[x/Z for x in ew]
    grad=[sum(prob[j]*hp[j][i] for j in range(N)) for i in range(3)]
    H=[[_I(0) for _ in range(3)] for __ in range(3)]
    for a in range(3):
        for b in range(3):
            H[a][b]=sum(prob[j]*hp[j][a]*hp[j][b] for j in range(N))-grad[a]*grad[b]
            if a==b: H[a][b]+=sum(prob[j]*hpp[j][a] for j in range(N))
    return grad,H

def center_blocks():
    j=np.arange(N,dtype=float);ks=np.array([3.,4.,5.]); ang=2*np.pi*j[:,None]*ks[None,:]/N+PH[None,:]
    q=np.column_stack([np.cos(ang)/math.sqrt(3), ((-1.0)**j)/math.sqrt(12)])
    h=q@FIX; m=h.max(); w=np.exp(h-m); p=w/w.sum(); mu=p@q; Y=q-mu
    Kaa=Y.T@(p[:,None]*Y)
    _,_,L,_,_,_=feature_spaces(); Q=np.diag([2/L[3],2/L[4],2/L[5],1/L[6]])
    Kpp=full_phase_value_grad_hess(*FIX,PH)[2]
    return mu,Kaa,Q,Kpp,L

def main():
    A=[_I([x-abs(x)*REL,x+abs(x)*REL]) for x in FIX]
    P=[_I(repr(float(x))) for x in PH]
    _,Kpp_iv=phase_hessian_iv(A,P)
    # Phi phase Hessian = -Kpp. Prove it positive definite on box.
    D,Lfac,ok=interval_ldlt([[-Kpp_iv[i][j] for j in range(3)] for i in range(3)])
    piv=[_bound(x) if x is not None else None for x in D]
    if not ok or any(x is None or float(x[0])<=0 for x in piv): raise RuntimeError('phase block failed')
    mu,Kaa,Q,Kpp,Lspec=center_blocks()
    # Since this phase lock is an exact translated aligned field for every amplitude,
    # grad_phi logZ=0 identically as a function of amplitudes. Hence K_a_phi=0.
    ratios=(Q@FIX)/mu
    out={
      'task':'MP7-041',
      'scientific_state':'PROVED_INTERVAL_ASSISTED_EXACT_SYMMETRY_CONTINUATION',
      'selected_full_phase_root_id':14,
      'fixture_amplitudes':FIX.tolist(),
      'phase_lock':PH.tolist(),
      'symmetry_reason':'odd label translation j->j+1 maps this negative-z6 phase lock to the aligned positive-z6 field; phase gradient therefore vanishes identically for all amplitudes preserving signs',
      'amplitude_box':{
        'relative_radius':REL,
        'bounds':[[float(x.a),float(x.b)] for x in A],
        'signs_preserved':True
      },
      'phase_map':'phi(a) is exactly constant at the displayed phase lock on this box',
      'phase_phi_hessian_positive_LDL_pivots':piv,
      'ift_nonsingular':True,
      'cross_block_K_a_phi':'exactly zero because grad_phi logZ is identically zero along the symmetry-locked amplitude family',
      'center_logZ_amplitude_hessian':Kaa.tolist(),
      'quadratic_metric_Q':Q.tolist(),
      'effective_amplitude_hessian_formula':'H_eff(g)=Q/g-Kaa at the center; more generally Q/g-Hess_aa logZ along the exact phase-locked branch',
      'center_effective_H_eigs':{str(g):np.linalg.eigvalsh(Q/g-Kaa).tolist() for g in [3.5,3.7,4.0,5.0]},
      'fixed_fixture_radial_gain_candidates':ratios.tolist(),
      'radial_gain_spread':float(ratios.max()-ratios.min()),
      'mp7_040_consistency':'The near-equality of radial gain candidates is diagnostic only; MP7-040 interval arithmetic proves they are not exactly equal at this frozen fixture.',
      'consequence':'Nearby full equilibria on this symmetry family can be searched as a four-amplitude radial problem without re-solving phases.',
      'nonconclusions':['does not prove a nearby full equilibrium exists','does not extend the fixed-fixture 60-root census globally in amplitude space','does not imply global minimality']
    }
    (ROOT/'results/MP7-041_phase_continuation.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
