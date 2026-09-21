"""R7P-090 validated local isolation of quartic phase roots.

Uses high-precision interval arithmetic (mpmath.iv) for a Krawczyk inclusion
on lifted R^3 phase boxes. The periodic equations make seam-crossing boxes valid.
"""
from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np
import mpmath as mp
from .phase_cumulants import _phase_term_data,k4_phase_value_grad_hess

FIXTURE=(0.1131879146,0.1698528641,0.2269339093,-0.3380663037)

def _terms():
    m2,base,cubic,quartic=_phase_term_data(*FIXTURE)
    return [(float(coef/6),np.asarray(n,float)) for coef,n in cubic]+[(float(coef/24),np.asarray(n,float)) for coef,n in quartic]

def _bounds(x):
    return float(x.a),float(x.b)

def _FH(vars,terms):
    z=mp.iv.mpf('0'); g=[z for _ in range(3)]; H=[[z for _ in range(3)] for __ in range(3)]
    for coef,n in terms:
        ang=sum(mp.iv.mpf(repr(float(n[i])))*vars[i] for i in range(3))
        sn=mp.iv.sin(ang); cs=mp.iv.cos(ang); cc=mp.iv.mpf(repr(float(coef)))
        for i in range(3):
            ni=mp.iv.mpf(repr(float(n[i]))); g[i] += -cc*ni*sn
            for j in range(3): H[i][j] += -cc*ni*mp.iv.mpf(repr(float(n[j])))*cs
    return g,H

def certify_record(rec,radius=1e-7,dps=50):
    mp.iv.dps=dps; terms=_terms(); c=np.asarray(rec['phase'],float)
    X=[mp.iv.mpf([repr(float(v-radius)),repr(float(v+radius))]) for v in c]
    C=[mp.iv.mpf(repr(float(v))) for v in c]
    Fc,Jc=_FH(C,terms); _,JX=_FH(X,terms)
    H0=k4_phase_value_grad_hess(*FIXTURE,c)[2]
    A=np.linalg.inv(H0)
    zero=mp.iv.mpf('0'); one=mp.iv.mpf('1'); D=mp.iv.mpf([repr(-radius),repr(radius)])
    b=[]
    for i in range(3):
        ss=zero
        for j in range(3): ss += mp.iv.mpf(repr(float(A[i,j])))*Fc[j]
        b.append(-ss)
    B=[]
    for i in range(3):
        row=[]
        for j in range(3):
            ss=zero
            for k in range(3): ss += mp.iv.mpf(repr(float(A[i,k])))*JX[k][j]
            row.append((one if i==j else zero)-ss)
        B.append(row)
    K=[]; maxoff=0.0
    for i in range(3):
        ss=b[i]
        for j in range(3): ss += B[i][j]*D
        lo,hi=_bounds(ss); K.append([lo,hi]); maxoff=max(maxoff,abs(lo),abs(hi))
    inclusion=maxoff < radius
    # Hessian inertia by interval perturbation from point Hessian.
    E=np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            lo,hi=_bounds(JX[i][j]); E[i,j]=max(abs(lo-H0[i,j]),abs(hi-H0[i,j]))
    err=float(np.linalg.norm(E,'fro')); eig=np.linalg.eigvalsh(H0); margin=float(np.min(np.abs(eig))-err)
    return {'id':rec['id'],'radius':radius,'krawczyk_offset':K,'max_abs_krawczyk_offset':maxoff,
            'inclusion':bool(inclusion),'hessian_mid_eigenvalues':eig.tolist(),
            'hessian_interval_spectral_error_bound':err,'inertia_margin':margin,
            'negative_index':int(np.sum(eig<0)) if margin>0 else None}

def certify_catalog(catalog_path,out_path,radius=1e-7):
    cat=json.load(open(catalog_path)); cert=[certify_record(r,radius) for r in cat['roots']]
    centers=[np.asarray(r['phase'],float) for r in cat['roots']]
    mind=10.0
    for i in range(len(centers)):
        for j in range(i):
            d=np.abs(centers[i]-centers[j]); d=np.minimum(d,2*math.pi-d)
            mind=min(mind,float(np.max(d))) # sufficient disjointness metric: max coordinate torus separation
    # For every pair, disjoint if some coordinate separation > 2r; record minimum of those max separations.
    out={'method':'mpmath.iv Krawczyk on lifted periodic phase boxes','radius':radius,'count':len(cert),
         'all_inclusions':all(x['inclusion'] for x in cert),'all_inertia_certified':all(x['inertia_margin']>0 for x in cert),
         'minimum_pairwise_max_coordinate_torus_separation':mind,
         'all_boxes_pairwise_disjoint':mind>2*radius,
         'min_inertia_margin':min(x['inertia_margin'] for x in cert),
         'max_krawczyk_offset':max(x['max_abs_krawczyk_offset'] for x in cert),'roots':cert}
    Path(out_path).write_text(json.dumps(out,indent=2)+'\n'); return out

if __name__=='__main__':
    import sys
    out=certify_catalog(sys.argv[1],sys.argv[2])
    print(json.dumps({k:v for k,v in out.items() if k!='roots'},indent=2))
