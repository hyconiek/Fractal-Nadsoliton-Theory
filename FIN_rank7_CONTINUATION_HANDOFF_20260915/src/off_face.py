"""R7P-065--067: stable off-face target, adversarial search, exact compactification.

Exact statements are kept separate from numerical search.  The compactification
uses x_k=exp(-J_k) and preserves the irrational k=5 exponents exactly.
"""
from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
import json, math
from functools import lru_cache
import numpy as np
import sympy as sp
from scipy.optimize import differential_evolution, minimize_scalar

from model import feature_spaces
from boundary_ising import strict_intervals
from intervals import sqrt_interval

ROOT=Path(__file__).resolve().parents[1]
KS=(3,4,5,6)


@lru_cache(maxsize=1)
def constants():
    _,_,L,_,C4,_=feature_spaces()
    sigma=(2*L[3]*(L[4]+L[5])-L[4]*L[5])/(24*L[3])
    scales=np.array([math.sqrt(L[3]/6),math.sqrt(L[4]/6),math.sqrt(L[5]/6),math.sqrt(L[6]/12)])
    j=np.arange(12,dtype=float)
    obs=np.column_stack([np.cos(2*np.pi*k*j/12) for k in (3,4,5)]+[(-1.)**j])
    return L,sigma,scales,obs,C4


def strict_target_data():
    L=strict_intervals(); l3,l4,l5,l6=L[3],L[4],L[5],L[6]
    sigma=(2*l3*(l4+l5)-l4*l5)/(24*l3)
    eta_lower=1-l6/(12*sigma)
    assert sigma.lo>0 and eta_lower.lo>0
    return {
      'sigma_interval':[str(sigma.lo),str(sigma.hi)],
      'eta_global_lower_interval':[str(eta_lower.lo),str(eta_lower.hi)],
      'eta_identity':'eta=1-lambda6*q*(1-q)/(3*sigma)',
      'eta_bound':'q(1-q)<=1/4 => eta>=1-lambda6/(12*sigma)>0',
    }


def stable_target_spec():
    return {
      'id':'R7P-065-off-face-target',
      'domain':'J3,J4,J5,J6>=0, equivalently s3,s4,s5,s6>=0 under positive spectral scalings',
      'target':'lambda2(M4)<=sigma_* where M4=Cov_p(C4), p_j proportional exp(J3 c3j+J4 c4j+J5 c5j+J6 c6j)',
      'parity_identity':'M4=W_par+b b^T with W_par fourth row/column zero and b6^2=lambda6*q*(1-q)/3',
      'stable_reduction':'Mtilde=W_par[0:3,0:3]+b[0:3]b[0:3]^T/eta',
      'inertia_equivalence':'n_-(sigma I4-M4)=n_-(sigma I3-Mtilde), because the scalar Schur pivot sigma-b6^2=sigma*eta is strictly positive',
      'equivalent_target':'lambda2(Mtilde)<=sigma_* (threshold equality retained); never invert the full sigma I-W_par near singularities',
      'route':'factored Schur/inertia route for proof work; direct M4 used independently for adversarial numerical checks',
      'dependencies':{'G':'R7P-055 global boundary theorem DONE','H':'R7P-063 intraparity theorem DONE','scope':'no full-7D transfer'},
      'strict_intervals':strict_target_data(),
    }


def p_from_fields(J):
    """Finite-field direct implementation in unscaled Fourier observables."""
    _,_,_,obs,_=constants(); J=np.asarray(J,float)
    h=obs@J; h-=np.max(h); w=np.exp(h); return w/w.sum()


def covariance_direct_fields(J):
    _,_,_,_,C4=constants(); p=p_from_fields(J); mu=p@C4; Y=C4-mu
    return Y.T@(p[:,None]*Y),p


def compact_exponents_symbolic():
    out=[]
    for j in range(12):
        out.append([sp.simplify(1-sp.cos(2*sp.pi*k*j/12)) for k in KS])
    return out


def p_from_compact(x):
    """Exact-chart numerical evaluation on x in [0,1]^4.

    The mathematical convention is x^0=1 also at x=0.  Positive powers of 0
    vanish.  This is the continuous compactified extension of finite fields.
    """
    x=np.asarray(x,float)
    if np.any(x<0) or np.any(x>1): raise ValueError('x must be in [0,1]^4')
    # Numerical cosine values are safe for evaluation; proof strings are symbolic.
    _,_,_,obs,_=constants(); exps=1-obs
    w=np.ones(12)
    for k in range(4):
        if x[k]==0:
            w*=np.where(np.abs(exps[:,k])<1e-13,1.0,0.0)
        else:
            w*=x[k]**exps[:,k]
    if not (w[0]>0 and w.sum()>0): raise AssertionError('compact denominator lost anchor state')
    return w/w.sum()


def covariance_direct_compact(x):
    _,_,_,_,C4=constants(); p=p_from_compact(x); mu=p@C4; Y=C4-mu
    return Y.T@(p[:,None]*Y),p


def lambda2(M): return float(np.linalg.eigvalsh(M)[-2])


def q0_from_three_fields(J3,J4,J5):
    _,_,_,obs,_=constants(); h=obs[:,:3]@np.array([J3,J4,J5],float)
    m=np.max(h); z=np.exp(h-m); return float(z[::2].sum()/z.sum())


def compactification_theorem():
    ex=compact_exponents_symbolic()
    assert all(sp.ask(sp.Q.nonnegative(a)) is not False for row in ex for a in row)
    assert all(a==0 for a in ex[0])
    return {
      'id':'R7P-067-exact-compactification',
      'finite_map':'x_k=exp(-J_k), k=3,4,5,6, maps [0,+infinity)^4 bijectively to (0,1]^4',
      'weights':'w_j(x)=prod_k x_k^(1-cos(2*pi*k*j/12)); p_j=w_j/sum_l w_l',
      'exponent_rows':[[str(a) for a in row] for row in ex],
      'anchor':'j=0 has all exponents zero, hence w_0=1 throughout the closed cube and the denominator never vanishes',
      'closure_theorem':'The formula extends continuously to [0,1]^4. Every unbounded field sequence has a compact-x convergent subsequence; conversely every closed-cube point is approached by finite fields using x_k^(n)=max(x_k,1/n). Therefore its image is exactly the probability-law closure of the full nonnegative four-field family.',
      'multiscale':'Simultaneous and arbitrarily different large-field rates are included automatically by arbitrary approaches to cube faces/corners.',
      'irrational_relation':'The k=5 rows retain exponents 1+-sqrt(3)/2; no independent polynomial variable is introduced, so the shared-field irrational relation is preserved.',
      'no_tail_error':'This is an exact compactification, not a finite cutoff plus discarded-probability perturbation bound.',
    }


def _record_candidate(x,tag,seed=None):
    L,sigma,_,_,_=constants(); M,p=covariance_direct_compact(x); eig=np.linalg.eigvalsh(M)
    rec={'tag':tag,'x':[float(v) for v in x],'lambda2':float(eig[-2]),'gap':float(eig[-2]-sigma),
         'eigenvalues':[float(v) for v in eig],'p':[float(v) for v in p]}
    if seed is not None: rec['seed']=int(seed)
    # Independent finite-field replay when all compact coordinates are safely positive.
    if min(x)>1e-12:
        J=-np.log(np.asarray(x,float)); M2,p2=covariance_direct_fields(J)
        rec['J']=[float(v) for v in J]
        rec['independent_covariance_max_abs_diff']=float(np.max(np.abs(M-M2)))
        rec['independent_p_max_abs_diff']=float(np.max(np.abs(p-p2)))
    return rec


def adversarial_search(seeds=(66066,66166,66266),maxiter=700,popsize=20):
    """Numerical-only R7P-066 search over the *entire compactified closure*.

    Differential evolution is not a proof.  It is supplemented by structured
    boundary/slice searches and direct finite-field replays where applicable.
    """
    _,sigma,_,_,_=constants(); records=[]
    for seed in seeds:
        res=differential_evolution(lambda x: -(lambda2(covariance_direct_compact(x)[0])-sigma),
                                   [(0.0,1.0)]*4,seed=seed,maxiter=maxiter,popsize=popsize,
                                   tol=1e-10,polish=True,workers=1,updating='immediate')
        records.append(_record_candidate(res.x,'DE_closed_cube',seed))
    # Structured extreme-face slice x4=x5=1, x6=0, optimize x3.
    opt=minimize_scalar(lambda x3: -(lambda2(covariance_direct_compact([x3,1,1,0])[0])-sigma),
                        bounds=(1e-6,1),method='bounded',options={'xatol':1e-14})
    records.append(_record_candidate([opt.x,1,1,0],'extreme_face_1D'))
    # Finite-J6 slices near the known saturation face, optimize x3 with x4=x5=1.
    for J6 in (2,4,8,12,20):
        x6=math.exp(-J6)
        opt=minimize_scalar(lambda x3: -(lambda2(covariance_direct_compact([x3,1,1,x6])[0])-sigma),
                            bounds=(1e-6,1),method='bounded',options={'xatol':1e-14})
        records.append(_record_candidate([opt.x,1,1,x6],f'face_finite_J6_{J6:g}'))
    # Off-face structured perturbations around the equality locator.
    Jstar=math.atanh(math.sqrt(1-sigma/(L_global()[3]/6)))
    for a in (1e-4,1e-3,1e-2,5e-2,.1,.25):
        for mode in ('J4','J5','both'):
            J=[Jstar,0.,0.,12.]
            if mode in ('J4','both'): J[1]=a
            if mode in ('J5','both'): J[2]=a
            x=np.exp(-np.array(J))
            records.append(_record_candidate(x,f'local_{mode}_{a:g}'))
    best=max(records,key=lambda r:r['gap'])
    return {
      'status':'NUMERICAL_SEARCH_ONLY',
      'search_domain':'entire exact compactified closure x in [0,1]^4, plus structured slices',
      'seeds':list(seeds),'DE_maxiter':maxiter,'DE_popsize':popsize,
      'best':best,'records':records,
      'interpretation':'No positive gap found is not a theorem. Any positive candidate would require R7P-070 certification.',
    }


def L_global(): return constants()[0]



def aligned_compact_chart():
    """Compact chart aligned with the already-certified boundary cube.

    r=e^(-2J3), s=e^(-3J4/2), t=e^(-J5/2), y=e^(-2J6).
    Reflection-related labels have identical cosine features, so seven aggregate
    states are exact for the four-cosine covariance.
    """
    return {
      'coordinates':{'r':'exp(-2J3)','s':'exp(-3J4/2)','t':'exp(-J5/2)','y':'exp(-2J6)'},
      'cube':'0<=r,s,t,y<=1',
      'even_aggregate_weights':['1','2*s*t^3','r*t^4','2*r*s*t'],
      'odd_aggregate_weights':['2*sqrt(r)*s*t^(2+sqrt(3))*y','2*sqrt(r)*t^2*y','2*sqrt(r)*s*t^(2-sqrt(3))*y'],
      'boundary_match':'At y=0 this is exactly the R7P-048/G boundary chart (r,s,t).',
      'irrational_constraint':'The two odd powers t^(2+-sqrt(3)) use the same t; they are not independent variables.',
      'anchor':'the j=0 even state has weight 1',
    }



def p_from_aligned_compact(r,s,t,y):
    """Twelve-label probability law in the boundary-aligned compact chart."""
    vals=(r,s,t,y)
    if any(v<0 or v>1 for v in vals): raise ValueError('r,s,t,y must lie in [0,1]')
    sr=math.sqrt(r)
    rt3=math.sqrt(3.0)
    w=np.array([
      1.0,
      sr*s*(t**(2+rt3))*y,
      r*s*t,
      sr*(t**2)*y,
      s*t**3,
      sr*s*(t**(2-rt3))*y,
      r*t**4,
      sr*s*(t**(2-rt3))*y,
      s*t**3,
      sr*(t**2)*y,
      r*s*t,
      sr*s*(t**(2+rt3))*y,
    ],dtype=float)
    return w/w.sum()

def reoptimized_parity_slope_certificate():
    """Upgrade the R7P-023 reoptimized extreme-face first-order slope.

    epsilon=1-q.  Allow J3=J3_*+d*epsilon while J4=J5=0.  The two
    threshold branches have first-order shifts delta_i+m_i' d.  Their lower
    envelope is maximized when the two affine shifts are equal.
    """
    L=strict_intervals();l3,l4,l5=L[3],L[4],L[5]
    sigma=(2*l3*(l4+l5)-l4*l5)/(24*l3); a=l3/6
    t2=1-sigma/a; t=sqrt_interval(t2,40)
    disc=(l5-l4)**2+4*l4*l5*t2; sd=sqrt_interval(disc,40)
    m3p=-2*a*t*(1-t2)
    m45p=l4*l5*t*(1-t2)/(6*sd)
    delta3=l3*(2*t2-1)/6
    delta45=-l4*l5*t2/(6*sd)
    dopt=(delta45-delta3)/(m3p-m45p)
    slope=(m3p*delta45-m45p*delta3)/(m3p-m45p)
    assert m3p.hi<0 and m45p.lo>0 and delta3.hi<0 and delta45.hi<0 and slope.hi<0
    q=lambda I:[str(I.lo),str(I.hi)]
    return {
      'id':'R7P-023-reoptimized-envelope-upgrade',
      'parameter':'epsilon=1-q; J3=J3_*+d*epsilon, J4=J5=0',
      'branch_derivatives':{
        'm3_prime':q(m3p),'m45_prime':q(m45p),'delta3':q(delta3),'delta45':q(delta45)},
      'optimal_d_interval':q(dopt),
      'envelope_slope_interval':q(slope),
      'display':[float(slope.lo),float(slope.hi)],
      'conclusion':'The reoptimized first-order lambda2 envelope slope is strictly negative; the imported ~-0.1312828584 value is now interval-certified on this extreme-face two-branch reduction.',
      'scope':'first-order extreme-face reoptimization only; not a finite-radius off-face theorem.'
    }


def local_first_order_cone_status():
    """Analytic R7P-068 progress, deliberately short of an explicit radius.

    Uses the global y=0 boundary theorem and R7P-023 negative-definite parity
    mixing on the threshold eigenspace.  It proves no first-order physical
    tangent can create two supercritical directions.  A validated mixed
    remainder is still required for the task's explicit-neighborhood gate.
    """
    slope=reoptimized_parity_slope_certificate()
    return {
      'id':'R7P-068-first-order-cone-partial',
      'status':'PARTIAL_ANALYTIC',
      'root':'(r_*,s,t,y)=(r_*,1,1,0), with r_*=exp(-2J3_*)=(1-tau_*)/(1+tau_*)',
      'boundary_input':'R7P-055 plus R7P-052: y=0 is globally safe and locally strict away from the unique double root.',
      'mixing_input':'R7P-023: on the two-dimensional threshold eigenspace, the epsilon=1-q mixing derivative has two strictly negative M-eigenvalue shifts (equivalently positive-definite K=sigma I-M derivative).',
      'directional_argument':'For any physical boundary tangent H_b, boundary safety implies lambda_min(P H_b P)<=0 for the M-shift cluster. Adding positive parity mixing adds a negative-definite projected M perturbation. Weyl therefore makes the right directional derivative of lambda2 strictly negative whenever the parity-mixing component is first order; pure boundary directions are nonincreasing and locally strict by R7P-052.',
      'reoptimized_extreme_face':slope,
      'qualitative_conclusion':'No first-order physical tangent direction through the equality point produces lambda2>sigma. By smoothness of the aligned chart near r_*>0,t=1 and the strict boundary tangent-cone certificate, this supports a local-safe neighborhood, but this package does not yet export a machine-checked explicit radius.',
      'remaining_atom':'validated second/higher-order mixed remainder bound (especially tangential paths with y=o(norm(boundary displacement))) sufficient to state an explicit dyadic 4D neighborhood.',
      'claim_level':'not promoted to R7P-068 DONE until an explicit radius/remainder is checked.'
    }

def build_results(run_search=True):
    out={'R7P-065':stable_target_spec(),'R7P-067':compactification_theorem()|{'aligned_chart':aligned_compact_chart()},'R7P-068_partial':local_first_order_cone_status()}
    if run_search: out['R7P-066']=adversarial_search()
    return out

if __name__=='__main__':
    out=build_results(True)
    p=ROOT/'results/R7P-065_067_off_face.json'; p.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
