"""R7P-075/078: paid-hypothesis Ginibre cooperativity and monotone iteration.

The external ingredient is Ginibre's Griffiths-inequality framework.  Everything
model-specific is checked here on the finite group G = Z4 x Z3:

* c3,c4,c5,c6 are real parts of characters;
* real parts of characters are positive-definite functions;
* for J_i >= 0, -H = sum_i J_i c_i is a nonnegative combination of them.

Ginibre (CMP 16 (1970), 310--328, DOI 10.1007/BF01646537), Example 4
and the general theorem then give nonnegative means and pair covariances for
real positive-definite observables under Haar/counting measure.

The monotone-iteration theorem below is a mathematical fixed-point algorithm,
not a claim about physical time evolution.
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

# Frequencies for chi_(k4,k3)(a,b)=exp(2*pi*i*(k4*a/4+k3*b/3)).
# beta=-2*pi*b/3, so c4=cos(beta) has frequency (0,-1) modulo 3.
FREQUENCIES = {
    'c3': (1, 0),
    'c4': (0, 2),   # -1 mod 3
    'c5': (1, 1),
    'c6': (2, 0),
}


def group_points():
    return [(a,b) for a in range(4) for b in range(3)]


def character(freq, g):
    k4,k3=freq; a,b=g
    return np.exp(2j*np.pi*(k4*a/4+k3*b/3))


def cosine_character(freq, g):
    return float(np.real(character(freq,g)))


def pd_kernel(freq):
    """Kernel K[g,h]=Re chi(g-h); PSD because K=(v v* + vbar vbar*)/2."""
    pts=group_points(); v=np.array([character(freq,g) for g in pts],complex)
    K=(np.outer(v,np.conjugate(v))+np.outer(np.conjugate(v),v))/2
    return np.real_if_close(K).astype(float)


def finite_group_checks(tol=1e-11):
    pts=group_points()
    rows={}
    for name,freq in FREQUENCIES.items():
        K=pd_kernel(freq)
        eig=np.linalg.eigvalsh(K)
        vals=np.array([cosine_character(freq,g) for g in pts])
        # direct difference-kernel comparison
        K2=np.empty((12,12))
        for i,(a,b) in enumerate(pts):
            for j,(c,d) in enumerate(pts):
                K2[i,j]=cosine_character(freq,((a-c)%4,(b-d)%3))
        rows[name]={
            'frequency':list(freq),
            'kernel_min_eigenvalue':float(eig[0]),
            'kernel_rank':int(np.sum(eig>tol)),
            'uniform_mean':float(vals.mean()),
            'kernel_identity_residual':float(np.max(np.abs(K-K2))),
        }
    return rows


def theorem_record():
    checks=finite_group_checks()
    assert all(v['kernel_min_eigenvalue']>-1e-10 for v in checks.values())
    assert all(abs(v['uniform_mean'])<1e-12 for v in checks.values())
    return {
      'id':'R7P-075-ginibre-cooperativity',
      'status':'EXACT_PROVED_WITH_APPLICABLE_EXTERNAL_THEOREM',
      'group':'Z4 x Z3 with uniform/Haar counting measure',
      'character_frequencies':{k:list(v) for k,v in FREQUENCIES.items()},
      'finite_group_checks':checks,
      'model_specific_exact_argument':[
        'For any character chi, K_chi(g,h)=chi(g-h)=v(g) conjugate(v(h)) is PSD.',
        'For c=Re chi=(chi+conjugate(chi))/2, K_c is the average of two PSD rank-one character kernels, hence PSD.',
        'Nonnegative linear combinations of real positive-definite functions are positive-definite.',
        'With H=-sum_i J_i c_i and J_i>=0, -H is real positive-definite.',
      ],
      'external_lemma':{
        'author':'J. Ginibre',
        'title':"General formulation of Griffiths' inequalities",
        'journal':'Communications in Mathematical Physics 16 (1970), 310-328',
        'doi':'10.1007/BF01646537',
        'use':'Ginibre general theorem plus Example 4: on a compact abelian group with Haar measure, the cone of real positive-definite functions satisfies Q3; if -H is in that cone, expectations and pair covariances of cone observables are nonnegative.',
      },
      'conclusions':{
        'means':'E[c_i] >= 0 for i in {3,4,5,6} and all J_i>=0',
        'covariances':'Cov(c_i,c_j) >= 0 for all i,j and all J_i>=0',
        'scaled_features':'Positive Fourier normalizations preserve the signs.',
        'fixed_point_jacobian':'For T(s)=g E_s[C4] on s>=0, DT=g Cov_s(C4) is entrywise nonnegative.',
        'isotonicity':'T is coordinatewise isotone on the nonnegative amplitude orthant.',
      },
      'nonconclusions':[
        'Strict positivity is not asserted at every boundary point.',
        'This theorem does not apply to the relaxed negative controls with independently split even/odd J5 fields.',
        'Isotonicity of this mathematical fixed-point map is not a physical time-evolution law.',
        'It does not imply the four-amplitude curvature ceiling or any full-seven-coordinate index bound.',
      ],
    }


def monotone_iteration_record():
    """Order-theoretic theorem for T(s)=g E_s[C4] at fixed g>0."""
    # Only the proof skeleton and exact bounds are needed.  For C4_i = scale_i*c_i,
    # |c_i|<=1, hence 0<=E[C4_i]<=scale_i by the Ginibre mean result.
    return {
      'id':'R7P-078-monotone-iteration',
      'status':'EXACT_PROVED_CONDITIONAL_ON_R7P_075',
      'map':'T_i(s)=g E_s[C4_i], s in R_+^4, fixed g>0',
      'premises':['R7P-075 gives T(s)>=0 and DT(s)>=0 entrywise on R_+^4','T is continuous (finite exponential family)'],
      'supersolution':'B=(g*sqrt(lambda3/6), g*sqrt(lambda4/6), g*sqrt(lambda5/6), g*sqrt(lambda6/12)). These are g times the coordinatewise feature maxima. Since E[C4_i]<=max C4_i, T(B)<=B and in fact T(s)<=B for every s>=0.',
      'subsolution':'0 is a fixed point because every retained nonconstant Fourier feature has zero uniform mean, so T(0)=0.',
      'trapping_interval':'[0,B] is invariant under T.',
      'upper_iteration':'s^(0)=B, s^(n+1)=T(s^(n)) is componentwise nonincreasing and bounded below by 0, hence converges coordinatewise.',
      'limit':'Continuity gives T(s^*)=s^*. For every fixed point y in [0,B], y<=s^(n) for all n, so y<=s^*: the limit is the greatest fixed point in [0,B].',
      'lower_iteration':'Starting at 0 stays at the minimal fixed point 0.',
      'algorithm_warning':'This is a convergence theorem for the declared fixed-point iteration only, not for unspecified physical dynamics.',
    }


def write_records():
    out={'R7P-075':theorem_record(),'R7P-078':monotone_iteration_record()}
    p=ROOT/'results/R7P-075_078_ginibre_monotonicity.json'
    p.write_text(json.dumps(out,indent=2)+'\n')
    return out


if __name__=='__main__':
    print(json.dumps(write_records(),indent=2))
