from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import ndtr
from observation_core import *
OUT=Path(__file__).resolve().parent

def main():
    results={"scope":"numerical tests of exact finite formulas and illustrative limiting integrals; no physical records", "q":Q,"chi":CHI,"gc_nominal":GC}
    # Exact hypergeom identity with two disjoint blocks, without probability conditioning mistakes.
    static=[]
    for N,M,k in ((12,6,4),(24,14,6),(120,66,12),(144,66,12)):
        ps=[hypergeom_pmf(N,M,k,h) for h in range(k+1)]
        empirical_mean=sum(h*p for h,p in enumerate(ps))
        empirical_var=sum((h-empirical_mean)**2*p for h,p in enumerate(ps))
        mom=hypergeom_moments(N,M,k)
        static.append(dict(N=N,M=M,k=k,support=[max(0,k-(N-M)),min(k,M)],sum=sum(ps),mean_error=empirical_mean-mom['mean_h'], variance_error=empirical_var-mom['var_h']))
    results['exact_hypergeom_checks']=static
    jchecks=[]
    for N,M,k1,k2 in ((12,6,3,4),(24,14,4,6),(120,66,12,11)):
        probs=np.array([[joint_hypergeom_pmf(N,M,k1,k2,x,y) for y in range(k2+1)] for x in range(k1+1)])
        h1=np.arange(k1+1);h2=np.arange(k2+1)
        mean1=sum(probs.sum(1)*h1);mean2=sum(probs.sum(0)*h2)
        cov_h=sum(probs[x,y]*(x-mean1)*(y-mean2) for x in range(k1+1) for y in range(k2+1))
        predicted=k1*k2*(-(M/N)*(1-M/N)/(N-1))
        jchecks.append(dict(N=N,M=M,k1=k1,k2=k2,total=float(probs.sum()),cov_h=cov_h,cov_error=cov_h-predicted))
    results['two_block_checks']=jchecks
    # Compare direct labelled generator and count projection / predicted tagged forcing bound.
    labelled=[]
    rng=np.random.default_rng(1187)
    for N,k in ((12,1),(24,2),(120,3),(120,12)):
        lab=rng.integers(0,Q,size=N)
        n=np.bincount(lab,minlength=12);g=GC
        a=labelled_stats(lab,k,g)
        a['count_projection_error']=a['Lm']-generator_drift_m(n,g)
        a['count_bracket_error']=a['bracket_m']-generator_bracket_m(n,g)
        a['forcing_times_N']=a['forcing']*N
        labelled.append(a)
    results['labelled_generator_checks']=labelled
    # Symmetric prior, median/sign decision probability for one frozen alpha=1.
    norm=quad(lambda z:float(quartic_density_unnorm(z)), -np.inf,np.inf,epsabs=2e-13,epsrel=2e-13)[0]
    out={}
    for alpha in (0.25,1.0,4.0):
        risk=2*quad(lambda z: math.exp(-z**4/12)*ndtr(-math.sqrt(alpha)*z),0,np.inf,epsabs=3e-13,epsrel=3e-13)[0]/norm
        out[str(alpha)]=risk
    results['bounded_sign_risk_nominal']={'distribution':'quartic lambda0', 'alpha':out,'risk_sign_tie':'zero in continuous limit', 'method':'scipy.quad + scipy.special.ndtr; NOT outward interval'}
    results['critical_private_ou_increment_variance_alpha1_delta_modeltime1']=2*(1-math.exp(-1))
    # Record-bracket target for fixed k at uniform after g0 matched.
    results['integrated_fixed_tag_bracket_target']={str(k):2/k for k in (1,2,3)}
    (OUT/'NUMERICAL_RESULTS.json').write_text(json.dumps(results,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    print('wrote NUMERICAL_RESULTS.json; gc',GC,'risk alpha=1',out['1.0'])
if __name__=='__main__':main()
