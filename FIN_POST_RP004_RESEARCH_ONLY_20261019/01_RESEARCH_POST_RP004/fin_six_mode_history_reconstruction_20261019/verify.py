#!/usr/bin/env python3
"""Bounded independent algebra/numeric checks for a scoped FIN Z12 inverse correlation theorem.
Does not generate a physics law or validate measurement data.
"""
import json
import math
from pathlib import Path
import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
r3 = sp.sqrt(3)
B = sp.Matrix([
    [sp.simplify(2*(1-sp.cos(sp.pi*sp.Rational(k*d,6))) if d < 6 else 1-(-1)**k)
     for d in range(1,7)] for k in range(1,7)
])
B_inv = B.inv().applyfunc(lambda x: sp.simplify(sp.radsimp(x)))
w = np.array([
    .46998567264502017, .19204355169010282, .09142861427792497,
    .04702916874565040, .02413122336363006, .011070817321442113
], dtype=float)
Bnum = np.array(B.evalf(17).tolist(),dtype=float)
inverse_num=np.array(B_inv.evalf(17).tolist(),dtype=float)
lams = Bnum@w
W=np.zeros((12,12))
for i in range(12):
    for j in range(12):
        if i!=j:
            d=min(abs(i-j),12-abs(i-j))
            W[i,j]=w[d-1]
A=np.diag(W.sum(axis=1))-W
normed=[]
for k in range(1,7):
    vec=np.array([math.sqrt(2)*math.cos(2*math.pi*k*j/12) if k<6 else (-1)**j for j in range(12)],dtype=float)
    normed.append(vec)
assert sp.simplify(sp.radsimp(B.det())+3456*r3)==0
assert (B*B_inv-sp.eye(6)).applyfunc(sp.simplify)==sp.zeros(6,6)
assert B[:5,:].rank()==5
null=sp.Matrix([-1,1,-1,1,-1,1])
assert B[:5,:]*null==sp.zeros(5,1)
assert B*null==sp.Matrix([0,0,0,0,0,-12])
assert min(w)>0
assert np.max(np.abs(A@np.ones(12)))<1e-13
assert np.min(np.linalg.eigvalsh(A))>-1e-12
assert np.max(np.abs(A-A.T))==0
assert np.max(np.abs(inverse_num@lams-w))<2e-15
assert np.max(np.abs(lams-np.linalg.eigvalsh(A)[[1,3,5,7,9,11]]))<1e-12
# The cosine k=6 is the alternating mode.
for k,vec in enumerate(normed,1):
    assert abs(vec@vec/12-1)<1e-13
    assert max(abs(vec@u/12) for u in normed if u is not vec) <1e-12
    assert np.max(np.abs(A@vec-lams[k-1]*vec))<2e-13

# Stationary two-time correlation diagonal formula, computed through matrix exponentials.
from scipy.linalg import expm
lags=(.15,.6,1.5)
correlation=[]
for lag in lags:
    Pt=expm(-lag*A)
    direct=np.array([v@Pt@v/12 for v in normed])
    exact=np.exp(-lag*lams)
    assert np.max(np.abs(direct-exact))<1e-13
    reconstructed = inverse_num@(-np.log(direct)/lag)
    assert np.max(np.abs(reconstructed-w))<2e-12
    correlation.append({'lag':lag,'values':direct.tolist(),'recovered_max_abs_weight_error':float(np.max(np.abs(reconstructed-w)))})

# Explicit nonidentifiability for all time lags of the first five modal observables.
eta=0.004
v=np.array([-1,1,-1,1,-1,1], dtype=float)
w_alt=w+eta*v
assert min(w_alt)>0
lam_alt=Bnum@w_alt
assert np.max(np.abs(lam_alt[:5]-lams[:5]))<2e-15
assert abs((lam_alt[-1]-lams[-1])+12*eta)<2e-15
assert np.max(np.abs(w_alt-w))>0

# Exact coefficient-free cone membership / incompatibility test (H4).
invalid_lam=sp.Matrix([3,1,1,1,1,1])
invalid_rates=(B_inv*invalid_lam).applyfunc(sp.simplify)
assert invalid_rates[0] == sp.Rational(1,12)-sp.sqrt(3)/6
assert invalid_rates[1] == -sp.Rational(1,12)
assert min(float(a) for a in invalid_rates)<0
assert min(float(a) for a in B_inv*sp.Matrix([float(x) for x in lams]))>0

# A deterministic synthetic measurement-error check and the proved uniform stability bound.
t=.6
c=np.exp(-lams*t)
noise=np.array([1.,-1.,.4,-.6,.8,-.25])*1e-5
c_hat=c+noise
assert np.min(c_hat)>0
w_hat=inverse_num@(-np.log(c_hat)/t)
err=float(np.max(np.abs(w_hat-w)))
norm_inf=float(np.max(np.sum(np.abs(inverse_num),axis=1)))
assert abs(norm_inf-11/12)<1e-14
max_noise=float(np.max(np.abs(noise)))
cmin=float(np.min(c))
bound=norm_inf*max_noise/(t*(cmin-max_noise))
assert err<=bound+1e-14

out={
    'type':'STRICT_Z12_REVERSIBLE_CIRCULANT_HISTORY_IDENTIFIABILITY_CONDITIONAL',
    'sympy_version':sp.__version__, 'numpy_version':np.__version__,
    'B':[[str(a) for a in row] for row in B.tolist()],
    'B_inverse':[[str(a) for a in row] for row in B_inv.tolist()],
    'det_exact':str(sp.simplify(sp.radsimp(B.det()))), 'rank':int(B.rank()),
    'source_strict_weights':w.tolist(), 'recovered_rates':lams.tolist(),
    'single_lag_correlation_checks':correlation,
    'log_correlation_cone_H4':{'criterion':'B_inverse @ (-log(C)) > 0 componentwise','counterexample_lambda':[3,1,1,1,1,1],'resulting_inadmissible_weights':[str(sp.simplify(sp.radsimp(a))) for a in invalid_rates],'passes_strict_FIN':True,'rejects_counterexample':True},
    'blind_five_mode_perturbation':{'eta':eta,'vector':v.tolist(),'alternative_weights':w_alt.tolist(),'max_error_observed_five_lambda':float(np.max(np.abs(lam_alt[:5]-lams[:5]))),'sixth_rate_delta':float(lam_alt[-1]-lams[-1])},
    'finite_noise_bound':{'lag':t,'absolute_correlation_noise_max':max_noise,'weight_error_max':err,'upper_bound':bound,'norm_Binv_inf':norm_inf},
    'controls':{'exact_determinant':True,'exact_inverse':True,'five_mode_nullspace':True,'time_correlation_matrix_exp':True,'all_six_mode_identification':True,'one_missing_mode_nonuniqueness':True,'noise_stability_bound':True,'six_correlation_cone_test':True,'no_physical_source_claim':True}
}
(ROOT/'RESULTS.json').write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
print('PASS | det=',out['det_exact'],' | eigenvalues=',lams.tolist())
print('PASS | 5-mode null=',null.T,'| sixth modal rate shift=',float(lam_alt[-1]-lams[-1]))
print('PASS | reconstruction bound',err,'<=',bound)
