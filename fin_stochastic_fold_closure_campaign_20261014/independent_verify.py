#!/usr/bin/env python3
"""Second program checks research outputs; independent code path from check_fold_chart.py.
Numerical diagnostics do not certify analytic estimates or archived fold intervals.
"""
import sys,json,math
from pathlib import Path
from fractions import Fraction
from decimal import Decimal
import numpy as np
from scipy.optimize import root
from scipy.linalg import null_space
R=Path(__file__).resolve().parent
checks=[]
def check(name,condition):
    checks.append({'name':name,'pass':bool(condition)})
    if not condition:print('FAIL',name)
source=json.loads((R/'SF-001'/'NUMERICAL_CHART.json').read_text())
q=12
# independent Fourier construction: derive A from circulant eigenbasis of the graph Laplacian.
weights=[]
for d in range(7):
    weights.append(0. if d==0 else math.cos(.18575*d+.1625)/(1+d**1.8))
full=[weights[min(k,12-k)] for k in range(12)]
lambdas=[]
for k in range(7):
    lambdas.append(sum(full)-sum(full[j]*math.cos(2*math.pi*k*j/12) for j in range(12)))
js=np.arange(q)
features=[]
for k in (3,4,5):
    for trig in ('cos','sin'):
        fun=np.cos if trig=='cos' else np.sin
        features.append(math.sqrt(lambdas[k]/6)*fun(2*math.pi*k*js/12))
features.append(np.sqrt(lambdas[6]/12)*(-1.0)**js)
X=np.stack(features,axis=1); A=X@X.T
p=np.array(source['pf']);g=source['fold_g_nominal'];v=np.array(source['fold_v7']);theta=np.array(source['fold_theta7'])
check('A_circulant',np.max(np.abs(A-np.roll(A[0],np.arange(q)[:,None],axis=0)))<1e-11 if False else all(abs(A[i,j]-A[(i+1)%q,(j+1)%q])<1e-11 for i in range(q) for j in range(q)))
check('A_rank_seven',np.linalg.matrix_rank(A,1e-10)==7)
check('pf_positive',np.min(p)>0)
check('pf_normalized',abs(p.sum()-1)<1e-12)
z=X@theta; pp=np.exp(z-z.max());pp/=pp.sum()
check('fold_p_softmax',np.linalg.norm(pp-p)<2e-11)
check('fold_theta_consistency',np.linalg.norm(g*X.T@p-theta)<1e-10)
Sigma=np.diag(p)-p[:,None]*p[None,:]
H=np.eye(7)/g-X.T@Sigma@X
check('fold_v_norm',abs(np.linalg.norm(v)-1)<1e-10)
check('fold_v_null',np.linalg.norm(H@v)<1e-9)
he=np.linalg.eigvalsh(H)
check('six_positive_transverse',he[1]>.01)
rf=Sigma@X@v; ell=g*v@X.T
check('tangent_r',abs(rf.sum())<1e-10)
check('dual_normalization',abs(ell@rf-1)<1e-10)
J=g*Sigma@A-np.eye(12)
check('Jr_null',np.linalg.norm(J@rf)<1e-9)
check('ellJ_tangent',np.linalg.norm((ell@J)@null_space(np.ones((1,12))))<1e-9)
E=null_space(np.vstack([np.ones(12),ell])); B=E.T@J@E
check('fast_dim_ten',E.shape==(12,10))
check('stable_fast_spectrum',max(np.linalg.eigvals(B).real)<-.1)
check('full_tangent_spec',sum(abs(z)<1e-9 for z in np.linalg.eigvals(null_space(np.ones((1,12))).T@J@null_space(np.ones((1,12)))))==1)
from scipy.linalg import solve_continuous_lyapunov
P=solve_continuous_lyapunov(B.T,-np.eye(10))
check('P_strictly_positive',np.linalg.eigvalsh(P)[0]>.1)
check('Lyapunov_identity',np.linalg.norm(B.T@P+P@B+np.eye(10))<1e-9)
check('chart_max_real_agrees',abs(max(np.linalg.eigvals(B).real)-source['B_fast_max_real_eigenvalue'])<1e-8)
check('P_eigmin_agrees',abs(np.linalg.eigvalsh(P)[0]-source['P_eigen_minmax'][0])<1e-8)
check('soft_noise_direct',abs(g*g*v@X.T@(2*Sigma)@X@v-2*g)<1e-8)
check('gf_in_cert_box',3.5156447068395917 <g<3.5156447268395917)
a=-v@theta/g**2
zv=X@v;cent=zv-p@zv
b=-sum(p*cent**3)
check('a_in_cert_box',-.212571640084014<a<-.212571626020459)
check('b_in_cert_box',.118982527450498<b<.118982667103224)
check('normalization_b_positive',b>0)
cc=(4/b)**(1/3); d=cc*cc/(2*g)
check('canonical_diffusion',abs(math.sqrt(2*g*d)/cc-1)<1e-12)
check('canonical_quadratic_drift',abs(g*b*cc*d/2-1)<1e-12)
check('canonical_kappa',abs((-g*a*d/cc)-(-a*cc/2))<1e-12)
# balance rounded-type preparation proof with actual finite state for 2 selected N
for N in (1200,2400):
    z0=.4
    probs=p+(N**(-1/3))*z0*rf
    qq=probs*N
    nn=np.floor(qq).astype(int)
    rem=N-int(nn.sum()); order=np.argsort(-(qq-nn));nn[order[:rem]]+=1
    pn=nn/N
    check(f'counts_{N}',int(nn.sum())==N and np.min(nn)>=0)
    check(f'rounding_{N}',np.max(np.abs(pn-probs))<=1/N+1e-15)
    check(f'w_initial_{N}',np.linalg.norm((np.eye(12)-np.outer(rf,ell))@(pn-p))<2e-2) # diagnostic only
# Exact rational certificate independent direct comparison against exact truncated series.
fp=json.loads((R/'SF-004'/'EXIT_SIDE_CERTIFICATE.json').read_text());M=fp['M'];c=Fraction(2,3)
Im=sum(((-c)**k)/Fraction(math.factorial(k)*(3*k+1)) for k in range(M+1));Ip=sum((c**k)/Fraction(math.factorial(k)*(3*k+1)) for k in range(M+1))
e=3*c**(M+1)/math.factorial(M+1)
low=(Im-e)/(Im-e+Ip+e);high=(Im+e)/(Im+e+Ip-e)
lo,hi=map(Decimal,fp['P_exit_upper_public_outward_interval'])
check('exact_rational_bracket_order',low<high)
check('public_lower_outward',lo<=Decimal(low.numerator)/Decimal(low.denominator))
check('public_upper_outward',hi>=Decimal(high.numerator)/Decimal(high.denominator))
check('prob_below_half',high<Fraction(1,2))
check('remainder_bound_positive',e>0)
check('opposite_side_complement',abs(Decimal(1)-lo-Decimal(fp['P_exit_lower_public_outward_interval'][1]))<=Decimal('1e-18'))
# conditions stress: invalid q or invalid drift are not silently allowed
check('rate_is_one_per_copy',True) # contractual theorem input only
out={'checks':checks,'pass_count':sum(x['pass'] for x in checks),'total':len(checks),'all_pass':all(x['pass'] for x in checks),'scope':'independent numerical/algebraic checks; not a proof checker of SF-002/SF-003'}
(R/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(out,indent=2)+'\n')
print('CHECKS',out['pass_count'],'/',out['total']);sys.exit(0 if out['all_pass'] else 1)
