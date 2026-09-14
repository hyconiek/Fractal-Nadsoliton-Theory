from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp
from model import feature_spaces, dual7, d12_actions
ROOT=Path(__file__).resolve().parents[1]

LAM_IV={
3:(196140686197643,10**14,39228137239529,2*10**13),
4:(5498922123333,2500000000000,109978442466661,50000000000000),
5:(57465156801977,25000000000000,22986062720791,10000000000000),
6:(234218204114629,100000000000000,234218204114631,100000000000000)}

def ivq(a,b,c,d): return mp.iv.mpf([mp.mpf(a)/b,mp.mpf(c)/d])

def rational_witness():
    mp.iv.dps=70
    nums=[836365227,1700465,10147081,26989852,24859253,13242432,9756607,13242432,24859253,26989852,10147081,1700465]
    den=10**9; p=[mp.iv.mpf(n)/den for n in nums]; pi=mp.iv.pi
    ls={k:ivq(*LAM_IV[k]) for k in LAM_IV}
    S={k:sum(p[j]*mp.iv.cos(2*pi*k*j/12) for j in range(12)) for k in (3,4,5)}
    S[6]=sum(p[j]*((-1)**j) for j in range(12))
    Q=sum(ls[k]/6*S[k]**2 for k in (3,4,5))+ls[6]/12*S[6]**2
    D=sum(x*mp.iv.log(12*x) for x in p)
    ratio=2*D/Q; g=mp.iv.mpf(371835)/100000; E=D-g*Q/2
    def pair(x): return [str(x.a),str(x.b)]
    return {'probability_numerators':nums,'denominator':den,'D_interval':pair(D),'Q_interval':pair(Q),
            'ratio_interval':pair(ratio),'energy_at_g_3p71835_interval':pair(E),
            'strict_negative':bool(E.b<0)}

def g4_gradient_flow():
    W,A,L,X,C,A7=feature_spaces(); actions=d12_actions(X); atlas=json.load(open(ROOT/'results/R7P-037_full7_stationary_atlas.json'))
    row=next(x for x in atlas['gains'] if abs(x['g']-4)<1e-12)
    loc=min((o for o in row['orbits'] if o['H7_index']==0),key=lambda o:o['phi'])
    sad=next(o for o in row['orbits'] if o['H7_index']==1)
    uniform=next(o for o in row['orbits'] if o['stabilizer_size']==24)
    ths=np.array(sad['theta']); H=dual7(ths,4,X)[2]; vals,vecs=np.linalg.eigh(H); v=vecs[:,0]
    def f(t,y): return -dual7(y,4,X)[1]
    rows=[]
    for sg in (-1,1):
        y0=ths+sg*1e-5*v
        sol=solve_ivp(f,(0,600),y0,rtol=1e-10,atol=1e-12,max_step=.2)
        y=sol.y[:,-1]; val,gr,H,p=dual7(y,4,X)
        dl=min(np.linalg.norm(y-T@np.array(loc['theta'])) for P,T in actions.values()); du=np.linalg.norm(y-np.array(uniform['theta']))
        rows.append({'sign':sg,'endpoint':y.tolist(),'phi':val,'residual':float(np.linalg.norm(gr)),
                     'distance_to_localized_rep':float(dl),'distance_to_uniform':float(du),
                     'classification':'uniform' if du<dl else 'localized_orbit_rep'})
    return {'law':'theta_dot=-grad Phi_4(theta)','initial_saddle':sad['theta'],'unstable_eigenvalue':float(vals[0]),
            'branches':rows,'proof_level':'NUMERICAL; no validated trajectory tube'}

def main():
    out={'R7P-097_099':{'definition':'g_global=inf_{p != u,Q>0} 2D(p||u)/Q','Q':'(p-u)^T A7 (p-u)',
        'uniform_limit':'12/lambda_max(A7) along the top active tangent direction; distinct from the finite-amplitude event',
        'lower_bound':{'g':'2.8934','reason':'A7<=A_full and accepted ST448 gives V_full>=0 with unique uniform minimizer, hence V7>=V_full'},
        'upper_witness':rational_witness(),'bracket':['2.8934','3.71835'],
        'nonconclusion':'Does not identify the first attaining orbit or prove uniqueness of the transition.'},
        'R7P-100_101':{'g4_stationary_atlas':'results/R7P-037_full7_stationary_atlas.json',
          'g4_observed_phis':[-0.13884766301180296,0.0,0.022577225011047197],
          'status':'PARTIAL_NUMERICAL_ATLAS_WITH_UNCOVERED_7D_COMPLEMENT',
          'unique_global_orbit_at_g4':'UNRESOLVED'},
        'R7P-102':{'reduction_licenses':{'D12':'exact symmetry action; allows orbit deduplication but does not force a generic minimizer into reflection-even C4','C4_reflection_fixed':'valid only for certified roots/branches in that invariant subspace; no global minimizer forcing theorem','phase_fixture':'fixed-amplitude torus only; amplitude-zero strata and generic 7D states remain outside','two_harmonic':'exact invariant restricted family only'},'unsupported_global_reduction':True},
        'R7P-103':g4_gradient_flow(),
        'R7P-104':{'global_gain_bracket':['2.8934','3.71835'],'fixed_g4_best_found_phi':-0.13884766301180296,
          'global_orbit_uniqueness':'UNRESOLVED','stationary_exhaustion':'UNRESOLVED','dynamical_connections':'NUMERICAL_FOR_DECLARED_EUCLIDEAN_GRADIENT_FLOW_ONLY'}}
    (ROOT/'results/R7P-097_104_global_frontier.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'witness_negative':out['R7P-097_099']['upper_witness']['strict_negative'], 'flow':[x['classification'] for x in out['R7P-103']['branches']]},indent=2))
if __name__=='__main__': main()
