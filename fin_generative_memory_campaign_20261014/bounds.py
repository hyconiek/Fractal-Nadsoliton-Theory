"""A-priori/reproducible, source-generator-based truncation envelope. Does NOT look at reserved true outputs."""
from model import *
from scipy.special import gammainc

def basis_matrix(m,B): return np.column_stack((m['U'].toarray(),B))
def exp_model(T,t,w): return expm(T*t)@w

def norm_error_analytic(m,V,T,z,w,t,nterms=20):
    """For spec(S), spec(T) subset [-3,0]: shifted Taylor plus Poisson remainder.
    ||exp(t S)z - V exp(tT)w||₂ ≤ e^-1.5t Σ_k=0^d t^k/k! ||(S+1.5I)^k z - V(T+1.5I)^k w||₂
          + (||z||₂+||w||₂) Pr{Pois(1.5t)>d}.
    Rounding allowances are documented separately.
    """
    zk=np.asarray(z).copy(); wk=np.asarray(w).copy(); coeff=math.exp(-1.5*t)
    residual=[];sm=0.
    for k in range(nterms+1):
        r=float(np.linalg.norm(zk-V@wk));residual.append(r)
        sm+=coeff*r
        zk=m['S']@zk+1.5*zk
        wk=T@wk+1.5*wk
        coeff*=t/(k+1)
    tail=(np.linalg.norm(z)+np.linalg.norm(w))*poisson_tail(1.5*t,nterms)
    return {'l2_error_upper_formula':float(sm+tail),'finite_sum':float(sm),'tail_formula':float(tail),'remainder_degree':nterms,'moment_residual_norms':residual}

def bound_onetime(m,V,T,z,t):
    w=V.T@z
    return norm_error_analytic(m,V,T,z,w,t)

def bound_twotime(m,V,T,z,t1,t2):
    w=V.T@z
    b1=norm_error_analytic(m,V,T,z,w,t1)
    w1=exp_model(T,t1,w);v1=V@w1
    ve=np.zeros(N);ve[:12]=v1[:12]
    wE=V.T@ve
    leak=float(np.linalg.norm(ve-V@wE))
    wN=w1-wE
    E=norm_error_analytic(m,V,T,V@wE,wE,t2-t1)
    NE=norm_error_analytic(m,V,T,V@wN,wN,t2-t1)
    return {'TV_upper_formula':float(b1['l2_error_upper_formula']+leak+(E['l2_error_upper_formula']+NE['l2_error_upper_formula'])/2),'first_stage_l2':b1['l2_error_upper_formula'],'filter_leakage_l2_each_branch':leak,'second_stage_l2_E':E['l2_error_upper_formula'],'second_stage_l2_other':NE['l2_error_upper_formula'],'readout':'at first time Y=(0,0) or complement, at second time all 144 Y'}

if __name__=='__main__':
    m=build(A7());B,n=orthogonal_basis(m,6,True);V=basis_matrix(m,B);T,wi,wr,val=reduced(m,B)
    print('basis',B.shape)
    for prep in prep_names()+['delta_2']:
        z=z_from_p(m,make_prep(m,prep))
        for t in (.35,.8):
            d=bound_onetime(m,V,T,z,t)
            print(prep,t,'one-time theoretical TV upper',.5*d['l2_error_upper_formula'],'tail',d['tail_formula'],'last residual',max(d['moment_residual_norms'][:8]),flush=True)
        d=bound_twotime(m,V,T,z,.35,.8)
        print('2time bound',prep,d['TV_upper_formula'],flush=True)
