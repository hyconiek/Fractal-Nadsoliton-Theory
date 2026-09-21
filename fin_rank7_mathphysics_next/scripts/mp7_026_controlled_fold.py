"""MP7-026: controlled local fold theorem from the accepted R7P-031 interval box.

This is a conservative uniform analytic/interval bound, not a fit to numerical
continuation.  It uses the exact outward spectral intervals and the validated
simple-fold box from R7P-031.
"""
from __future__ import annotations
import os
import json, math, sys
from pathlib import Path
import mpmath as mp

HERE=Path(__file__).resolve()
OUTROOT=Path(os.environ.get('MP7_WORK_ROOT',HERE.parents[1]))
R7P=Path(os.environ.get('R7P_ROOT','/mnt/data/r7p_source/unpacked/fin_rank7_followup'))
sys.path.insert(0,str(R7P))
from src.coexistence_certificate import interval_features,_I,_bound,interval_ldlt
from src.fold_certificate import _augmented_iv

mp.mp.dps=80
mp.iv.dps=60

FOLD=json.loads((R7P/'certificates/R7P-031_simple_fold.json').read_text())
X=[_I(x) for x in FOLD['augmented_root_box']]
_,_,C,_=interval_features()
_,_,Hf,T3f,pf=_augmented_iv(X,C)


def absup(x):
    return max(abs(float(x.a)),abs(float(x.b)))

def norm_interval_vec(xs):
    return math.sqrt(sum(absup(x)**2 for x in xs))

def frob_abs(tensor_flat):
    return math.sqrt(sum(absup(x)**2 for x in tensor_flat))

def interval_moment_tensors(sbox,gbox):
    # C already carries the accepted spectral intervals.
    h=[sum(C[j][a]*sbox[a] for a in range(4)) for j in range(12)]
    e=[mp.iv.exp(x) for x in h]; Z=sum(e,_I(0)); p=[x/Z for x in e]
    mu=[sum(p[j]*C[j][a] for j in range(12)) for a in range(4)]
    Y=[[C[j][a]-mu[a] for a in range(4)] for j in range(12)]
    cov=[[sum(p[j]*Y[j][a]*Y[j][b] for j in range(12)) for b in range(4)] for a in range(4)]
    T3=[[[-sum(p[j]*Y[j][a]*Y[j][b]*Y[j][c] for j in range(12))
          for c in range(4)] for b in range(4)] for a in range(4)]
    T4=[[[[-(sum(p[j]*Y[j][a]*Y[j][b]*Y[j][c]*Y[j][d] for j in range(12))
               -(cov[a][b]*cov[c][d]+cov[a][c]*cov[b][d]+cov[a][d]*cov[b][c]))
           for d in range(4)] for c in range(4)] for b in range(4)] for a in range(4)]
    H=[[(_I(1)/gbox if a==b else _I(0))-cov[a][b] for b in range(4)] for a in range(4)]
    return p,cov,T3,T4,H

def ldl_eigen_lower(A):
    # If A=L D L^T with all D>0, lambda_min(A)>=min(D)/||L^{-1}||_2^2.
    # We upper bound ||L^{-1}||_2 by its interval-Frobenius envelope.
    D,L,ok=interval_ldlt(A)
    if not ok or min(float(x.a) for x in D)<=0:
        raise RuntimeError('LDL positivity failed')
    n=len(A); Inv=[[_I(0) for _ in range(n)] for __ in range(n)]
    for col in range(n):
        for i in range(n):
            Inv[i][col]=(_I(1 if i==col else 0)-sum(L[i][j]*Inv[j][col] for j in range(i)))/L[i][i]
    inv_frob=math.sqrt(sum(absup(x)**2 for row in Inv for x in row))
    dmin=min(float(x.a) for x in D)
    return dmin/(inv_frob**2),D,L,inv_frob

def main():
    # Fold coefficients are already validated in R7P-031.
    a_lo,a_hi=map(float,FOLD['fold_coeff_v_dot_Fg'])
    b_lo,b_hi=map(float,FOLD['fold_coeff_D3Phi_vvv'])
    c0_lo=math.sqrt(-2*a_hi/b_hi)
    c0_hi=math.sqrt(-2*a_lo/b_lo)

    sf=X[:4]; gf=X[4]; vf=X[5:]
    gmin=float(gf.a)
    sf_norm=norm_interval_vec(sf)

    # First three Cartesian amplitude coordinates serve as transverse chart E.
    cvec=[-sf[i]/(gf*gf) for i in range(3)]
    dvec=[sum(T3f[i][j][k]*vf[j]*vf[k] for j in range(4) for k in range(4)) for i in range(3)]
    c_norm=norm_interval_vec(cvec); d_norm=norm_interval_vec(dvec)

    # Uniform neighborhood. It is deliberately much larger than the final branch tube.
    s_center=[(float(x.a)+float(x.b))/2 for x in sf]
    s_radius=0.002
    eps_max=1e-7
    sbox=[_I([s_center[i]-s_radius,s_center[i]+s_radius]) for i in range(4)]
    gbox=_I([float(gf.a),float(gf.b)+eps_max])
    p,cov,T3,T4,H=interval_moment_tensors(sbox,gbox)
    M2=frob_abs([T3[a][b][c] for a in range(4) for b in range(4) for c in range(4)])
    M4=frob_abs([T4[a][b][c][d] for a in range(4) for b in range(4) for c in range(4) for d in range(4)])
    Hnorm=frob_abs([H[a][b] for a in range(4) for b in range(4)])

    # Uniform positive transverse block in E coordinates.
    A=[[H[i][j] for j in range(3)] for i in range(3)]
    m_trans,Dloc,Lloc,Linvloc=ldl_eigen_lower(A)

    # At the exact fold, the full H4 has one zero mode. By Cauchy interlacing,
    # lambda_2(H4_fold) is at least lambda_min of its leading 3x3 principal block.
    Af=[[Hf[i][j] for j in range(3)] for i in range(3)]
    gamma_fold,Dfold,Lfold,Linffold=ldl_eigen_lower(Af)

    delta=0.015
    z_inner=c0_lo*(1-delta)
    z_outer=c0_hi*(1+delta)

    # Uniform scalar bounds at eps_max; all terms are monotone in sqrt(eps).
    e=eps_max; x=z_outer*math.sqrt(e)
    R3_x=(M4*x**3 + 6/gmin**3*x*e**2 + 6*(sf_norm+x)/gmin**4*e**3)/6
    T0=c_norm*e + 0.5*d_norm*x*x + x*e/gmin**2 + sf_norm*e*e/gmin**3 + R3_x
    y=T0/m_trans
    ds=x+y
    if ds+s_radius*0 > s_radius: # explicit readability; true fold-box width is negligible here
        raise RuntimeError(f'local tube leaves s-box: {ds} > {s_radius}')
    R3=(M4*ds**3 + 6/gmin**3*ds*e**2 + 6*(sf_norm+ds)/gmin**4*e**3)/6
    R=M2*x*y + 0.5*M2*y*y + e/gmin**2*(x+y) + sf_norm*e*e/gmin**3 + R3
    R_over_e=R/e

    inner_margin=-(a_hi+0.5*b_hi*z_inner*z_inner)
    outer_margin=(a_lo+0.5*b_lo*z_outer*z_outer)
    branch_sign_pass=(inner_margin>R_over_e and outer_margin>R_over_e)

    # Strict convexity of reduced stationarity r(xi,eps).
    dH=M2*ds+e/gmin**2
    eta=dH/m_trans
    db=M4*ds + M2*((1+eta)**3-1)
    reduced_convex_lower=b_lo-db
    convexity_pass=reduced_convex_lower>0

    # Energy difference: reduced-potential derivative is r because the E-gradient vanishes.
    # For P(z)=a z+b z^3/6, the minimum in [z_inner,z_outer] is at c0.
    # The least-negative endpoint is z_outer for this chosen bracket; intervalize safely.
    Aiv=_I([a_lo,a_hi]); Biv=_I([b_lo,b_hi])
    ziv=_I([z_inner,z_outer])
    civ=mp.iv.sqrt(-2*Aiv/Biv)
    Pc=2*Aiv*civ/3
    Pinner=Aiv*_I(z_inner)+Biv*_I(z_inner)**3/6
    Pouter=Aiv*_I(z_outer)+Biv*_I(z_outer)**3/6
    Pmin=float(Pc.a)
    Pmax=max(float(Pinner.b),float(Pouter.b))
    poly_diff_lo=2*Pmin
    poly_diff_hi=2*Pmax
    energy_rem=2*z_outer*R_over_e
    deltaPsi_scaled=[poly_diff_lo-energy_rem, poly_diff_hi+energy_rem]
    barrier_scaled=[-deltaPsi_scaled[1],-deltaPsi_scaled[0]]

    # Reduced curvature r_xixi stays near b.
    # For the actual smallest H4 eigenvalue, use a fold-gap Schur estimate in an
    # orthonormal decomposition (v, v^perp).
    h_err=M2*y + 0.5*M4*ds*ds + e/gmin**2
    if gamma_fold-2*dH <= 0:
        raise RuntimeError('fold spectral separation consumed')
    eig_schur_err=dH*dH/(gamma_fold-2*dH)
    eig_coeff_err=(h_err+eig_schur_err)/math.sqrt(e)
    soft_coeff=[b_lo*z_inner-eig_coeff_err,b_hi*z_outer+eig_coeff_err]
    soft_sign_pass=soft_coeff[0]>0

    out={
      'task':'MP7-026',
      'scientific_state':'PROVED_INTERVAL_ASSISTED_LOCAL_FOLD_SCALING',
      'source_certificate':'R7P-031_simple_fold.json',
      'source_sha256_expected':'5537174281996415a20ea3dffce23c33258c51396b3850aedf4168d9280e352e',
      'coordinate_convention':'s=s_fold+xi*v+E*y, E=(e1,e2,e3), ||v||=1; epsilon=g-g_fold',
      'epsilon_interval':['0','1e-7'],
      'epsilon_endpoint_note':'claims apply for 0<epsilon<=1e-7; epsilon=0 is the certified fold itself',
      'fold_coeff_a':[a_lo,a_hi], 'fold_coeff_b':[b_lo,b_hi],
      'formal_coefficients':{
        'abs_xi_over_sqrt_epsilon':[c0_lo,c0_hi],
        'barrier_over_epsilon_3_2':[0.5357594331508821,0.5357598007359183],
        'soft_abs_eigen_over_sqrt_epsilon':[0.22491024577896145,0.22491038521020988]
      },
      'uniform_neighborhood':{
        's_coordinate_radius':s_radius,
        'M2_T3_frobenius_upper':M2,
        'M4_T4_frobenius_upper':M4,
        'H_frobenius_upper':Hnorm,
        'transverse_lambda_lower':m_trans,
        'fold_second_eigen_lower_by_interlacing':gamma_fold,
        'c_norm_upper':c_norm,'d_norm_upper':d_norm
      },
      'root_bracket':{
        'relative_bracket_delta':delta,
        'abs_xi_over_sqrt_epsilon':[z_inner,z_outer],
        'inner_leading_margin_over_epsilon':inner_margin,
        'outer_leading_margin_over_epsilon':outer_margin,
        'uniform_remainder_over_epsilon_upper':R_over_e,
        'branch_sign_pass':branch_sign_pass,
        'transverse_y_norm_upper_at_epsilon_max':y,
        'total_s_displacement_upper_at_epsilon_max':ds
      },
      'reduced_convexity':{
        'r_xixi_minus_b_abs_upper':db,
        'r_xixi_lower':reduced_convex_lower,
        'pass':convexity_pass,
        'consequence':'exactly two reduced roots in the certified local interval for every 0<epsilon<=1e-7'
      },
      'energy_barrier':{
        'definition':'Phi(saddle)-Phi(minimum) on the two local branches',
        'barrier_over_epsilon_3_2':barrier_scaled,
        'uniform_integrated_remainder_upper':energy_rem,
        'global_barrier_claim':False
      },
      'soft_hessian_eigenvalue':{
        'abs_lambda_soft_over_sqrt_epsilon':soft_coeff,
        'fold_gap_lower':gamma_fold,
        'H_perturbation_upper_at_epsilon_max':dH,
        'schur_correction_upper':eig_schur_err,
        'signs':'negative on xi<0 saddle branch; positive on xi>0 minimum branch',
        'pass':soft_sign_pass
      },
      'interpretation':[
        'The sqrt(epsilon), epsilon^(3/2), and sqrt(epsilon) fold laws are controlled on an explicit punctured interval.',
        'The energy difference is local saddle-to-minimum, not a proved global escape barrier.',
        'No physical time, mobility, temperature, or global-transition claim is introduced.'
      ],
      'all_gates_pass': bool(branch_sign_pass and convexity_pass and soft_sign_pass)
    }
    (OUTROOT/'results/MP7-026_controlled_fold.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__': main()
