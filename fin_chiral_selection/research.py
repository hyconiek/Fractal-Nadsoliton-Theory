"""Programme ST8651: an exact intrinsic Fourier-spin sector.

The plus-i state convention is the conjugate of the previously audited
minus-i convention and matches exp(-it(sI-K)) up to a scalar phase.
No physical derivation or global mathematical priority is claimed.
"""
from fractions import Fraction as F
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


def strict():
    return np.array([[0. if i==j else math.cos(.18575*min(abs(i-j),12-abs(i-j))+.1625)/
        (1+min(abs(i-j),12-abs(i-j))**1.8) for j in range(12)] for i in range(12)])


def off(M):
    value=np.array(M,copy=True);np.fill_diagonal(value,0);return value


def sector():
    root=math.sqrt(3)/2
    co=[1,root,.5,0,-.5,-root,-1,-root,-.5,0,.5,root]
    si=[0,.5,root,1,root,.5,0,-.5,-root,-1,-root,-.5]
    E=np.array([[co[(i+j)%12]/3 if (i+j)%2 else 0. for j in range(12)] for i in range(12)])
    Fm=np.array([[si[(i+j)%12]/3 if (i+j)%2 else 0. for j in range(12)] for i in range(12)])
    G=(E@Fm-Fm@E)/(2j)
    return E,Fm,G


def rhs(t,state,kappa):
    x,y,u,v,z=state
    return np.array([kappa*(u-x),kappa*(v-y),-2*y*z,2*x*z,-2*(x*v-y*u)])


def energy(state):
    x,y,u,v,z=np.asarray(state)
    return x*x+y*y-2*(x*u+y*v)


def unstable_rate(kappa):
    return (-kappa+np.sqrt(complex(kappa*kappa,8*kappa)))/2


def integrate(seed,kappa=1.,margin=60.):
    rate=unstable_rate(kappa)
    horizon=math.log(1/seed)/rate.real+margin
    initial=np.array([0.,0.,seed,0.,math.sqrt(1-seed*seed)])
    result=solve_ivp(lambda t,y:rhs(t,y,kappa),(0,horizon),initial,
        method='DOP853',rtol=2e-11,atol=min(1e-15,seed*1e-8),max_step=.2)
    if not result.success:raise RuntimeError(result.message)
    x,y,u,v,z=result.y
    radius=np.hypot(x,y)
    valid=radius>seed*1e-5
    angles=np.unwrap(np.angle(x[valid]+1j*y[valid]))
    last=result.y[:,-1]
    return dict(seed=seed,kappa=kappa,horizon=horizon,endpoint=last.tolist(),
        final_radius=float(radius[-1]),final_angle=float(angles[-1]),
        spin_invariant_error=float(np.max(abs(u*u+v*v+z*z-1))),
        final_learning_lag=float(math.hypot(last[2]-last[0],last[3]-last[1])),
        energy_initial=float(energy(initial)),energy_final=float(energy(last)),
        sampled_energy_increase=float(np.max(np.diff(x*x+y*y-2*(x*u+y*v)))))


def exact_algebra():
    import sympy as sp
    rt=sp.sqrt(3)/2
    co=[1,rt,sp.Rational(1,2),0,-sp.Rational(1,2),-rt,-1,-rt,-sp.Rational(1,2),0,sp.Rational(1,2),rt]
    si=[0,sp.Rational(1,2),rt,1,rt,sp.Rational(1,2),0,-sp.Rational(1,2),-rt,-1,-rt,-sp.Rational(1,2)]
    E=sp.Matrix(12,12,lambda i,j:sp.sympify(co[(i+j)%12])/3 if (i+j)%2 else sp.S.Zero)
    Fm=sp.Matrix(12,12,lambda i,j:sp.sympify(si[(i+j)%12])/3 if (i+j)%2 else sp.S.Zero)
    G=(E*Fm-Fm*E)/(2*sp.I)
    def zero(M):return all(sp.simplify(v)==0 for v in M)
    checks=dict(diagonal=all(E[i,i]==Fm[i,i]==G[i,i]==0 for i in range(12)),
        rows=zero(E*sp.ones(12,1)) and zero(Fm*sp.ones(12,1)) and zero(G*sp.ones(12,1)),
        EF=zero(E*Fm-Fm*E-2*sp.I*G),FG=zero(Fm*G-G*Fm-2*sp.I*E),
        GE=zero(G*E-E*G-2*sp.I*Fm),same_square=zero(E*E-Fm*Fm) and zero(E*E-G*G),
        support_projector=zero(E**4-E**2))
    if not all(checks.values()):raise AssertionError(checks)
    weights=sp.symbols('w1:7',real=True)
    W=sp.Matrix(12,12,lambda i,j:0 if i==j else weights[min(abs(i-j),12-abs(i-j))-1])
    checks['commutes_with_any_real_C12_circulant']=all(zero(W*M-M*W) for M in [E,Fm,G])
    checks['pairwise_anticommutators']=all(zero(A*B+B*A) for A,B in [(E,Fm),(Fm,G),(G,E)])
    checks['orthogonal_trace_products']=all(sp.trace(A*B)==0 for A,B in [(E,Fm),(Fm,G),(G,E)])
    checks['support_rank']=int(sp.trace(E*E))
    assert checks['support_rank']==4 and all(checks.values())
    return checks


def exact_scalar_identities():
    import sympy as sp
    x,y,u,v,z,k=sp.symbols('x y u v z k',real=True)
    state=[x,y,u,v,z]
    field=[k*(u-x),k*(v-y),-2*y*z,2*x*z,-2*(x*v-y*u)]
    def D(f):return sp.expand(sum(sp.diff(f,s)*f_s for s,f_s in zip(state,field)))
    radius=x*x+y*y;V=radius-2*(x*u+y*v);spin=u*u+v*v+z*z
    checks=dict(spin=D(spin)==0,
        lyapunov=sp.expand(D(V)+2*k*((u-x)**2+(v-y)**2))==0,
        energy_gap=sp.expand(V+spin-((x-u)**2+(y-v)**2+z*z))==0,
        radius_power=sp.expand(D(D(radius))+k*D(radius)-2*k*k*((u-x)**2+(v-y)**2))==0,
        angular_momentum=sp.expand(k*(x*v-y*u)+k*field[-1]/2)==0)
    assert all(checks.values())
    return checks


def run():
    K0=strict();E,Fm,G=sector();gamma=.05;R=.0005;scale=R/gamma;kappa=1.;eta=kappa*R/gamma**2
    base=np.eye(12)/12+gamma*K0
    residuals=[]
    for state in [[.2,-.1,.3,.4,math.sqrt(.75)],[-.3,.2,.4,-.5,math.sqrt(.59)]]:
        x,y,u,v,z=state
        K=K0+scale*(x*E+y*Fm);rho=base+R*(u*E+v*Fm+z*G)
        a,b,c,d,e=rhs(0,state,kappa)
        reduced_K=scale*scale*(a*E+b*Fm)
        reduced_rho=R*scale*(c*E+d*Fm+e*G)
        full_K=eta*(off(rho.real)-gamma*K)
        full_rho=1j*(K@rho-rho@K)
        residuals.append(float(np.linalg.norm(reduced_K-full_K)+np.linalg.norm(reduced_rho-full_rho)))
    assert max(residuals)<1e-13
    rate=unstable_rate(kappa);slope=-rate.imag/rate.real
    seeds=[10.**(-k) for k in range(2,11)]
    runs=[integrate(seed,kappa) for seed in seeds]
    observed=[(runs[i+1]['final_angle']-runs[i]['final_angle'])/
              math.log(runs[i+1]['seed']/runs[i]['seed']) for i in range(len(runs)-1)]
    assert max(x['spin_invariant_error'] for x in runs)<1e-8
    assert max(abs(x['final_radius']-1) for x in runs)<1e-6
    return dict(status='Exact sector and scalar identities checked; global convergence and endpoint obstruction proved in PROOF.md. Fixed-K0 phase asymptotic remains open.',
        exact_sector_algebra=exact_algebra(),exact_scalar_identities=exact_scalar_identities(),parameters=dict(gamma=gamma,R=R,kappa=kappa,eta=eta,
            final_kernel_radius=scale),full_model_residuals=residuals,
        uniform_population_error=float(np.max(abs(np.diag(base+R*G)-1/12))),
        minimum_base_density_eigenvalue=float(np.linalg.eigvalsh(base)[0]),
        positive_odd_edge_margin=float(K0[0,5]-scale/3),
        unstable_exponent=dict(real=float(rate.real),imaginary=float(rate.imag)),
        asymptotic_phase_slope_prediction=float(slope),observed_log_seed_slopes=observed,
        seed_scale_ratio_prediction=float(math.exp(2*math.pi*rate.real/abs(rate.imag))),runs=runs)


if __name__=='__main__':
    data=run();Path(__file__).with_name('results.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({k:v for k,v in data.items() if k!='runs'},indent=2))
