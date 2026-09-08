"""Final research batch: event laws, quantum counterexample and clock records.

Analytic proofs are in REPORT.md. No PDF or plotting pipeline.
"""
from fractions import Fraction as F
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from research import strict,lap,noise_curve,wasserstein_one,selfsimilar_moments


def tensor(matrix,N):
    value=matrix
    for _ in range(1,N):value=np.kron(value,matrix)
    return value


def intensity_generator(atoms,N):
    n=len(atoms[0][0]);out=np.zeros((n**N,n**N))
    for R,weight in atoms:out+=weight*(tensor(R,N)-np.eye(n**N))
    return out


def independent_generator(Q,N):
    n=len(Q);out=np.zeros((n**N,n**N))
    for k in range(N):
        term=np.array([[1.]])
        for j in range(N):term=np.kron(term,Q if k==j else np.eye(n))
        out+=term
    return out


def finite_recursive_law(ratio,depth):
    ratio=F(ratio);a=1-ratio;law={F(0):F(1)}
    for _ in range(depth):
        nxt={}
        for u,p in law.items():
            for sign in [-1,1]:
                v=ratio*u+sign*a;nxt[v]=nxt.get(v,F(0))+p/2
        law=nxt
    return sorted(law.items())


def recover_two_atom_measure(moments):
    m0,m1,m2,m3,m4=moments
    determinant=m0*m2-m1*m1
    s1=(m0*m3-m1*m2)/determinant
    s2=(m1*m3-m2*m2)/determinant
    discriminant=s1*s1-4*s2
    num=math.isqrt(discriminant.numerator);den=math.isqrt(discriminant.denominator)
    assert F(num*num,den*den)==discriminant
    low=(s1-F(num,den))/2;high=(s1+F(num,den))/2
    highweight=(m1-low*m0)/(high-low);lowweight=m0-highweight
    nullnorm=m4-2*s1*m3+(s1*s1+2*s2)*m2-2*s1*s2*m1+s2*s2*m0
    assert nullnorm==0
    return [(low,lowweight),(high,highweight)],nullnorm


def phase_coefficient(k,harmonic,epsilon,sign):
    if k==0:return F(1)
    if abs(k)==harmonic:return sign*F(epsilon)/2
    return F(0)


def phase_channel_matrix(N,harmonic,epsilon,sign):
    weights=[int(i).bit_count() for i in range(2**N)]
    return np.array([[float(phase_coefficient(i-j,harmonic,epsilon,sign))
                      for j in weights] for i in weights])


def legacy_cover():
    V=np.zeros((12,12))
    for i in range(12):
        for j in range(12):
            if i!=j:
                d=min(abs(i-j),12-abs(i-j))
                V[i,j]=4*math.log(2)*math.cos(math.pi*d/4+math.pi/6)/(1+.01*d)
    positive=np.maximum(V,0);negative=np.maximum(-V,0)
    W=np.block([[positive,negative],[negative,positive]])
    s=float(W[0].sum());R=W/s;D=np.zeros_like(R)
    h=min(R[0,1]/2,R[0,6]/4)
    for sheet in [0,12]:
        for i in range(12):
            for j in range(12):
                d=min(abs(i-j),12-abs(i-j))
                if d==1:D[sheet+i,sheet+j]=h
                elif d==6:D[sheet+i,sheet+j]=-2*h
    for u in [-1.,1.]:
        assert (R+u*D).min()>=0
        assert np.max(abs((R+u*D).sum(axis=1)-1))<1e-13
    return W,s,R,D


def run():
    W=strict();Q=-lap(W);s,R,D=noise_curve(W)
    base=[(R,s)];with_idle=base+[(np.eye(12),5.)]
    idle_error=float(np.max(abs(intensity_generator(base,2)-intensity_generator(with_idle,2))))
    assert idle_error==0
    smallQ=np.array([[-1.,1.],[1.,-1.]])
    erosion=[]
    Gind=independent_generator(smallQ,3)
    for eps in [.1,.03,.01,.003]:
        G=intensity_generator([(np.eye(2)+eps*smallQ,1/eps)],3)
        error=float(np.linalg.norm(G-Gind,2))
        c=float(np.linalg.norm(smallQ,2))
        upper=sum(math.comb(3,k)*eps**(k-1)*c**k for k in range(2,4))
        assert error<=upper+1e-12
        erosion.append(dict(epsilon=eps,event_rate=1/eps,error=error,upper_bound=upper,
                            pair_joint_rate=eps))
    # Finite weighted-measure reconstruction witness for the general theorem.
    points=[F(1,4),F(3,4)];masses=[F(2),F(1)]
    weighted=[rate*(2*u)**2 for u,rate in zip(points,masses)]
    moments=[sum(rate*(2*u)**2*u**k for u,rate in zip(points,masses)) for k in range(5)]
    assert weighted==[F(1,2),F(9,4)]
    assert [w/(2*u)**2 for u,w in zip(points,weighted)]==masses
    recovered,nullnorm=recover_two_atom_measure(moments)
    assert recovered==list(zip(points,weighted))
    recovered_masses=[weight/(2*u)**2 for u,weight in recovered]
    total_rate=F(13,8)
    independent_rate=total_rate-sum(weight*u for u,weight in zip(points,recovered_masses))
    assert independent_rate==F(3,8)
    eps=F(3,4);harmonic=3
    quantum=[]
    for N in [1,2,3]:
        plus=phase_channel_matrix(N,harmonic,eps,1)
        minus=phase_channel_matrix(N,harmonic,eps,-1)
        gap=float(np.max(abs(plus-minus)))
        assert gap==0 if N<3 else gap==float(eps)
        quantum.append(dict(N=N,channel_multiplier_difference=gap))
    # Nondegenerate preparation/detector ambiguity with identical heat records.
    u=np.ones(12)/12;e=np.eye(12)[0];p1=.6*e+.4*u;p2=.75*e+.25*u
    M1=.8*np.eye(12)+.2*np.ones((12,12))/12
    M2=.64*np.eye(12)+.36*np.ones((12,12))/12
    detector_errors=[float(np.linalg.norm(M1@expm(t*Q)@p1-M2@expm(t*Q)@p2))
                     for t in [0.,.1,1.,3.]]
    assert max(detector_errors)<1e-14
    # Legacy lift: same moment mechanism, no identity with strict assumed.
    WL,sL,RL,DL=legacy_cover()
    lawA=[(F(-1,2),F(1,2)),(F(1,2),F(1,2))]
    lawB=[(F(-1,4),F(4,5)),(F(1),F(1,5))]
    atomsA=[(RL+float(u)*DL,sL*float(p)) for u,p in lawA]
    atomsB=[(RL+float(u)*DL,sL*float(p)) for u,p in lawB]
    legacy_error=float(np.max(abs(intensity_generator(atomsA,2)-intensity_generator(atomsB,2))))
    legacy_gap=sL*float(F(3,16))*DL[0,1]**3
    assert legacy_error<1e-13 and legacy_gap>0
    ratio=F(3,5);reference=finite_recursive_law(ratio,8)
    approximations=[]
    for depth in range(1,7):
        mu=finite_recursive_law(ratio,depth);nextmu=finite_recursive_law(ratio,depth+1)
        residual=wasserstein_one(mu,nextmu);distance=wasserstein_one(mu,reference)
        residual_bound=residual/(1-ratio)
        assert distance<=residual_bound and distance<=ratio**depth
        var=sum(p*u*u for u,p in mu)
        assert F(1,4)-var==F(1,4)*ratio**(2*depth)
        approximations.append(dict(depth=depth,atoms=len(mu),recursion_residual=str(residual),
            finite_reference_W1=str(distance),residual_bound=str(residual_bound),
            tail_bound=str(ratio**depth),variance_bias=str(F(1,4)-var)))
    omega=.8;steps=3
    clock=omega*np.linalg.inv(omega*np.eye(12)-Q)
    invariances=[]
    for scale in [.2,3.,10.]:
        other=scale*omega*np.linalg.inv(scale*omega*np.eye(12)-scale*Q)
        invariances.append(float(np.max(abs(np.linalg.matrix_power(clock,steps)-
                                           np.linalg.matrix_power(other,steps)))))
    assert max(invariances)<1e-14
    return dict(completed_computational_rounds=list(range(21,30)),
        idle_clock_gauge=dict(original_total_rate=s,with_idle_total_rate=s+5,pair_generator_error=idle_error),
        erosion_limit=erosion,
        intensity_recovery_example=dict(points=[str(x) for x in points],masses=[str(x) for x in masses],
            weighted_masses=[str(x) for x in weighted],weighted_moments=[str(x) for x in moments],
            recovered_weighted_atoms=[(str(u),str(w)) for u,w in recovered],
            null_polynomial_norm=str(nullnorm),recovered_independent_rate=str(independent_rate)),
        quantum_finite_hierarchy=dict(harmonic=harmonic,epsilon=str(eps),checks=quantum,
            first_distinguishing_GHZ_trace_distance=str(eps/2)),
        preparation_detector_alias=dict(preparation_contrasts=[.6,.75],detector_contrasts=[.8,.64],
                                       output_errors=detector_errors),
        legacy_kernel_robustness=dict(cover_size=24,cover_rate=sL,pair_error=legacy_error,
            triple_transition_gap=legacy_gap,noise_amplitude=float(DL[0,1])),
        recursive_source_stability=approximations,
        internal_clock=dict(reference_rate=omega,ticks=steps,resolvent_row_sum_error=float(np.max(abs(clock.sum(axis=1)-1))),
            common_rescaling_errors=invariances))


if __name__=='__main__':
    output=run()
    Path(__file__).with_name('completion_results.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))
