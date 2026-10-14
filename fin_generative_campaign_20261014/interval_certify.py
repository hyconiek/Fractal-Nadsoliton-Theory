#!/usr/bin/env python3
"""Outward mpmath.iv source enclosure + EXACT rational ANOVA of nominal logs.
Certified result is conditional on mpmath interval transcendentals correctly enclosing the exact decimal kernel.
"""
from __future__ import annotations
import itertools, json, hashlib, time, math
from pathlib import Path
from fractions import Fraction
import numpy as np
import mpmath as mp
from run_campaign import build_rank7, marginal_F3, Q, G, ROOT

def frac_mpf_interval_endpoint(z):
    # mpmath interval stores exact binary rationals (sign, mantissa, binary exp, bitcount)
    sign, man, exp, bc=z
    numerator= -man if sign else man
    return Fraction(numerator*(2**exp),1) if exp>=0 else Fraction(numerator,2**(-exp))

def interval_fin():
    mp.iv.dps=52
    I=mp.iv
    omega=I.mpf('0.18575'); phi=I.mpf('0.1625'); eta=I.mpf('1.8')
    gg=I.mpf('5.145228719489142'); ivpi=I.pi
    W=[]
    for d in range(12):
        if d==0: W.append(I.mpf(0))
        else:
            dist=min(d,12-d)
            W.append(I.cos(omega*dist+phi)/(1+I.power(I.mpf(dist),eta)))
    lam={k:sum([W[d]*(1-I.cos(2*ivpi*k*d/12)) for d in range(1,Q)],I.mpf(0)) for k in (3,4,5,6)}
    A=[]
    for d in range(Q):
        A.append(sum([lam[k]/6*I.cos(2*ivpi*k*d/Q) for k in (3,4,5)],I.mpf(0))+(lam[6]/12)*((-1)**d))
    return gg,A,lam

def exact_pair_component_nominal(F):
    # For floating F treated as exact dyadic rationals: third ANOVA component
    fr=np.empty((12,12,12),dtype=object)
    for idx in np.ndindex(fr.shape): fr[idx]=Fraction.from_float(float(F[idx]))
    m0=np.empty((12,12),dtype=object);m1=m0.copy();m2=m0.copy()
    for a in range(12):
        for b in range(12):
            m0[a,b]=sum((fr[k,a,b] for k in range(12)),Fraction())/12
            m1[a,b]=sum((fr[a,k,b] for k in range(12)),Fraction())/12
            m2[a,b]=sum((fr[a,b,k] for k in range(12)),Fraction())/12
    # Exact rational means over pairs of active coordinates
    m01=[sum((fr[a,b,k] for a in range(12) for b in range(12)),Fraction())/144 for k in range(12)]
    m02=[sum((fr[a,k,c] for a in range(12) for c in range(12)),Fraction())/144 for k in range(12)]
    m12=[sum((fr[k,b,c] for b in range(12) for c in range(12)),Fraction())/144 for k in range(12)]
    globalmean=sum(fr.flat,Fraction())/1728
    vals=[]
    for a,b,c in np.ndindex(fr.shape):
        vals.append(fr[a,b,c]-m0[b,c]-m1[a,c]-m2[a,b]+m01[c]+m02[b]+m12[a]-globalmean)
    return max(vals)-min(vals),vals

def main():
    t=time.time()
    _,_,_,A_double=build_rank7()
    F_double=marginal_F3(A_double,G)
    gg,A,lam=interval_fin()
    # A7 is circulant, reduce every many-spin evaluation to 12 diagonal differences.
    E=Fraction()
    h0=Fraction()
    for i,j in itertools.product(range(Q),repeat=2):
        iv=A[(j-i)%Q]
        for endpoint in iv._mpi_:
            E=max(E,abs(frac_mpf_interval_endpoint(endpoint)-Fraction.from_float(float(A_double[i,j]))))
    # Tight interval for each exact effective 3-label logweight and its difference from binary64 producer F.
    maxdiff=Fraction()
    max_width=Fraction()
    for a,b,c in itertools.product(range(Q),repeat=3):
        states=(a,b,c)
        # sum over active 3x3 entries and logsum 12 terms
        pair=sum((A[(z-x)%Q] for x in states for z in states), mp.iv.mpf(0))*(gg/8)
        vals=[]
        for j in range(Q):
            vals.append(sum((A[(x-j)%Q] for x in states),mp.iv.mpf(0))*(gg/4))
        val=pair+mp.iv.ln(sum((mp.iv.exp(v) for v in vals),mp.iv.mpf(0)))
        fd=Fraction.from_float(float(F_double[a,b,c]))
        low=frac_mpf_interval_endpoint(val._mpi_[0]); high=frac_mpf_interval_endpoint(val._mpi_[1]);
        assert low<=high
        maxdiff=max(maxdiff, abs(low-fd),abs(high-fd));max_width=max(max_width,high-low)
    omega_float_rational,values=exact_pair_component_nominal(F_double)
    # Three active coordinates; each (I-E_i) has sup operator norm <=2.
    # r_exact=(I-E1)(I-E2)(I-E3)F_exact; ||r_exact-r_nom||_inf <= 8 epsF;
    # osc difference <= 16 epsF, with epsF from interval endpoints including the actual binary64 float.
    delta=16*maxdiff
    lo=max(Fraction(0),omega_float_rational-delta)
    hi=omega_float_rational+delta
    # Global TV bound monotonic tanh; calculate outward using mpmath.iv from decimal rational enclosures
    mp.iv.dps=50
    def to_iv(fr):
        return mp.iv.mpf(fr.numerator)/mp.iv.mpf(fr.denominator)
    # interval rational num/den conversion is outward-enclosed.
    e=mp.iv.exp(to_iv(hi)/2)
    B=(e-1)/(e+1)
    def dec(z):
        mp.mp.dps=65
        return mp.nstr(mp.mpf(z.numerator)/z.denominator, 45)
    res={'execution_status':'EXECUTED','scientific_verdict':'VALIDATED_INTERVAL_ENCLOSURE_CONDITIONAL_ON_MPMATH_IV_BACKEND','mathematical_input':'exact decimal omega=0.18575 phi=0.1625 eta=1.8 g=5.145228719489142, exact pi','source_certifier':'independent circulant Fourier A7 calculation with outward interval arithmetic','source_A7_max_discrepancy_from_binary64':dec(E),'F_exact_max_discrepancy_from_binary64_upper':dec(maxdiff),'F_exact_max_interval_width':dec(max_width),'nominal_exact_dyadic_ANOVA_osc':dec(omega_float_rational),'validated_residual_osc_lower':dec(lo),'validated_residual_osc_upper':dec(hi),'validated_TV_bound_upper':mp.nstr(mp.make_mpf(B._mpi_[1]),40),'third_contrast_finite_numerical':None,'method':'1728 outward interval logweights + rational (binary64-as-exact) ANOVA + projection norm bound ||Pi3||∞≤8; no uncertified floating argmin used','runtime_seconds':time.time()-t,'limitations':['mpmath.iv interval arithmetic implementation is assumed sound, not formally verified','certifies the exact decimal mathematical tuple, not a separate implementation whose gain is interpreted as an IEEE binary64 constant','physical source for kernel is outside this proof']}
    p=ROOT/'RG-C02/INTERVAL_CERTIFICATE.json';p.write_text(json.dumps(res,indent=2)+'\n')
    (ROOT/'RG-C02/INTERVAL_PROOF.md').write_text('# Exact-kernel numerical certificate\n\nThe input values are **exact decimals**; `pi` is the mathematical circle constant enclosed by mpmath.iv. The k=3,4,5,6 Fourier coefficients of the finite q12 kernel are evaluated with interval functions, then A7 is built by the exact circulant Fourier sum. All 1728 logweights are outward evaluated. The maximal exact discrepancy between each interval endpoint and the binary64 producer value (regarded as an exact dyadic rational) is ε. Computing the third Hoeffding component using **exact rational arithmetic** on those dyadics yields ω0. Since ||(I-E1)(I-E2)(I-E3)||∞ ≤ 8, true residual oscillation belongs to [ω0−16ε,ω0+16ε]. The ensuing tanh(ω/4) upper bound uses outward interval tanh. This is a mathematically validated bound conditional on the mpmath.iv implementation (not a formal Lean/Coq proof).\n\n'+json.dumps(res,indent=2)+'\n')
    print(json.dumps(res,indent=2))
if __name__=='__main__':main()
