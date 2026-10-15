#!/usr/bin/env python3
"""Positive-series Gaussian event enclosures from inherited rational spectra.

Assumes correctness of mpmath.iv. Public endpoints are rounded outward using
exact Fractions of raw dyadic endpoints, not nearest decimal serialization.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
SPECTRA=ROOT/"fin_rank7_followup/inputs/fin_handoff_audit/results.json"


def dyadic(t):
    sign,man,exponent,bits=t
    assert bits>=0
    value=F(int(man))*2**int(exponent) if exponent>=0 else F(int(man),1<<(-int(exponent)))
    return -value if sign else value


def endpoints(x):return [dyadic(t) for t in x._mpi_]


def public(x,digits=18):
    lo,hi=endpoints(x); scale=10**digits
    a=lo.numerator*scale//lo.denominator
    b=-((-hi.numerator*scale)//hi.denominator)
    def decimal(n):
        prefix="-" if n<0 else ""; n=abs(n)
        return prefix+str(n//scale)+"."+str(n%scale).zfill(digits)
    result=[decimal(a),decimal(b)]
    assert F(result[0])<=lo<=hi<=F(result[1])
    return result


def build():
    mp.iv.dps=80; iv=mp.iv
    def rational(f):
        f=F(f); return iv.mpf(f.numerator)/f.denominator
    boxes=json.loads(SPECTRA.read_text())["exact"]["laplacian_intervals"]
    lam={k:iv.mpf([rational(boxes[k][0]).a,rational(boxes[k][1]).b]) for k in (3,4,5,6)}
    # Independent full first-row DFT, not the producer's paired-distance sum.
    w=[iv.mpf(0)]
    for j in range(1,12):
        d=min(j,12-j)
        w.append(iv.cos(rational(F(743,4000))*d+rational(F(13,80)))/(1+iv.mpf(d)**rational(F(9,5))))
    first=[sum(w)]+[-w[j] for j in range(1,12)]
    direct={k:sum(first[j]*iv.cos(2*iv.pi*k*j/12) for j in range(12)) for k in (3,4,5,6)}
    for k in direct:
        lo,hi=endpoints(direct[k]); assert F(boxes[k][0])<=lo<=hi<=F(boxes[k][1]),k
    # Positive identity: erf x=2/sqrt(pi)*exp(-x^2)*sum 2^n*x^(2n+1)/(2n+1)!!.
    def event(gamma):
        variance=1+1/gamma; x=1/iv.sqrt(2*variance)
        assert endpoints(x)[0]>0 and endpoints(x)[1]<F(1,3)
        term=x; total=iv.mpf(0); terms=18
        for n in range(terms):
            total+=term
            term*=2*x*x/(2*n+3)
        ratio=2*x*x/(2*terms+3)
        assert endpoints(ratio)[1]<1
        tail=term/(1-ratio)
        erf=2/iv.sqrt(iv.pi)*iv.exp(-x*x)*(total+iv.mpf([0,1])*tail)
        return 1-erf
    g0=1-lam[3]/lam[6]; ge=1-(lam[3]+rational(F(1,100)))/lam[6]
    p0,pe=event(g0),event(ge); gap=pe-p0; risk=(1-gap)/2
    assert endpoints(gap)[0]>0 and endpoints(risk)[1]<F(1,2)
    # Recheck the comparator transverse maximum and old exact geometry payment.
    assert F(boxes[3][1])+F(1,100)<F(boxes[5][0])
    assert F(boxes[4][1])+F(1,100)<F(boxes[5][0])
    gamma5=1-F(boxes[5][1])/F(boxes[6][0])
    assert gamma5>F(93,5000)
    pmax=F(101,1200)+F(1,10000)
    lambda0_floor=1/pmax-12*(1-F(93,5000))
    hessian_floor=lambda0_floor-F(1,1000)*F(12,5)  # retain the older gain-box reserve
    assert hessian_floor==F(8679,98750)>F(2,25)
    return {"date":"2026-10-15","status":"VALIDATED_IDEAL_LIMIT_EVENT_UNDER_INHERITED_INPUT_AND_IV_BACKEND",
            "spectral_source":str(SPECTRA.relative_to(ROOT)),"spectral_sha256":hashlib.sha256(SPECTRA.read_bytes()).hexdigest(),
            "method":"accepted rational spectral boxes; fresh interval DFT input binding;18 positive erf-series terms and geometric tail; exact dyadic-to-decimal outward rounding",
            "probability_A7":public(p0),"probability_comparator":public(pe),"Bernoulli_TV_gap":public(gap),"equal_prior_single_event_risk":public(risk),
            "gamma3_A7":public(g0),"gamma3_comparator":public(ge),"same_transverse_maximum":"Lambda5",
            "gamma5_lower_exact":str(gamma5),"fast_Hessian_floor_exact":str(hessian_floor),"lambda0_Hessian_floor_exact":str(lambda0_floor),
            "not_certified":["finite-N error or usable N0","physical readout calibration","all65 producer decimal endpoint digits","fast-time process theorem"]}


if __name__=="__main__":
    parser=argparse.ArgumentParser(); parser.add_argument("--output",type=Path); args=parser.parse_args(); result=build()
    if args.output:args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
