#!/usr/bin/env python3
"""Ideal Gaussian sign-event enclosure, not a finite-N error certificate."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
SPECTRA=ROOT/"fin_rank7_followup/inputs/fin_handoff_audit/results.json"


def endpoint(t):
    s,m,e,b=t; assert b>=0
    x=F(int(m))*2**int(e) if e>=0 else F(int(m),1<<(-int(e)))
    return -x if s else x


def bounds(x):return [endpoint(t) for t in x._mpi_]


def public(x,digits=18):
    lo,hi=bounds(x); scale=10**digits
    a=lo.numerator*scale//lo.denominator; b=-((-hi.numerator*scale)//hi.denominator)
    def dec(n):
        sign="-" if n<0 else ""; n=abs(n)
        return sign+str(n//scale)+"."+str(n%scale).zfill(digits)
    result=[dec(a),dec(b)]; assert F(result[0])<=lo<=hi<=F(result[1]); return result


def build():
    mp.iv.dps=80; iv=mp.iv
    def rat(f):
        f=F(f); return iv.mpf(f.numerator)/f.denominator
    raw=json.loads(SPECTRA.read_text())["exact"]["laplacian_intervals"]
    lam={k:iv.mpf([rat(raw[k][0]).a,rat(raw[k][1]).b]) for k in (3,6)}
    gamma=[1-lam[3]/lam[6],1-(lam[3]+rat(F(1,100)))/lam[6]]
    def arcsin(x):
        assert 0<bounds(x)[0]<=bounds(x)[1]<F(9,10)
        term=x; total=iv.mpf(0); terms=110
        for n in range(terms):
            total+=term
            term*=x*x*rat(F((2*n+1)**2,2*(n+1)*(2*n+3)))
        # All remaining ratios are <=x^2<1, so this encloses the positive tail.
        return total+iv.mpf([0,1])*term/(1-x*x)
    q=iv.exp(-1)
    rho=[(iv.exp(-g)+g*q)/(1+g) for g in gamma]
    prob=[rat(F(1,2))+arcsin(r)/iv.pi for r in rho]
    gap=prob[1]-prob[0]; risk=(1-gap)/2
    assert bounds(gap)[0]>0 and bounds(risk)[1]<F(1,2)
    return {"date":"2026-10-15","status":"VALIDATED_IDEAL_GAUSSIAN_SIGN_EVENT_UNDER_INPUT_AND_IV_ASSUMPTIONS",
            "spectral_path":str(SPECTRA.relative_to(ROOT)),"spectral_sha256":hashlib.sha256(SPECTRA.read_bytes()).hexdigest(),
            "event":"V3(0)*V3(1)>0; alpha1/2, epsilon0,+1/100, native unit clock",
            "method":"rational spectral boxes; iv exp/pi;110 positive arcsin terms with geometric tail; exact dyadic outward decimal rounding",
            "rho_A7":public(rho[0]),"rho_control":public(rho[1]),"probability_A7":public(prob[0]),"probability_control":public(prob[1]),
            "Bernoulli_TV_gap":public(gap),"equal_prior_single_event_risk":public(risk),
            "not_certified":["finite-N accuracy or N0","physical sign detector/clock","Gaussian cost transfer to shrinking blocks or growing finite-N records"]}


if __name__=="__main__":
    parser=argparse.ArgumentParser(); parser.add_argument("--output",type=Path); args=parser.parse_args(); result=build()
    if args.output:args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
