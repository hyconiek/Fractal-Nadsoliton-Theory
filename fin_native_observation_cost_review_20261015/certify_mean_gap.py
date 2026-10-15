#!/usr/bin/env python3
"""Outward ideal mean-gap and explicit corrector-tail arithmetic, not native N0."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
SPECTRA=ROOT/"fin_rank7_followup/inputs/fin_handoff_audit/results.json"
def dyadic(t):
    s,m,e,b=t; assert b>=0
    value=F(int(m))*2**int(e) if e>=0 else F(int(m),1<<(-int(e)))
    return -value if s else value
def bounds(x):return [dyadic(t) for t in x._mpi_]
def public(x,digits=18):
    lo,hi=bounds(x); scale=10**digits
    a=lo.numerator*scale//lo.denominator; b=-((-hi.numerator*scale)//hi.denominator)
    def dec(n):
        prefix="-" if n<0 else ""; n=abs(n)
        return prefix+str(n//scale)+"."+str(n%scale).zfill(digits)
    out=[dec(a),dec(b)]; assert F(out[0])<=lo<=hi<=F(out[1]); return out

def build():
    mp.iv.dps=80; iv=mp.iv
    def rat(value):
        f=F(value); return iv.mpf(f.numerator)/f.denominator
    raw=json.loads(SPECTRA.read_text())["exact"]["laplacian_intervals"]
    lam={k:iv.mpf([rat(raw[k][0]).a,rat(raw[k][1]).b]) for k in (3,6)}
    gs=[1-lam[3]/lam[6],1-(lam[3]+rat(F(1,100)))/lam[6]]
    d=iv.exp(-rat(F(1,2)))*(1/gs[1]-1/gs[0])/2
    assert bounds(d)[0]>0
    # MVT: gap/alpha >=d*exp(-alpha*beta_control).
    beta=(1/gs[1]-1)/2; reserve=iv.exp(-beta/4)
    assert bounds(reserve)[0]>F(1,2)
    tails=[]
    for g in gs:
        for alpha in (rat(F(1,100)),rat(F(1,2))):
            K=240
            bound=(alpha/(2*g)*iv.exp(-2*g*K)/(1-iv.exp(-2*g))
                   +(1-alpha)/2*iv.exp(-2*K)/(1-iv.exp(-2))
                   +iv.sqrt(alpha)*rat(F(6,5))*iv.exp(-g*K)/(1-iv.exp(-g))
                   +rat(F(7,20))*iv.exp(-K)/(1-iv.exp(-1)))
            assert bounds(bound)[1]<F(5,10**15)
            tails.append(public(bound,digits=24))
    return {"date":"2026-10-15","status":"VALIDATED_IDEAL_GAP_AND_POINT_TAILS_UNDER_INPUT_AND_IV_ASSUMPTIONS",
            "spectral_path":str(SPECTRA.relative_to(ROOT)),"spectral_sha256":hashlib.sha256(SPECTRA.read_bytes()).hexdigest(),
            "positive_linear_mean_gap_coefficient":public(d),"uniform_alpha_interval":["0 (excluded)","1/4"],
            "uniform_ideal_mean_gap_bound":"mu_A7(alpha)-mu_control(alpha)>= (d/2)*alpha for0<alpha<=1/4",
            "positive_lower_coefficient":public(d/2),"tail_fixture":"sum from ell240 onward, |f|=6/5,|r|=7/20,alpha=.01,.5,each frozen gamma",
            "corrector_tail_upper_enclosures":tails,
            "rounding":"exact raw-dyadic Fraction floor/ceiling to public decimals",
            "not_certified":["all source70-digit decimals","native constants or useful N0","finite-N probability at an illustrative N","physical observation/calibration"]}

if __name__=="__main__":
    parser=argparse.ArgumentParser(); parser.add_argument("--output",type=Path); args=parser.parse_args(); result=build()
    if args.output:args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
