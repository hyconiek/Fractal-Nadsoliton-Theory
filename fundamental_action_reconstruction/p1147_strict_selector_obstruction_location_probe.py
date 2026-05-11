#!/usr/bin/env python3
"""P1147: locate first sign-flip obstruction for strict selector candidate family.

Family:
S_gamma(d) = exp(-alpha*d) * exp(-gamma*d^2) * cos(omega*d+phi)/(1+beta*d^eta), gamma>=0.

Because positive weights cannot change cosine zeros, first sign-flip location is
weight-invariant; this is exported as a nonclosure obstruction certificate proxy.
"""
from __future__ import annotations
import json, math
from pathlib import Path

ALPHA=4.0*math.log(2.0)
OMEGA=0.18575
PHI=0.16250
BETA=1.0
ETA=1.8


def s_gamma(d: float, gamma: float) -> float:
    return math.exp(-ALPHA*d-gamma*d*d)*math.cos(OMEGA*d+PHI)/(1.0+BETA*d**ETA)


def first_sign_flip(gamma: float, d_max: float=40.0, step: float=1e-3):
    prev = s_gamma(0.0, gamma)
    d=step
    while d<=d_max:
        cur = s_gamma(d,gamma)
        if prev==0:
            return d-step,d-step,d-step
        if cur==0 or prev*cur<0:
            a,b=d-step,d
            for _ in range(60):
                m=(a+b)/2
                fm=s_gamma(m,gamma)
                fa=s_gamma(a,gamma)
                if fa==0:
                    b=a;break
                if fm==0:
                    a=b=m;break
                if fa*fm<=0: b=m
                else: a=m
            r=(a+b)/2
            return a,b,r
        prev=cur
        d+=step
    return None


def main():
    gammas=[0.0,0.25,0.5,1.0,2.0,4.0]
    roots={}
    for g in gammas:
        res=first_sign_flip(g)
        roots[str(g)]={"bracket": [res[0],res[1]], "root_estimate": res[2]} if res else None

    expected=(math.pi/2-PHI)/OMEGA
    out={
        "packet":"P1147",
        "as_of":"2026-05-10",
        "family":"S_gamma(d)=exp(-alpha*d-gamma*d^2)*cos(omega*d+phi)/(1+beta*d^eta)",
        "params":{"alpha":ALPHA,"omega":OMEGA,"phi":PHI,"beta":BETA,"eta":ETA},
        "analytic_first_zero_from_phase":expected,
        "numeric_first_sign_flip_by_gamma":roots,
        "audit":{
            "weight_invariant_flip": True,
            "qw_2191_nonclosure_proxy_verdict":"BLOCKED",
            "reason":"positive monotone weights cannot remove cosine-induced first sign flip"
        }
    }
    p=Path(__file__).resolve().parent/"generated"/"p1147_strict_selector_obstruction_location_probe_summary.json"
    p.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(f"[P1147] wrote {p}")

if __name__=='__main__':
    main()
