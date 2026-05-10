#!/usr/bin/env python3
"""P1165 asymmetric selector-term probe.

Tests a strict-side asymmetric additive numerator term:
N_sigma(d)=cos(omega*d+phi) + sigma*(1-exp(-kappa*d)), sigma>=0.
K_sigma(d)=exp(-alpha*d)*N_sigma(d)/(1+beta*d^eta).
Evaluates whether local sign-change burden can be removed.
"""
from __future__ import annotations
import json, math
from pathlib import Path

ROOT=Path(__file__).resolve().parent
GEN=ROOT/'generated'
ALPHA=4.0*math.log(2.0)


def eval_metrics(omega,phi,beta,eta,sigma,kappa):
    vals=[]
    for i in range(241):
        d=i*0.1
        num=math.cos(omega*d+phi)+sigma*(1.0-math.exp(-kappa*d))
        den=1.0+beta*(d**eta)
        vals.append(math.exp(-ALPHA*d)*num/den)
    pos=sum(v>0 for v in vals); neg=sum(v<0 for v in vals)
    sc=0; prev=0
    for v in vals:
        s=1 if v>0 else (-1 if v<0 else 0)
        if s!=0 and prev!=0 and s!=prev: sc+=1
        if s!=0: prev=s
    return pos,neg,sc


def main():
    b=json.loads((GEN/'p1163_targeted_phase_corridor_sweep_summary.json').read_text())['best']
    omega,phi,beta,eta=b['omega'],b['phi'],b['beta'],b['eta']
    kappa=0.12
    sigmas=[0.0,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0]
    rows=[]
    for s in sigmas:
        pos,neg,sc=eval_metrics(omega,phi,beta,eta,s,kappa)
        rows.append({"sigma":s,"kappa":kappa,"positive_count":pos,"negative_count":neg,
                     "sign_change_count":sc,
                     "qw_2191_proxy":"BLOCKED" if (neg>0 or sc>0) else "CANDIDATE_PASS_ONLY"})
    rows.sort(key=lambda r:(r['sign_change_count'], r['negative_count'], -r['positive_count']))
    out={"packet":"P1165","as_of":"2026-05-10",
         "fixed_params":{"omega":omega,"phi":phi,"beta":beta,"eta":eta},
         "rows":rows,"best":rows[0],
         "zero_sign_change_count":len([r for r in rows if r['sign_change_count']==0]),
         "note":"Asymmetric-term candidate probe only; no closure/QW-2191 discharge claim."}
    outp=GEN/'p1165_asymmetric_selector_term_probe_summary.json'
    outp.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(f"[P1165] sigma_cases={len(rows)} zero_sc={out['zero_sign_change_count']} wrote {outp}")

if __name__=='__main__':
    main()
