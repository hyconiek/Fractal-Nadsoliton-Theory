#!/usr/bin/env python3
"""P1164 non-phase damping class probe.

Tests alternative positive damping class at fixed phase:
K_gamma(d)=cos(omega*d+phi)/(1+beta*d^eta+gamma*d^(eta+1)), gamma>=0.
Checks whether sign_change_count can reach 0 on [0,24].
"""
from __future__ import annotations
import json, math
from pathlib import Path

ROOT=Path(__file__).resolve().parent
GEN=ROOT/'generated'
ALPHA=4.0*math.log(2.0)


def eval_metrics(omega,phi,beta,eta,gamma):
    vals=[]
    for i in range(241):
        d=i*0.1
        den=1.0+beta*(d**eta)+gamma*(d**(eta+1.0))
        vals.append(math.exp(-ALPHA*d)*math.cos(omega*d+phi)/den)
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
    gammas=[0.0,0.05,0.1,0.2,0.4,0.8,1.2,1.6,2.0]
    rows=[]
    for g in gammas:
        pos,neg,sc=eval_metrics(omega,phi,beta,eta,g)
        rows.append({"gamma":g,"positive_count":pos,"negative_count":neg,"sign_change_count":sc,
                     "qw_2191_proxy":"BLOCKED" if (neg>0 or sc>0) else "CANDIDATE_PASS_ONLY"})
    zero_sc=[r for r in rows if r['sign_change_count']==0]
    out={"packet":"P1164","as_of":"2026-05-10",
         "fixed_params":{"omega":omega,"phi":phi,"beta":beta,"eta":eta},
         "rows":rows,"zero_sign_change_count":len(zero_sc),
         "no_go_local_slice": len(zero_sc)==0,
         "note":"Alternative damping-class test only; no closure/QW-2191 discharge claim."}
    outp=GEN/'p1164_nonphase_damping_class_probe_summary.json'
    outp.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(f"[P1164] gamma_cases={len(rows)} zero_sc={len(zero_sc)} wrote {outp}")

if __name__=='__main__':
    main()
