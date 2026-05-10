#!/usr/bin/env python3
"""P1166 2D sweep for asymmetric selector term class.

Scans (sigma, kappa) for P1165 class at fixed base params and reports whether
any local combination achieves sign_change_count=0.
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
    sigmas=[0.2,0.4,0.6,0.8,1.0,1.2]
    kappas=[0.04,0.08,0.12,0.16,0.2,0.3]
    rows=[]
    for s in sigmas:
        for k in kappas:
            pos,neg,sc=eval_metrics(omega,phi,beta,eta,s,k)
            rows.append({"sigma":s,"kappa":k,"positive_count":pos,"negative_count":neg,
                         "sign_change_count":sc,
                         "qw_2191_proxy":"BLOCKED" if (neg>0 or sc>0) else "CANDIDATE_PASS_ONLY"})
    rows.sort(key=lambda r:(r['sign_change_count'], r['negative_count'], -r['positive_count']))
    zero=[r for r in rows if r['sign_change_count']==0]
    out={"packet":"P1166","as_of":"2026-05-10",
         "fixed_params":{"omega":omega,"phi":phi,"beta":beta,"eta":eta},
         "grid_size":len(rows),"zero_sign_change_count":len(zero),
         "best":rows[0],"top10":rows[:10],
         "note":"2D sigma-kappa sweep only; no closure/QW-2191 discharge claim."}
    outp=GEN/'p1166_sigma_kappa_sweep_summary.json'
    outp.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(f"[P1166] grid={len(rows)} zero_sc={len(zero)} wrote {outp}")

if __name__=='__main__':
    main()
