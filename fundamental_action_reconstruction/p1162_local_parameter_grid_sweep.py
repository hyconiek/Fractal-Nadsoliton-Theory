#!/usr/bin/env python3
"""P1162 local parameter grid sweep around top P1161 candidate.

Searches local grid around candidate A and ranks by:
1) minimal sign_change_count
2) minimal negative_count
3) maximal positive_count
Methodological only.
"""
from __future__ import annotations
import json, math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
ALPHA = 4.0 * math.log(2.0)


def metrics(omega, phi, beta, eta):
    vals=[]
    for i in range(241):
        d=i*0.1
        vals.append(math.exp(-ALPHA*d)*math.cos(omega*d+phi)/(1.0+beta*(d**eta)))
    pos=sum(v>0 for v in vals); neg=sum(v<0 for v in vals)
    sc=0; prev=0
    for v in vals:
        s=1 if v>0 else (-1 if v<0 else 0)
        if s!=0 and prev!=0 and s!=prev: sc+=1
        if s!=0: prev=s
    return {"positive_count":pos,"negative_count":neg,"sign_change_count":sc,
            "qw_2191_proxy":"BLOCKED" if (neg>0 or sc>0) else "CANDIDATE_PASS_ONLY"}


def main():
    base=json.loads((GEN/'p1161_candidate_physical_a.json').read_text())
    bo, bp, bb, be = base['omega_hint'], base['phi_hint'], base['beta_hint'], base['eta_hint']

    omegas=[bo-0.01, bo, bo+0.01]
    phis=[bp-0.02, bp, bp+0.02]
    betas=[bb-0.05, bb, bb+0.05]
    etas=[be-0.05, be, be+0.05]

    rows=[]
    for o in omegas:
        for p in phis:
            for b in betas:
                for e in etas:
                    m=metrics(o,p,b,e)
                    rows.append({"omega":o,"phi":p,"beta":b,"eta":e,**m})

    rows.sort(key=lambda r:(r['sign_change_count'], r['negative_count'], -r['positive_count']))
    best=rows[0] if rows else None
    out={
      "packet":"P1162","as_of":"2026-05-10","base_candidate":"p1161_candidate_physical_a.json",
      "grid_size":len(rows),"best":best,"top10":rows[:10],
      "note":"Local sweep only; no closure/QW-2191 discharge claim."
    }
    outp=GEN/'p1162_local_parameter_grid_sweep_summary.json'
    outp.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(f"[P1162] grid={len(rows)} wrote {outp}")

if __name__=='__main__':
    main()
