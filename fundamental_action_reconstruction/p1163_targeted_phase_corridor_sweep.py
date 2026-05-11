#!/usr/bin/env python3
"""P1163 targeted phase corridor sweep.

Fix beta/eta at P1162-best and sweep omega/phi corridor to test whether
sign_change_count can reach zero on d in [0,24].
"""
from __future__ import annotations
import json, math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
ALPHA = 4.0 * math.log(2.0)


def eval_metrics(omega, phi, beta, eta):
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
    return pos,neg,sc


def main():
    base = json.loads((GEN/'p1162_local_parameter_grid_sweep_summary.json').read_text())['best']
    beta, eta = base['beta'], base['eta']
    o0, p0 = base['omega'], base['phi']

    omegas=[o0 + k*0.005 for k in range(-4,5)]
    phis=[p0 + k*0.01 for k in range(-4,5)]

    rows=[]
    for o in omegas:
        for p in phis:
            pos,neg,sc = eval_metrics(o,p,beta,eta)
            rows.append({"omega":o,"phi":p,"beta":beta,"eta":eta,
                         "positive_count":pos,"negative_count":neg,
                         "sign_change_count":sc,
                         "qw_2191_proxy":"BLOCKED" if (neg>0 or sc>0) else "CANDIDATE_PASS_ONLY"})

    rows.sort(key=lambda r:(r['sign_change_count'], r['negative_count'], -r['positive_count']))
    zero_sc = [r for r in rows if r['sign_change_count']==0]
    out={
      "packet":"P1163","as_of":"2026-05-10",
      "fixed_beta_eta":{"beta":beta,"eta":eta},
      "center":{"omega":o0,"phi":p0},
      "grid_size":len(rows),
      "zero_sign_change_count":len(zero_sc),
      "best":rows[0],
      "top10":rows[:10],
      "note":"Corridor sweep only; no closure/QW-2191 discharge claim."
    }
    outp=GEN/'p1163_targeted_phase_corridor_sweep_summary.json'
    outp.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(f"[P1163] grid={len(rows)} zero_sc={len(zero_sc)} wrote {outp}")

if __name__=='__main__':
    main()
