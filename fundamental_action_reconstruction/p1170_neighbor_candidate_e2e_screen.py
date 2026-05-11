#!/usr/bin/env python3
"""P1170: neighbor candidate E2E screen around P1166 top point."""
from __future__ import annotations
import json, math, subprocess, sys
from pathlib import Path

ROOT=Path(__file__).resolve().parent
GEN=ROOT/'generated'
ALPHA=4.0*math.log(2.0)


def proxy_metrics(omega,phi,beta,eta,sigma,kappa):
    vals=[]
    for i in range(241):
        d=i*0.1
        num=math.cos(omega*d+phi)+sigma*(1-math.exp(-kappa*d))
        den=1+beta*(d**eta)
        vals.append(math.exp(-ALPHA*d)*num/den)
    pos=sum(v>0 for v in vals); neg=sum(v<0 for v in vals)
    sc=0; prev=0
    for v in vals:
        s=1 if v>0 else (-1 if v<0 else 0)
        if s!=0 and prev!=0 and s!=prev: sc+=1
        if s!=0: prev=s
    return pos,neg,sc


def run_pipeline(candidate_path: Path):
    cmd=[sys.executable,str(ROOT/'p1151_strict_selector_pipeline_runner.py'),str(candidate_path)]
    p=subprocess.run(cmd,capture_output=True,text=True)
    summ=json.loads((GEN/'p1151_strict_selector_pipeline_runner_summary.json').read_text())
    return p.returncode, summ.get('overall_pass')


def main():
    src=json.loads((GEN/'p1166_sigma_kappa_sweep_summary.json').read_text())
    fp=src['fixed_params']; best=src['best']
    omega,phi,beta,eta=fp['omega'],fp['phi'],fp['beta'],fp['eta']
    s0,k0=best['sigma'],best['kappa']
    neighbors=[(s0,k0),(s0-0.1,k0),(s0,k0+0.04),(s0-0.1,k0+0.04),(s0-0.2,k0+0.08)]

    rows=[]
    for i,(s,k) in enumerate(neighbors, start=1):
        c={"premise_id":f"P1170_NEIGHBOR_{i}","strict_side_provenance":True,
           "noncyclic_anchor_declared":True,"no_legacy_role_transfer":True,
           "no_closure_claim":True,"no_qw2191_discharge_claim":True,
           "omega_hint":omega,"phi_hint":phi,"beta_hint":beta,"eta_hint":eta,
           "sigma_hint":s,"kappa_hint":k}
        cpath=GEN/f"p1170_neighbor_{i}.json"
        cpath.write_text(json.dumps(c,indent=2,sort_keys=True)+"\n")
        pos,neg,sc=proxy_metrics(omega,phi,beta,eta,s,k)
        rc,pass_ok=run_pipeline(cpath)
        rows.append({"candidate":str(cpath),"sigma":s,"kappa":k,"positive_count":pos,
                     "negative_count":neg,"sign_change_count":sc,
                     "proxy_status":"CANDIDATE_PASS_ONLY" if (sc==0 and neg==0) else "BLOCKED",
                     "pipeline_pass":bool(pass_ok),"pipeline_returncode":rc})

    rows.sort(key=lambda r:(-int(r['pipeline_pass']), r['sign_change_count'], r['negative_count']))
    out={"packet":"P1170","as_of":"2026-05-10","rows":rows,
         "top_recommendation":rows[0]['candidate'] if rows else None,
         "note":"Neighbor E2E screen only; no closure/QW-2191 discharge claim."}
    outp=GEN/'p1170_neighbor_candidate_e2e_screen_summary.json'
    outp.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(f"[P1170] screened {len(rows)} candidates wrote {outp}")

if __name__=='__main__':
    main()
