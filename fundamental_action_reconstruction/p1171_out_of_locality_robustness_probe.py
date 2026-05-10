#!/usr/bin/env python3
"""P1171 out-of-locality robustness probe for top P1170 candidate.

Checks proxy stability on extended domain and broader perturbation ranges in
(omega,phi,sigma,kappa).
"""
from __future__ import annotations
import json, math, sys
from pathlib import Path

ROOT=Path(__file__).resolve().parent
GEN=ROOT/'generated'
ALPHA=4.0*math.log(2.0)


def metrics(omega,phi,beta,eta,sigma,kappa,dmax=72.0,step=0.1):
    n=int(round(dmax/step))+1
    vals=[]
    for i in range(n):
        d=i*step
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
    args = sys.argv[1:]
    candidate_path = Path(args[0]).resolve() if args else (GEN/'p1170_neighbor_1.json').resolve()
    top = json.loads(candidate_path.read_text())
    omega,phi,beta,eta = top['omega_hint'],top['phi_hint'],top['beta_hint'],top['eta_hint']
    sigma,kappa = top['sigma_hint'],top['kappa_hint']

    base72=metrics(omega,phi,beta,eta,sigma,kappa,72.0)

    d_omega=[-0.02,0.0,0.02]
    d_phi=[-0.03,0.0,0.03]
    d_sigma=[-0.2,0.0,0.2]
    d_kappa=[-0.04,0.0,0.04]

    cases=[]
    for a in d_omega:
        for b in d_phi:
            for c in d_sigma:
                for d in d_kappa:
                    pos,neg,sc=metrics(omega+a,phi+b,beta,eta,sigma+c,kappa+d,72.0)
                    cases.append({"domega":a,"dphi":b,"dsigma":c,"dkappa":d,
                                  "positive_count":pos,"negative_count":neg,
                                  "sign_change_count":sc,
                                  "stable_proxy":(neg==0 and sc==0)})
    robust=sum(1 for r in cases if r['stable_proxy'])
    out={"packet":"P1171","as_of":"2026-05-10",
         "candidate_path": str(candidate_path),
         "base_candidate":{"omega":omega,"phi":phi,"beta":beta,"eta":eta,"sigma":sigma,"kappa":kappa},
         "base72":{"positive_count":base72[0],"negative_count":base72[1],"sign_change_count":base72[2]},
         "cases":len(cases),"robust_cases":robust,"robust_fraction":robust/len(cases),
         "note":"Out-of-locality robustness probe only; no closure/QW-2191 discharge claim."}
    outp=GEN/'p1171_out_of_locality_robustness_probe_summary.json'
    outp.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(f"[P1171] robust={robust}/{len(cases)} wrote {outp}")

if __name__=='__main__':
    main()
