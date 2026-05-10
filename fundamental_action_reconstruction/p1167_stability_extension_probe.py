#!/usr/bin/env python3
"""P1167 stability extension probe for P1166 local candidates.

Tests top P1166 candidate under:
1) extended domain [0,48]
2) small parameter perturbations
Reports robustness of sign_change_count=0 proxy.
"""
from __future__ import annotations
import json, math
from pathlib import Path

ROOT=Path(__file__).resolve().parent
GEN=ROOT/'generated'
ALPHA=4.0*math.log(2.0)


def metrics(omega,phi,beta,eta,sigma,kappa,dmax=24.0,step=0.1):
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
    s=json.loads((GEN/'p1166_sigma_kappa_sweep_summary.json').read_text())
    fp=s['fixed_params']; best=s['best']
    omega,phi,beta,eta=fp['omega'],fp['phi'],fp['beta'],fp['eta']
    sigma,kappa=best['sigma'],best['kappa']

    base24=metrics(omega,phi,beta,eta,sigma,kappa,24.0)
    ext48=metrics(omega,phi,beta,eta,sigma,kappa,48.0)

    deltas=[-0.01,0.0,0.01]
    perturb=[]
    for do in deltas:
        for dp in deltas:
            for ds in deltas:
                po,pn,psc=metrics(omega+do,phi+dp,beta,eta,sigma+ds,kappa,48.0)
                perturb.append({"domega":do,"dphi":dp,"dsigma":ds,
                                "positive_count":po,"negative_count":pn,
                                "sign_change_count":psc,
                                "stable_proxy": (psc==0 and pn==0)})
    robust=sum(1 for r in perturb if r['stable_proxy'])

    out={"packet":"P1167","as_of":"2026-05-10",
         "base_candidate":{"omega":omega,"phi":phi,"beta":beta,"eta":eta,"sigma":sigma,"kappa":kappa},
         "base24":{"positive_count":base24[0],"negative_count":base24[1],"sign_change_count":base24[2]},
         "extended48":{"positive_count":ext48[0],"negative_count":ext48[1],"sign_change_count":ext48[2]},
         "perturbation_cases":len(perturb),"robust_cases":robust,
         "robust_fraction": robust/len(perturb) if perturb else 0.0,
         "note":"Stability probe only; no closure/QW-2191 discharge claim."}
    outp=GEN/'p1167_stability_extension_probe_summary.json'
    outp.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(f"[P1167] robust={robust}/{len(perturb)} wrote {outp}")

if __name__=='__main__':
    main()
