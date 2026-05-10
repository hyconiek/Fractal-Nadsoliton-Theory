#!/usr/bin/env python3
"""P1172 safe operating region extraction for top candidate.

Computes robust region from P1171 perturbation grid and exports practical
parameter bounds where proxy remains stable.
"""
from __future__ import annotations
import json, math
from pathlib import Path

ROOT=Path(__file__).resolve().parent
GEN=ROOT/'generated'
ALPHA=4.0*math.log(2.0)


def stable(omega,phi,beta,eta,sigma,kappa,dmax=72.0,step=0.1):
    n=int(round(dmax/step))+1
    prev=0
    neg=0; sc=0
    for i in range(n):
        d=i*step
        num=math.cos(omega*d+phi)+sigma*(1.0-math.exp(-kappa*d))
        den=1.0+beta*(d**eta)
        v=math.exp(-ALPHA*d)*num/den
        if v<0: neg+=1
        s=1 if v>0 else (-1 if v<0 else 0)
        if s!=0 and prev!=0 and s!=prev: sc+=1
        if s!=0: prev=s
    return neg==0 and sc==0


def main():
    top=json.loads((GEN/'p1170_neighbor_1.json').read_text())
    o,p,b,e,s,k=top['omega_hint'],top['phi_hint'],top['beta_hint'],top['eta_hint'],top['sigma_hint'],top['kappa_hint']

    domega=[-0.02,-0.01,0.0,0.01,0.02]
    dphi=[-0.03,-0.015,0.0,0.015,0.03]
    dsigma=[-0.2,-0.1,0.0,0.1,0.2]
    dkappa=[-0.04,-0.02,0.0,0.02,0.04]

    points=[]
    for a in domega:
        for bb in dphi:
            for c in dsigma:
                for d in dkappa:
                    O,P,S,K=o+a,p+bb,s+c,k+d
                    ok=stable(O,P,b,e,S,K)
                    points.append({"omega":O,"phi":P,"sigma":S,"kappa":K,"stable":ok})

    stable_pts=[x for x in points if x['stable']]
    def bounds(key):
        vals=[x[key] for x in stable_pts]
        return {"min":min(vals),"max":max(vals)} if vals else None

    out={"packet":"P1172","as_of":"2026-05-10","grid_points":len(points),"stable_points":len(stable_pts),
         "stable_fraction":len(stable_pts)/len(points),
         "safe_bounds":{"omega":bounds('omega'),"phi":bounds('phi'),"sigma":bounds('sigma'),"kappa":bounds('kappa')},
         "note":"Safe operating region extraction only; no closure/QW-2191 discharge claim."}
    outp=GEN/'p1172_safe_operating_region_summary.json'
    outp.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(f"[P1172] stable={len(stable_pts)}/{len(points)} wrote {outp}")

if __name__=='__main__':
    main()
