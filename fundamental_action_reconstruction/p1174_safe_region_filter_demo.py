#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent
GEN=ROOT/'generated'
B=json.loads((GEN/'p1172_safe_operating_region_summary.json').read_text())['safe_bounds']

def check(c):
    return {
      'omega_in_bounds': B['omega']['min'] <= c.get('omega_hint',0) <= B['omega']['max'],
      'phi_in_bounds': B['phi']['min'] <= c.get('phi_hint',0) <= B['phi']['max'],
      'sigma_in_bounds': B['sigma']['min'] <= c.get('sigma_hint',0) <= B['sigma']['max'],
      'kappa_in_bounds': B['kappa']['min'] <= c.get('kappa_hint',0) <= B['kappa']['max'],
    }

def main():
    inside=GEN/'p1170_neighbor_1.json'
    outside=GEN/'p1174_candidate_outside_safe_region.json'
    rows=[]
    for p in [inside,outside]:
        c=json.loads(p.read_text())
        ch=check(c)
        rows.append({'candidate':str(p),'checks':ch,'safe_region_pass':all(ch.values())})
    out={'packet':'P1174','as_of':'2026-05-10','rows':rows,'note':'Safe-region prefilter demo only; no closure claim.'}
    outp=GEN/'p1174_safe_region_filter_demo_summary.json'
    outp.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(f'[P1174] wrote {outp}')

if __name__=='__main__':
    main()
