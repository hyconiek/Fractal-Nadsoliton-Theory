#!/usr/bin/env python3
"""P1149: reproducibility audit for P1146-P1148 generated summaries.

Recomputes core observables from probe definitions and checks equality/tolerance
against committed JSON artifacts. Produces one machine-readable audit verdict.
"""
from __future__ import annotations
import json, math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

OMEGA=0.18575
PHI=0.1625
BETA=1.0
ETA=1.8
ALPHA=4.0*math.log(2.0)


def load(name:str):
    return json.loads((GEN/name).read_text(encoding='utf-8'))


def approx(a,b,tol=1e-9):
    return abs(a-b)<=tol


def audit_p1146():
    j=load('p1146_strict_shannon_selector_candidate_grid_nonclosure_probe_summary.json')
    vals=[]
    for i in range(j['domain']['points']):
        d=j['domain']['d_min']+i*j['domain']['step']
        v=math.exp(-ALPHA*d)*math.cos(OMEGA*d+PHI)/(1+BETA*d**ETA)
        vals.append(v)
    pos=sum(v>0 for v in vals); neg=sum(v<0 for v in vals); zero=len(vals)-pos-neg
    return {
      'positive_count_ok': pos==j['observables']['positive_count'],
      'negative_count_ok': neg==j['observables']['negative_count'],
      'zero_count_ok': zero==j['observables']['zero_count'],
      'max_ok': approx(max(vals),j['observables']['max_s_cand'],1e-12),
      'min_ok': approx(min(vals),j['observables']['min_s_cand'],1e-12),
      'verdict_ok': j['audit']['qw_2191_nonclosure_verdict']=='BLOCKED'
    }


def audit_p1147():
    j=load('p1147_strict_selector_obstruction_location_probe_summary.json')
    expected=(math.pi/2-PHI)/OMEGA
    roots=j['numeric_first_sign_flip_by_gamma']
    ok_analytic=approx(j['analytic_first_zero_from_phase'],expected,1e-12)
    ok_roots=all(approx(v['root_estimate'],expected,1e-9) for v in roots.values())
    return {
      'analytic_ok': ok_analytic,
      'roots_ok': ok_roots,
      'invariant_flag_ok': bool(j['audit']['weight_invariant_flip']) is True,
      'verdict_ok': j['audit']['qw_2191_nonclosure_proxy_verdict']=='BLOCKED'
    }


def audit_p1148():
    j=load('p1148_strict_phase_shifted_selector_family_probe_summary.json')
    threshold=(math.pi/2)-PHI-OMEGA*j['domain']['d_max']
    ok_threshold=approx(j['delta_threshold_for_no_flip_on_domain'],threshold,1e-12)
    ok_empty=(len(j['admissible_no_flip_deltas'])==0)
    ok_verdict=(j['audit']['verdict']=='BLOCKED')
    return {'threshold_ok':ok_threshold,'empty_ok':ok_empty,'verdict_ok':ok_verdict}


def main():
    r={'p1146':audit_p1146(),'p1147':audit_p1147(),'p1148':audit_p1148()}
    passed=all(all(v.values()) for v in r.values())
    out={"packet":"P1149","as_of":"2026-05-10","checks":r,"overall_pass":passed}
    out_path=GEN/'p1149_strict_selector_probe_reproducibility_audit_summary.json'
    out_path.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding='utf-8')
    print(f"[P1149] overall_pass={passed} wrote {out_path}")
    if not passed:
        raise SystemExit(1)

if __name__=='__main__':
    main()
