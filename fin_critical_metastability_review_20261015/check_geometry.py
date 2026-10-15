#!/usr/bin/env python3
"""Exact MT geometric payments, conditional on the admitted spectral boxes."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parents[1]
SPECTRUM=ROOT/"fin_rank7_followup/inputs/fin_handoff_audit/results.json"

def build(epsilon=F(1,100)):
    expected=json.loads((ROOT/"fin_intrinsic_theory_review_20261014/verification.json").read_text())["input_hashes"][str(SPECTRUM.relative_to(ROOT))]
    assert hashlib.sha256(SPECTRUM.read_bytes()).hexdigest()==expected
    raw=json.loads(SPECTRUM.read_text())["exact"]["laplacian_intervals"]
    boxes={k:tuple(map(F,raw[k])) for k in (3,4,5,6)}
    assert F(2)<boxes[6][0]<=boxes[6][1]<F(12,5)
    gamma=F(93,5000)
    assert 1-max(boxes[k][1] for k in (3,4,5))/boxes[6][0]>gamma
    delta=epsilon**2;h=F(1,10**10);b=h/2
    pmin=(1-epsilon)/12-delta;pmax=(1+epsilon)/12+delta
    hessian=1/pmax-12*(1-gamma)-F(1,1000)*F(12,5)
    assert hessian>F(2,25)
    cm=1/(24*pmin*pmin);gradient=F(2,25)-cm*delta
    assert cm<F(307,50) and gradient>F(793,10000)
    soft=epsilon**4/24;fast=F(1,25)*delta**2
    assert soft>h and fast>h and b>0
    return {"date":"2026-10-15","status":"EXACT_RATIONAL_PAYMENTS_UNDER_INHERITED_SPECTRA",
            "spectral_path":str(SPECTRUM.relative_to(ROOT)),"spectral_sha256":expected,
            "epsilon":str(epsilon),"delta":str(delta),"h":str(h),"b":str(b),
            "gain_deviation_cap":"1/1000","gamma_floor":str(gamma),"p_min":str(pmin),"p_max":str(pmax),
            "fast_hessian_floor":str(hessian),"entropy_m_derivative_constant":str(cm),
            "gradient_fast_margin":str(gradient),"limiting_soft_boundary_floor":str(soft),
            "limiting_fast_boundary_floor":str(fast),
            "scope":"geometric constants and inherited premise payments only; N0 and root-path theorem are analytic, not certified by this script"}

def main():
    p=argparse.ArgumentParser();p.add_argument("--output",type=Path);args=p.parse_args()
    result=build()
    if args.output:args.output.write_text(json.dumps(result,indent=2)+"\n")
    else:assert json.loads(Path(__file__).with_name("GEOMETRY_CERTIFICATE.json").read_text())==result
    try:build(F(1,10))
    except AssertionError:rejected=True
    else:rejected=False
    assert rejected
    print(json.dumps({"status":result["status"],"H_lower":result["fast_hessian_floor"],"reject_oversized_cylinder":rejected}))
if __name__=="__main__":main()
