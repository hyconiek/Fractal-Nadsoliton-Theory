#!/usr/bin/env python3
"""P1148: strict phase-shifted selector family probe.

Tests whether adding a phase offset delta (as explicit extra premise) can move
first sign flip beyond a finite operational domain, without claiming theorem
closure or QW-2191 discharge.
"""
from __future__ import annotations
import json, math
from pathlib import Path

OMEGA=0.18575
PHI=0.16250
BETA=1.0
ETA=1.8
ALPHA=4.0*math.log(2.0)
D_MAX=24.0
STEP=0.1


def s_delta(d: float, delta: float) -> float:
    return math.exp(-ALPHA*d)*math.cos(OMEGA*d+PHI+delta)/(1.0+BETA*d**ETA)


def first_flip_on_domain(delta: float):
    prev = s_delta(0.0, delta)
    d=STEP
    while d<=D_MAX:
        cur=s_delta(d,delta)
        if prev*cur<0:
            return d
        if cur!=0:
            prev=cur
        d+=STEP
    return None


def main():
    deltas=[-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0]
    results=[]
    for delta in deltas:
        flip=first_flip_on_domain(delta)
        results.append({"delta":delta,"first_flip_d":flip,"no_flip_on_domain":flip is None})

    # Threshold to push first cosine zero beyond D_MAX: omega*d+phi+delta < pi/2 for d<=D_MAX
    delta_threshold=(math.pi/2)-PHI-OMEGA*D_MAX

    admissible=[r for r in results if r["no_flip_on_domain"]]
    out={
        "packet":"P1148",
        "as_of":"2026-05-10",
        "domain":{"d_max":D_MAX,"step":STEP},
        "base_params":{"alpha":ALPHA,"omega":OMEGA,"phi":PHI,"beta":BETA,"eta":ETA},
        "delta_threshold_for_no_flip_on_domain":delta_threshold,
        "scan_results":results,
        "admissible_no_flip_deltas":admissible,
        "audit":{
            "finite_domain_no_flip_possible": len(admissible)>0,
            "strict_core_closure_claim":False,
            "qw_2191_discharge":False,
            "verdict":"CANDIDATE_PREMISE_ONLY" if len(admissible)>0 else "BLOCKED"
        }
    }
    p=Path(__file__).resolve().parent/"generated"/"p1148_strict_phase_shifted_selector_family_probe_summary.json"
    p.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(f"[P1148] wrote {p}")

if __name__=='__main__':
    main()
