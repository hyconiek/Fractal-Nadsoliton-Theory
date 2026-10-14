#!/usr/bin/env python3
"""Stdlib exact checker of serialized coverage and TV arithmetic.

Probability-box validity is supplied by the independently replayable
mpmath.iv input construction; this checker does not verify transcendentals.
"""
import argparse
import copy
import json
from fractions import Fraction as F
from pathlib import Path


def check(c):
    g=c["gain_grid"]
    assert F(g["origin"])==F(297,100) and F(g["step"])==F(1,1000)
    n=g["cells_per_axis"]
    assert n==60 and g["closed_boxes"]==n*n
    assert F(g["origin"])+n*F(g["step"])==F(g["last_endpoint"])==F(303,100)
    laws={p["model"]:p["boxes"] for p in c["probability_boxes"]}
    assert set(laws)=={"FIN_A7","FLAT_P7_TRACE"}
    for vv in laws.values():
        assert len(vv)==n
        for row in vv:
            assert len(row)==12
            assert all(0<=F(lo)<=F(hi)<=1 for lo,hi in row)
            assert sum(F(x[0]) for x in row)<=1<=sum(F(x[1]) for x in row)
    assert len(c["cell_lower_bounds"])==n
    lower=[]
    for i,p in enumerate(laws["FIN_A7"]):
        assert len(c["cell_lower_bounds"][i])==n
        for j,q in enumerate(laws["FLAT_P7_TRACE"]):
            exact=sum((max(F(0),F(a[0])-F(b[1]),F(b[0])-F(a[1])) for a,b in zip(p,q)),F())/2
            saved=F(c["cell_lower_bounds"][i][j])
            assert 0<=saved<=exact
            lower.append(saved)
    reported=F(c["global_ideal_TV_lower"])
    assert reported<=min(lower) and reported>F(1,100)
    assert F(c["robust_observed_law_lower_after_two_detector_allowances"])<=reported-2*F(c["per_hypothesis_detector_output_TV_allowance"])
    assert F(c["robust_observed_law_lower_after_two_detector_allowances"])>0
    wp,wq=c["point_witness_probability_boxes"]
    assert len(wp)==len(wq)==12
    upper=sum((max(abs(F(a[0])-F(b[1])),abs(F(a[1])-F(b[0]))) for a,b in zip(wp,wq)),F())/2
    assert upper<=F(c["global_ideal_TV_upper_from_witness"])
    return len(lower)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("certificate",type=Path)
    ap.add_argument("--output",type=Path)
    args=ap.parse_args()
    c=json.loads(args.certificate.read_text())
    count=check(c)
    bad=copy.deepcopy(c)
    bad["cell_lower_bounds"].pop()
    try:
        check(bad)
    except AssertionError:
        rejected=True
    else:
        raise AssertionError("missing coverage row was accepted")
    result={"exact_fraction_checks_pass":True,"covered_boxes":count,
            "missing_coverage_mutation_rejected":rejected,
            "limit":"input probability enclosures assume audited mpmath.iv construction; no physical calibration"}
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result))


if __name__=="__main__":
    main()
