#!/usr/bin/env python3
"""Exact rational payments for an equilibrium exclusion bound, not mixing."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json,math

ROOT=Path(__file__).resolve().parents[1]
SPECTRUM=ROOT/"fin_rank7_followup/inputs/fin_handoff_audit/results.json"

def check_public(certificate, expected):
    assert certificate["spectral_input_sha256"]==expected["spectral_input_sha256"]
    for field in ("active_intervals","A00_lower","A00_upper","exp_5_over_2_partial_sum",
                  "log12_upper","critical_delta0_upper","gain_interval","occupation_radius",
                  "uniform_ball_F_lower","one_colour_F_upper","exclusion_exponent","bound"):
        assert certificate[field]==expected[field],field

def build():
    inherited=json.loads((ROOT/"fin_intrinsic_theory_review_20261014/verification.json").read_text())
    assert hashlib.sha256(SPECTRUM.read_bytes()).hexdigest()==inherited["input_hashes"][str(SPECTRUM.relative_to(ROOT))]
    raw=json.loads(SPECTRUM.read_text())["exact"]["laplacian_intervals"]
    intervals={k:tuple(map(F,raw[k])) for k in (3,4,5,6)}
    assert all(0<lo<=hi for lo,hi in intervals.values())
    assert all(intervals[k][1]<intervals[6][0] for k in (3,4,5))
    diaglo=sum(intervals[k][0] for k in (3,4,5))/6+intervals[6][0]/12
    diaghi=sum(intervals[k][1] for k in (3,4,5))/6+intervals[6][1]/12
    assert diaglo>F(5,4)
    assert intervals[6][1]<F(5,2)
    # Positive Taylor terms prove exp(5/2)>12, hence log(12)<5/2.
    exp_lower=sum(F(5,2)**k/math.factorial(k) for k in range(13))
    assert exp_lower>12
    assert 4*sum(intervals[k][0] for k in (3,4,5))>11*intervals[6][1]
    # At gc: F(e0)=log12-1/2-(lambda3+lambda4+lambda5)/lambda6 < -3/4.
    # For every 5<=g<=6: F(e0)<5/2-5*(5/4)/2=-5/8.
    # On ||p-u||<=1/10: D>=0 and ||A||<5/2 give F>=-3/40.
    one_colour=F(5,2)-F(5,2)*F(5,4)
    ball=-F(6,2)*F(5,2)*F(1,10)**2
    exponent=ball-one_colour
    assert one_colour==F(-5,8) and ball==F(-3,40) and exponent==F(11,20)
    return {"date":"2026-10-14","status":"EXACT_FRACTION_PAYMENTS_CHECKED",
            "spectral_input":str(SPECTRUM.relative_to(ROOT)),
            "spectral_input_sha256":hashlib.sha256(SPECTRUM.read_bytes()).hexdigest(),
            "active_intervals":{str(k):[str(x) for x in v] for k,v in intervals.items()},
            "A00_lower":str(diaglo),"A00_upper":str(diaghi),
            "exp_5_over_2_partial_sum":str(exp_lower),"log12_upper":"5/2",
            "critical_delta0_upper":"-3/4","gain_interval":["5","6"],
            "occupation_radius":"1/10","uniform_ball_F_lower":str(ball),
            "one_colour_F_upper":str(one_colour),"exclusion_exponent":str(exponent),
            "bound":"pi_N,g(||p-u||_2<=1/10) <= min(1,(N+1)^12 exp(-11 N/20))",
            "scope":"all positive integer N and g in [5,6], supplied Gibbs model with accepted spectral boxes; not dynamical error or mixing"}

def main():
    parser=argparse.ArgumentParser();parser.add_argument("--output",type=Path);args=parser.parse_args()
    expected=build()
    path=Path(__file__).with_name("GIBBS_EXCLUSION_CERTIFICATE.json")
    if args.output:
        args.output.write_text(json.dumps(expected,indent=2)+"\n")
    else:
        check_public(json.loads(path.read_text()),expected)
    altered=dict(expected);altered["exclusion_exponent"]="3/5"
    try:check_public(altered,expected)
    except AssertionError:rejected=True
    else:rejected=False
    assert rejected
    print(json.dumps({"status":expected["status"],"exponent":expected["exclusion_exponent"],
                      "reject_altered_exponent":rejected,"bound":expected["bound"]},indent=2))

if __name__=="__main__":main()
