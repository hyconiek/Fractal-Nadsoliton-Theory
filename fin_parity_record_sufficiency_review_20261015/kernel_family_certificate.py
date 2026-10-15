#!/usr/bin/env python3
"""Exact diagnostic-family algebra under inherited spectral boxes; not KO execution."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPECTRA = ROOT/"fin_rank7_followup/inputs/fin_handoff_audit/results.json"


def product(a,b):
    return [[sum((a[i][k]*b[k][j] for k in range(12)),F(0)) for j in range(12)] for i in range(12)]


def build():
    # P_k(i,j)=cos(2*pi*k*(i-j)/12)/6 for these two rational cosine tables.
    c3 = [F(1),F(0),F(-1),F(0)]
    c4 = [F(1),F(-1,2),F(-1,2)]
    P3 = [[c3[(i-j)%4]/6 for j in range(12)] for i in range(12)]
    P4 = [[c4[(i-j)%3]/6 for j in range(12)] for i in range(12)]
    B = [[P3[i][j]-P4[i][j] for j in range(12)] for i in range(12)]
    zero = [[F(0)]*12 for _ in range(12)]
    assert product(P3,P3)==P3 and product(P4,P4)==P4 and product(P3,P4)==zero
    assert sum(P3[i][i] for i in range(12))==sum(P4[i][i] for i in range(12))==2
    assert all(sum(row)==0 for row in B) and all(B[i][i]==0 for i in range(12))
    assert all(sum(B[i][j]*(-1)**j for j in range(12))==0 for i in range(12))
    assert B!=zero and B==[list(row) for row in zip(*B)]
    intervals = json.loads(SPECTRA.read_text())["exact"]["laplacian_intervals"]
    box = {str(k):[F(x) for x in intervals[k]] for k in (3,4,5,6)}
    eps = F(1,100)
    floors = {}
    for k in (3,4,5):
        lo,hi = box[str(k)]
        shift = eps if k in (3,4) else F(0)
        assert lo-shift>0 and hi+shift<box["6"][0]
        floors[str(k)] = 1-(hi+shift)/box["6"][0]
    return {
        "date":"2026-10-15","status":"EXACT_RATIONAL_FAMILY_SEED_UNDER_INHERITED_SPECTRA",
        "spectra_path":str(SPECTRA.relative_to(ROOT)),"spectra_sha256":hashlib.sha256(SPECTRA.read_bytes()).hexdigest(),
        "family":"A_epsilon=A7+epsilon*(P3-P4)","epsilon_interval":["-1/100","1/100"],
        "P3":[[str(x) for x in row] for row in P3],"P4":[[str(x) for x in row] for row in P4],
        "paid":["rational orthogonal rank-two Fourier projectors","symmetric centered nonzero perturbation","zero diagonal and trace change","parity direction annihilated","rank-seven positivity and unique top parity gap for the declared interval"],
        "gamma_fast_lower_bounds":{k:str(v) for k,v in floors.items()},
        "not_paid":["full LI/MT theorem across the family","common leading/fixed-resolution parity law across the family","new fast block-observation path theorem","physical source of A7"],
        "critical_scope_warning":"LI single-target forecast-risk equality does not prove inability of full parity histories to discriminate different kernels"}


if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser(); parser.add_argument("--output",type=Path); args=parser.parse_args()
    result=build()
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":result["status"],"family":result["family"],"gamma_floors":result["gamma_fast_lower_bounds"]},indent=2))
