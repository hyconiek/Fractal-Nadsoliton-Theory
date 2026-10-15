#!/usr/bin/env python3
"""Ideal slope only: inherited rational spectra plus mpmath.iv, not native N0."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
SPECTRA = ROOT / "fin_rank7_followup/inputs/fin_handoff_audit/results.json"

def dyadic(t):
    sign, mantissa, exponent, _ = t
    value = F(int(mantissa)) * (F(2) ** int(exponent))
    return -value if sign else value

def bounds(x):
    return [dyadic(t) for t in x._mpi_]

def public(x, digits=18):
    lo, hi = bounds(x)
    scale = 10 ** digits
    left = lo.numerator * scale // lo.denominator
    right = -((-hi.numerator * scale) // hi.denominator)
    def decimal(n):
        sign = "-" if n < 0 else ""
        n = abs(n)
        return sign + str(n // scale) + "." + str(n % scale).zfill(digits)
    result = [decimal(left), decimal(right)]
    assert F(result[0]) <= lo <= hi <= F(result[1])
    return result

def build():
    mp.iv.dps = 80
    iv = mp.iv
    def rational(x):
        x = F(x)
        return iv.mpf(x.numerator) / x.denominator
    raw = json.loads(SPECTRA.read_text())["exact"]["laplacian_intervals"]
    lam = {k: iv.mpf([rational(raw[k][0]).a, rational(raw[k][1]).b]) for k in (3, 6)}
    gs = [1-lam[3]/lam[6], 1-(lam[3]+rational(".01"))/lam[6]]
    beta = [1/g-1 for g in gs]
    d = (beta[1]-beta[0])/(2*iv.pi)
    assert bounds(d)[0] > 0
    # MVT in v: dp/dv = 1/(pi sqrt(v)(1+v)).  Bound the complete
    # v interval, rather than subtract nearly equal atan evaluations.
    vmax = 1 + beta[1]/16
    relative_lower = 2/(iv.sqrt(vmax)*(1+vmax))
    assert bounds(relative_lower)[0] > F(1, 2)
    return {
        "date": "2026-10-15",
        "status": "VALIDATED_IDEAL_SLOPE_UNDER_INHERITED_SPECTRA_AND_IV_ASSUMPTIONS",
        "spectral_path": str(SPECTRA.relative_to(ROOT)),
        "spectral_sha256": hashlib.sha256(SPECTRA.read_bytes()).hexdigest(),
        "d_ang": public(d),
        "uniform_alpha_interval": ["0 (excluded)", "1/16"],
        "relative_gap_lower_factor": public(relative_lower),
        "uniform_ideal_gap": "p_control(alpha)-p_A7(alpha)>=(d_ang/2)*alpha",
        "rounding": "exact Fraction floor/ceiling of raw dyadic interval endpoints",
        "not_certified": ["native constants or usable N0", "physical gain or detector", "all source printed digits"]
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = build()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
