#!/usr/bin/env python3
"""P1174 safe-region gate for candidate pre-filtering.

Checks candidate hint parameters against P1172 safe bounds.
"""
from __future__ import annotations
import json, sys
from pathlib import Path

ROOT=Path(__file__).resolve().parent
GEN=ROOT/'generated'


def main():
    if len(sys.argv)!=2:
        print('usage: p1174_safe_region_gate.py <candidate_json>')
        raise SystemExit(2)
    cand=Path(sys.argv[1]).resolve()
    c=json.loads(cand.read_text())
    b=json.loads((GEN/'p1172_safe_operating_region_summary.json').read_text())['safe_bounds']

    checks={
      'omega_in_bounds': b['omega']['min'] <= c.get('omega_hint',0) <= b['omega']['max'],
      'phi_in_bounds': b['phi']['min'] <= c.get('phi_hint',0) <= b['phi']['max'],
      'sigma_in_bounds': b['sigma']['min'] <= c.get('sigma_hint',0) <= b['sigma']['max'],
      'kappa_in_bounds': b['kappa']['min'] <= c.get('kappa_hint',0) <= b['kappa']['max'],
    }
    passed=all(checks.values())
    out={"packet":"P1174","as_of":"2026-05-10","candidate":str(cand),"checks":checks,"safe_region_pass":passed}
    outp=GEN/'p1174_safe_region_gate_summary.json'
    outp.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(f"[P1174] safe_region_pass={passed} wrote {outp}")
    if not passed:
        raise SystemExit(1)

if __name__=='__main__':
    main()
