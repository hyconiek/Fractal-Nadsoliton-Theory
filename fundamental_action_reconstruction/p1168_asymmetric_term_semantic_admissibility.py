#!/usr/bin/env python3
"""P1168 semantic admissibility check for asymmetric selector term candidate.

Assesses whether top P1166 candidate satisfies explicit strict-side semantic
criteria before any stronger theoretical interpretation.
"""
from __future__ import annotations
import json, math
from pathlib import Path

ROOT=Path(__file__).resolve().parent
GEN=ROOT/'generated'


def main():
    top=json.loads((GEN/'p1166_sigma_kappa_sweep_summary.json').read_text())['best']
    fp=json.loads((GEN/'p1166_sigma_kappa_sweep_summary.json').read_text())['fixed_params']
    sigma,kappa=top['sigma'],top['kappa']

    # semantic checks for A(d)=sigma*(1-exp(-kappa d)) on d>=0
    d_test=[0.0,0.1,1.0,5.0,24.0,48.0]
    vals=[sigma*(1.0-math.exp(-kappa*d)) for d in d_test]

    checks={
      "sigma_nonnegative": sigma>=0,
      "kappa_positive": kappa>0,
      "A_zero_at_origin": abs(vals[0])<1e-15,
      "A_bounded_by_sigma": all((0.0 <= v <= sigma+1e-12) for v in vals),
      "A_monotone_in_d": all(vals[i]<=vals[i+1]+1e-12 for i in range(len(vals)-1)),
      "strict_side_nonclaim_context": True,
    }
    overall=all(checks.values())

    out={
      "packet":"P1168","as_of":"2026-05-10",
      "candidate":{"sigma":sigma,"kappa":kappa,**fp},
      "A_samples": [{"d":d,"A":v} for d,v in zip(d_test,vals)],
      "semantic_checks":checks,
      "semantic_admissible":overall,
      "note":"Semantic admissibility only; no closure/QW-2191 discharge claim."
    }
    outp=GEN/'p1168_asymmetric_term_semantic_admissibility_summary.json'
    outp.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(f"[P1168] semantic_admissible={overall} wrote {outp}")

if __name__=='__main__':
    main()
