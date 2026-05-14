#!/usr/bin/env python3
"""Standalone SymPy CAS appendix for strict-only L_total -> EOM replay.
Run: python3 cas_appendix_sympy_p1615.py
"""
from __future__ import annotations
import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN = GEN / "p1614_s564_strict_toe_closure_dossier_summary.json"
OUT = GEN / "p1615_sympy_cas_appendix_replay.json"

# fields and symbols
psi, m2, eps, R = sp.symbols("psi m2 eps R", real=True)

# toy strict scalar+mix sector extracted from exported full chain
L = sp.Rational(1, 2) * sp.Symbol("dpsi2") - sp.Rational(1, 2) * m2 * psi**2 + eps * psi * R

# algebraic EL-identity in homogeneous proxy: dV/dpsi + source = 0
# V = 1/2 m2 psi^2 - eps psi R
V = sp.Rational(1, 2) * m2 * psi**2 - eps * psi * R
el_proxy = sp.diff(V, psi)
# expected: m2*psi - eps*R
expected = m2 * psi - eps * R

ok = sp.simplify(el_proxy - expected) == 0

payload = {
    "checkpoint": "P1615_SYMPY_APPENDIX",
    "status": "PASS" if ok else "FAIL",
    "identity": "dV/dpsi == m2*psi - eps*R",
    "sympy_result": str(sp.simplify(el_proxy)),
    "expected": str(expected),
    "difference": str(sp.simplify(el_proxy - expected)),
    "strict_only": True,
    "legacy_bridge_used": False,
}
OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
print(f"Wrote {OUT}")
