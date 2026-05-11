#!/usr/bin/env python3
"""P1202: first CAS bridge probe for W4 'expand' operation."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def main() -> None:
    p1201 = json.loads((GEN / "p1201_w4_trace_hash_probe_summary.json").read_text(encoding="utf-8"))

    sympy_available = False
    expand_ok = False
    expanded_expr = None
    error = None

    try:
        import sympy as sp  # type: ignore

        psi4, a = sp.symbols("psi4 a")
        expr = a * psi4**2 - a * psi4**2
        expanded = sp.expand(expr)
        sympy_available = True
        expand_ok = True
        expanded_expr = str(expanded)
    except Exception as exc:  # runtime probe only
        error = str(exc)

    out = {
        "packet": "P1202",
        "as_of": "2026-05-10",
        "input_auditability_ready": bool(p1201.get("auditability_ready", False)),
        "sympy_available": sympy_available,
        "expand_ok": expand_ok,
        "expanded_expr": expanded_expr,
        "cas_error": error,
        "partial_symbolic_execution_ready": bool(sympy_available and expand_ok),
        "strict_closure_claim_allowed": False,
        "note": "Only first expand bridge is probed; no W4 discharge claim.",
    }

    out_path = GEN / "p1202_w4_cas_expand_bridge_probe_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1202] sympy_available={sympy_available} expand_ok={expand_ok} wrote {out_path}")


if __name__ == "__main__":
    main()
