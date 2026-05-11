#!/usr/bin/env python3
"""P1205: external SymPy CAS runner for W4 pre-discharge artifacts.

Run this in an environment where sympy is installed.
It writes: fundamental_action_reconstruction/generated/p1205_w4_sympy_cas_runner_summary.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Run W4 SymPy CAS pre-discharge checks and export JSON artifact.")
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parent / "generated" / "p1205_w4_sympy_cas_runner_summary.json",
        help="Output JSON path.",
    )
    args = parser.parse_args()

    try:
        import sympy as sp  # type: ignore
    except Exception as exc:
        payload = {
            "packet": "P1205",
            "as_of": "2026-05-10",
            "sympy_available": False,
            "error": str(exc),
            "strict_closure_claim_allowed": False,
            "note": "SymPy missing in runtime environment.",
        }
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"[P1205] sympy unavailable, wrote {args.out}")
        return

    psi4, a, b = sp.symbols("psi4 a b")

    # Control expressions for auditable CAS readiness.
    expr_cancel = a * psi4**2 - a * psi4**2
    expr_expand = (a + b) ** 2

    expanded_cancel = sp.expand(expr_cancel)
    expanded_binom = sp.expand(expr_expand)
    grouped_binom = sp.collect(expanded_binom, [a, b])

    reduced_zero_ok = bool(sp.simplify(expanded_cancel) == 0)

    trace_payload = {
        "cancel_expression": str(expr_cancel),
        "expanded_cancel": str(expanded_cancel),
        "expand_expression": str(expr_expand),
        "expanded_binom": str(expanded_binom),
        "grouped_binom": str(grouped_binom),
        "reduced_zero_ok": reduced_zero_ok,
    }
    trace_hash = hashlib.sha256(json.dumps(trace_payload, sort_keys=True).encode("utf-8")).hexdigest()

    payload = {
        "packet": "P1205",
        "as_of": "2026-05-10",
        "sympy_available": True,
        "sympy_version": sp.__version__,
        "trace_payload": trace_payload,
        "trace_hash_sha256": trace_hash,
        "partial_symbolic_execution_ready": True,
        "w4_discharge_performed": False,
        "strict_closure_claim_allowed": False,
        "note": "CAS pre-discharge artifact only; no W4 discharge claim.",
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1205] sympy OK, wrote {args.out}")


if __name__ == "__main__":
    main()
