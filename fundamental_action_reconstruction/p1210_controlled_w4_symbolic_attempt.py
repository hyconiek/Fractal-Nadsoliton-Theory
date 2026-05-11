#!/usr/bin/env python3
"""P1210: controlled symbolic attempt gate (no false pass)."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1205", type=Path, default=GEN / "p1205_w4_sympy_cas_runner_summary.json")
    parser.add_argument("--p1209", type=Path, default=GEN / "p1209_symbolic_attempt_gate_open_summary.json")
    parser.add_argument("--p1212", type=Path, default=GEN / "p1212_p1205_trace_provenance_seal_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1210_controlled_w4_symbolic_attempt_summary.json")
    args = parser.parse_args()

    p1205 = json.loads(args.p1205.read_text(encoding="utf-8"))
    p1209 = json.loads(args.p1209.read_text(encoding="utf-8"))
    p1212 = json.loads(args.p1212.read_text(encoding="utf-8"))

    gate_open = bool(p1209.get("symbolic_attempt_gate_open", False))
    provenance_seal_pass = bool(p1212.get("provenance_seal_pass", False))
    trace = p1205.get("trace_payload") if isinstance(p1205.get("trace_payload"), dict) else {}

    reduced_zero_ok = bool(trace.get("reduced_zero_ok", False))
    has_hash = isinstance(p1205.get("trace_hash_sha256"), str)

    attempt_executed = gate_open and provenance_seal_pass and bool(trace)
    w4_discharged = bool(attempt_executed and reduced_zero_ok and has_hash)

    out = {
        "packet": "P1210",
        "as_of": "2026-05-11",
        "gate_open": gate_open,
        "provenance_seal_pass": provenance_seal_pass,
        "attempt_executed": attempt_executed,
        "reduced_zero_ok": reduced_zero_ok,
        "has_trace_hash": has_hash,
        "w4_status": "DISCHARGED" if w4_discharged else "OPEN",
        "w4_discharge_pass": w4_discharged,
        "strict_closure_claim_allowed": False,
        "note": "Controlled symbolic attempt result; closure remains forbidden regardless of local witness state.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1210] attempt_executed={attempt_executed} w4_discharge_pass={w4_discharged} wrote {args.out}")


if __name__ == "__main__":
    main()
