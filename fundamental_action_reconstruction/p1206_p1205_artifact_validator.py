#!/usr/bin/env python3
"""P1206: validate external P1205 SymPy artifact before downstream use."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate p1205_w4_sympy_cas_runner_summary.json schema/readiness.")
    parser.add_argument(
        "--in",
        dest="inp",
        type=Path,
        default=GEN / "p1205_w4_sympy_cas_runner_summary.json",
        help="Input P1205 JSON artifact path.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=GEN / "p1206_p1205_artifact_validator_summary.json",
        help="Output validation summary path.",
    )
    args = parser.parse_args()

    data = json.loads(args.inp.read_text(encoding="utf-8"))

    required_base = {"packet", "as_of", "sympy_available", "strict_closure_claim_allowed"}
    missing = sorted(x for x in required_base if x not in data)

    sympy_available = bool(data.get("sympy_available", False))
    has_trace = isinstance(data.get("trace_payload"), dict)
    has_hash = isinstance(data.get("trace_hash_sha256"), str)

    schema_ok = not missing
    cas_ready_artifact = bool(schema_ok and sympy_available and has_trace and has_hash)

    out = {
        "packet": "P1206",
        "as_of": "2026-05-11",
        "input_path": str(args.inp),
        "schema_ok": schema_ok,
        "missing_required_fields": missing,
        "sympy_available": sympy_available,
        "has_trace_payload": has_trace,
        "has_trace_hash": has_hash,
        "cas_ready_artifact": cas_ready_artifact,
        "strict_closure_claim_allowed": False,
        "note": "Validation gate for externally generated P1205 CAS artifact.",
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1206] schema_ok={schema_ok} cas_ready_artifact={cas_ready_artifact} wrote {args.out}")


if __name__ == "__main__":
    main()
