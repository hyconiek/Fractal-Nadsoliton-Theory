#!/usr/bin/env python3
"""P1209: open symbolic-attempt gate when validated CAS artifact is ready."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1205", type=Path, default=GEN / "p1205_w4_sympy_cas_runner_summary.json")
    parser.add_argument("--p1206", type=Path, default=GEN / "p1206_p1205_artifact_validator_summary.json")
    parser.add_argument("--p1213", type=Path, default=GEN / "p1213_external_cas_certificate_contract_summary.json")
    parser.add_argument("--allow-uncertified-artifact", action="store_true")
    parser.add_argument("--assume-strict-artifact", action="store_true")
    parser.add_argument("--p1205-only", action="store_true")
    parser.add_argument("--out", type=Path, default=GEN / "p1209_symbolic_attempt_gate_open_summary.json")
    args = parser.parse_args()

    p1205 = json.loads(args.p1205.read_text(encoding="utf-8"))
    p1206 = json.loads(args.p1206.read_text(encoding="utf-8"))
    p1213 = json.loads(args.p1213.read_text(encoding="utf-8"))

    p1205_base_ready = bool(
        p1205.get("sympy_available")
        and isinstance(p1205.get("trace_payload"), dict)
        and isinstance(p1205.get("trace_hash_sha256"), str)
    )

    if args.p1205_only:
        gate_open = p1205_base_ready
        gate_mode = "P1205_ONLY"
        base_ready = p1205_base_ready
        certificate_ready = False
    else:
        base_ready = bool(
            p1206.get("schema_ok")
            and p1206.get("sympy_available")
            and p1206.get("has_trace_payload")
            and p1206.get("has_trace_hash")
            and p1206.get("cas_ready_artifact")
        )
        certificate_ready = bool(p1213.get("external_certificate_contract_pass", False))

        strict_gate_open = bool(base_ready and certificate_ready)
        assumed_strict_gate_open = bool(args.assume_strict_artifact and base_ready)
        provisional_gate_open = bool(args.allow_uncertified_artifact and base_ready)
        gate_open = bool(strict_gate_open or assumed_strict_gate_open or provisional_gate_open)

        gate_mode = (
            "STRICT_CERTIFIED" if strict_gate_open else (
                "STRICT_ASSUMED" if assumed_strict_gate_open else (
                    "PROVISIONAL_UNCERTIFIED" if provisional_gate_open else "BLOCKED"
                )
            )
        )

    out = {
        "packet": "P1209",
        "as_of": "2026-05-11",
        "input_packets": ["P1205", "P1206", "P1213"],
        "p1205_only": bool(args.p1205_only),
        "p1205_base_ready": p1205_base_ready,
        "base_ready": base_ready,
        "certificate_ready": certificate_ready,
        "allow_uncertified_artifact": bool(args.allow_uncertified_artifact),
        "assume_strict_artifact": bool(args.assume_strict_artifact),
        "strict_assumption_applied": bool(args.assume_strict_artifact and not certificate_ready and base_ready),
        "gate_mode": gate_mode,
        "symbolic_attempt_gate_open": gate_open,
        "next_action": "RUN_CONTROLLED_W4_SYMBOLIC_ATTEMPT" if gate_open else "BLOCKED_WAIT_FOR_CAS_READY_ARTIFACT",
        "theory_closure_status": "OPEN",
        "strict_closure_claim_allowed": False,
        "note": "Gate opening is procedural and does not imply W4 discharge or ToE closure.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1209] symbolic_attempt_gate_open={gate_open} gate_mode={gate_mode} wrote {args.out}")


if __name__ == "__main__":
    main()
