#!/usr/bin/env python3
"""P1258: strict-only operational lane commitment with no implicit ontology transfer.

Implements the requested strategy decision: continue computations on K_strict_gate only,
while preserving guardrails that forbid claiming legacy-equivalence or global closure
without B1/NB1 discharge.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1256", type=Path, default=GEN / "p1256_bridge_nonbridge_theorem_spec_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1258_strict_only_operational_lane_commitment_summary.json")
    args = parser.parse_args()

    p1256 = _read(args.p1256)

    policy_ok = p1256.get("closure_policy") == "STRICT_CLOSURE_FORBIDDEN_UNTIL_B1_OR_NB1_DISCHARGED"
    if not policy_ok:
        raise SystemExit("P1258 requires P1256 closure policy to remain active.")

    out = {
        "packet": "P1258",
        "as_of": "2026-05-11",
        "input": {"p1256": str(args.p1256)},
        "strategy_decision": {
            "operational_lane": "STRICT_ONLY",
            "kernel_in_active_use": "K_strict_gate",
            "reason": "legacy lane already extensively analyzed; proceed with strict operational program",
        },
        "hard_constraints": [
            "NO_IMPLICIT_LEGACY_EQUIVALENCE_CLAIM",
            "NO_LEGACY_PHYSICAL_ROLE_TRANSFER",
            "NO_STRICT_CORE_CLOSURE_UNTIL_B1_OR_NB1",
        ],
        "scientific_status": "ADVANCE_STRICT_COMPUTATIONAL_PROGRAM_WITH_EXPLICIT_SCOPE_LIMIT",
        "honest_limit_statement": (
            "Strict-only progress is operational. It does not by itself settle the legacy->strict identity question "
            "nor authorize full theory closure."
        ),
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1258] wrote {args.out}")


if __name__ == "__main__":
    main()
