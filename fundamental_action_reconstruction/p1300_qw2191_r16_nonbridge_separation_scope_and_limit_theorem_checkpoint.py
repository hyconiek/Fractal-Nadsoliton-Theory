#!/usr/bin/env python3
"""P1300: R16 non-bridge separation scope and limit theorem checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1299", type=Path, default=GEN / "p1299_qw2191_r15_axiom_tagged_closure_policy_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1300_qw2191_r16_nonbridge_separation_scope_and_limit_theorem_summary.json")
    args = parser.parse_args()

    p1299 = _read(args.p1299)
    if p1299.get("next_priority") != "R16_NONBRIDGE_SEPARATION_SCOPE_AND_LIMIT_THEOREM":
        raise SystemExit("P1300 requires next_priority=R16_NONBRIDGE_SEPARATION_SCOPE_AND_LIMIT_THEOREM from P1299.")

    policy = p1299.get("policy", {})
    if policy.get("preferred_resolution_path") != "NB1_NONBRIDGE":
        raise SystemExit("P1300 requires preferred_resolution_path=NB1_NONBRIDGE.")
    if policy.get("bridge_to_legacy_allowed") is not False:
        raise SystemExit("P1300 requires bridge_to_legacy_allowed=false.")

    theorem = {
        "id": "NB_SCOPE_LIMIT_R16",
        "statement": "If legacy kernel is non-strict historical, strict closure claims cannot be transferred via bridge-to-legacy.",
        "admissible_claims": [
            "strict_operational_predictions_under_declared_assumptions",
            "nonstrict_axiom_tagged_global_narratives",
            "nb1_nonbridge_consistency_claims",
        ],
        "prohibited_claims": [
            "strict_core_closure_via_legacy_bridge",
            "implicit_role_transfer_from_legacy_parameters",
            "global_strict_closure_without_selector_source",
        ],
        "status": "SCOPE_LIMIT_DRAFTED",
    }

    out = {
        "packet": "P1300",
        "as_of": "2026-05-11",
        "lane": "NB1_NONBRIDGE_TRACK",
        "input": {"p1299": str(args.p1299)},
        "r16_nonbridge_scope_limit_theorem": theorem,
        "next_priority": "R17_NB1_FORMAL_NONTRANSFER_THEOREM_DRAFT",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1300] wrote {args.out}; status={theorem['status']}")


if __name__ == "__main__":
    main()
