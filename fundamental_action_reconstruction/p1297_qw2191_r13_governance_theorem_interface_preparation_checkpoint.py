#!/usr/bin/env python3
"""P1297: R13 governance-theorem interface preparation checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1296", type=Path, default=GEN / "p1296_qw2191_r12_strict_selector_closure_motion_review_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1297_qw2191_r13_governance_theorem_interface_preparation_summary.json")
    args = parser.parse_args()

    p1296 = _read(args.p1296)
    if p1296.get("next_priority") != "R13_GOVERNANCE_THEOREM_INTERFACE_PREPARATION":
        raise SystemExit("P1297 requires next_priority=R13_GOVERNANCE_THEOREM_INTERFACE_PREPARATION from P1296.")

    if p1296.get("r12_motion_review", {}).get("review_result") != "CONDITIONAL_HOLD":
        raise SystemExit("P1297 requires CONDITIONAL_HOLD review result from P1296.")

    interface = {
        "governance_targets": ["B1_BRIDGE", "NB1_NONBRIDGE"],
        "strict_lane_exports": [
            "THM_R9_STRICT_SELECTOR_SOURCE_A",
            "R10_COUNTERMODEL_SWEEP_PASS",
            "R11_PEER_REPLAY_PASS",
        ],
        "required_mapping": [
            "strict_result_to_governance_domain_map",
            "role_transfer_permission_boundary",
            "closure_claim_eligibility_predicate",
        ],
        "status": "DRAFT_READY",
    }

    out = {
        "packet": "P1297",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1296": str(args.p1296)},
        "r13_interface": interface,
        "closure_policy": {
            "strict_core_selector_closure_allowed": False,
            "global_qw2191_closure_allowed": False,
        },
        "next_priority": "R14_B1_NB1_INTERFACE_THEOREM_DRAFT",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1297] wrote {args.out}; interface_status={interface['status']}")


if __name__ == "__main__":
    main()
