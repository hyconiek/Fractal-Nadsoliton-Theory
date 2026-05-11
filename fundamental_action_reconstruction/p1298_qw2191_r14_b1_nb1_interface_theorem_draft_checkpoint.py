#!/usr/bin/env python3
"""P1298: R14 B1/NB1 interface theorem draft checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1297", type=Path, default=GEN / "p1297_qw2191_r13_governance_theorem_interface_preparation_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1298_qw2191_r14_b1_nb1_interface_theorem_draft_summary.json")
    args = parser.parse_args()

    p1297 = _read(args.p1297)
    if p1297.get("next_priority") != "R14_B1_NB1_INTERFACE_THEOREM_DRAFT":
        raise SystemExit("P1298 requires next_priority=R14_B1_NB1_INTERFACE_THEOREM_DRAFT from P1297.")

    if p1297.get("r13_interface", {}).get("status") != "DRAFT_READY":
        raise SystemExit("P1298 requires DRAFT_READY interface status from P1297.")

    theorem_interface_draft = {
        "bridge_branch": {
            "id": "B1",
            "claim": "Export explicit legacy->strict bridge map with role-transfer admissibility constraints.",
            "required_obligations": [
                "kernel_identification_or_bounded_deviation_theorem",
                "alpha_geo_beta_tors_role_transfer_justification",
                "selector_closure_eligibility_after_bridge",
            ],
        },
        "nonbridge_branch": {
            "id": "NB1",
            "claim": "Prove persistent non-equivalence with explicit prohibition boundary for role transfer.",
            "required_obligations": [
                "non_identifiability_or_obstruction_certificate",
                "legacy_role_nontransfer_theorem",
                "strict_lane_operational_scope_preservation",
            ],
        },
        "global_policy": {
            "strict_core_closure_claim": "NOT_ALLOWED",
            "global_closure_claim": "NOT_ALLOWED",
            "status_until_b1_or_nb1": "OPEN",
        },
        "status": "THEOREM_INTERFACE_DRAFTED",
    }

    out = {
        "packet": "P1298",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1297": str(args.p1297)},
        "r14_interface_theorem": theorem_interface_draft,
        "next_priority": "R15_B1_NB1_OBLIGATION_MATRIX_AND_PROOF_PLAN",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1298] wrote {args.out}; status={theorem_interface_draft['status']}")


if __name__ == "__main__":
    main()
