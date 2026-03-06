from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    return json.loads((GEN / name).read_text())


def main() -> None:
    t4 = load("t4_strict_core_export_completeness_principle_theorem_spec_summary.json")
    t5 = load("t5_export_completeness_principle_discharge_attempt_summary.json")

    out = {
        "status": "T6_PACKET_READY_ROUTE_FAMILY_CLOSURE_CERTIFICATE_THEOREM_SPEC_NO_FALSE_PASS",
        "goal": "Write a packet-ready theorem spec for the missing route-family closure certificate isolated by T5.",
        "supports": {
            "T4": t4["target_theorem"],
            "T5": t5["residual_blockers"]["T5_B1"]
        },
        "route_family_under_certificate": [
            "C32",
            "C33",
            "C34",
            "C49",
            "C50",
            "C51"
        ],
        "target_theorem": "route_family_closure_certificate_for_current_selector_track",
        "residual_blockers": {
            "T6_B1": "the_route_family_closure_certificate_is_specified_but_not_discharged_for_the_current_strict_core_selector_track"
        },
        "forbidden_claims": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_the_closure_certificate_is_already_proved",
            "no_claim_that_t4_is_already_discharged",
            "no_claim_that_t1_is_already_discharged",
            "no_claim_that_qw_2191_is_discharged"
        ]
    }

    out_path = GEN / "t6_route_family_closure_certificate_theorem_spec_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True))
    print(out_path)


if __name__ == "__main__":
    main()
