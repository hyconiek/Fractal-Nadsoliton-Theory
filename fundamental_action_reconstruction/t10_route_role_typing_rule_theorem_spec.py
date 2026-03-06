from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    return json.loads((GEN / name).read_text())


def main() -> None:
    t9 = load("t9_route_admissibility_grammar_discharge_attempt_summary.json")

    out = {
        "status": "T10_PACKET_READY_ROUTE_ROLE_TYPING_RULE_THEOREM_SPEC_NO_FALSE_PASS",
        "goal": "Write a packet-ready theorem spec for the missing route-role typing rule isolated by T9.",
        "supports": {
            "T9": t9["discharge_target"],
            "T9_B1": t9["residual_blockers"]["T9_B1"],
        },
        "target_theorem": "route_role_typing_rule_for_current_selector_track",
        "route_roles_under_typing": [
            "raw_overlap",
            "formula",
            "representative",
            "downstream_schema",
            "source_skeleton",
            "strict_to_axiom_bridge"
        ],
        "residual_blockers": {
            "T10_B1": "the_route_role_typing_rule_is_specified_but_not_discharged_for_the_current_selector_track"
        },
        "forbidden_claims": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_the_route_role_typing_rule_is_already_proved",
            "no_claim_that_t8_t6_t4_or_t1_are_already_discharged",
            "no_claim_that_qw_2191_is_discharged"
        ]
    }

    out_path = GEN / "t10_route_role_typing_rule_theorem_spec_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True))
    print(out_path)


if __name__ == "__main__":
    main()
