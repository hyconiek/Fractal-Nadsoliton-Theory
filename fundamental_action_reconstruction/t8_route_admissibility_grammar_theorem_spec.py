from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    return json.loads((GEN / name).read_text())


def main() -> None:
    t7 = load("t7_route_family_closure_certificate_discharge_attempt_summary.json")

    out = {
        "status": "T8_PACKET_READY_ROUTE_ADMISSIBILITY_GRAMMAR_THEOREM_SPEC_NO_FALSE_PASS",
        "goal": "Write a packet-ready theorem spec for the missing route admissibility grammar isolated by T7.",
        "supports": {
            "T7": t7["discharge_target"],
            "T7_B1": t7["residual_blockers"]["T7_B1"],
        },
        "target_theorem": "route_admissibility_grammar_for_current_selector_track",
        "constructor_families_under_grammar": [
            "raw_overlap",
            "formula",
            "representative",
            "downstream_schema",
            "source_skeleton",
            "strict_to_axiom_bridge"
        ],
        "residual_blockers": {
            "T8_B1": "the_route_admissibility_grammar_is_specified_but_not_discharged_for_the_current_selector_track"
        },
        "forbidden_claims": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_the_admissibility_grammar_is_already_proved",
            "no_claim_that_t6_is_already_discharged",
            "no_claim_that_t4_or_t1_are_already_discharged",
            "no_claim_that_qw_2191_is_discharged"
        ]
    }

    out_path = GEN / "t8_route_admissibility_grammar_theorem_spec_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True))
    print(out_path)


if __name__ == "__main__":
    main()
