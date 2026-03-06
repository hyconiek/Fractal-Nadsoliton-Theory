from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    return json.loads((GEN / name).read_text())


def main() -> None:
    t8 = load("t8_route_admissibility_grammar_theorem_spec_summary.json")

    out = {
        "status": "T9_EXECUTED_T8_DISCHARGE_ATTEMPT_BLOCKED_BY_NO_FORMAL_ROUTE_ROLE_TYPING_RULE",
        "goal": "Attempt the first discharge of T8 by checking whether the current selector-track audit vocabulary already defines route admissibility by role.",
        "supports": {
            "T8": t8["target_theorem"],
            "constructor_families_under_grammar": t8["constructor_families_under_grammar"],
        },
        "discharge_target": "route_admissibility_grammar_for_current_selector_track",
        "established_partial_results": {
            "finite_explicit_route_role_vocabulary_present": True,
            "currently_audited_routes_have_stable_role_labels": True,
            "formal_route_role_typing_rule_present": False,
            "admissibility_exhausted_by_role_labels": False,
        },
        "residual_blockers": {
            "T9_B1": "no_formal_route_role_typing_rule_or_admissibility_by_role_declaration_showing_that_every_current_strict_core_theta_export_route_must_instantiate_exactly_one_of_the_six_named_route_roles"
        },
        "forbidden_claims": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_t8_is_discharged",
            "no_claim_that_route_role_labels_already_form_a_proved_admissibility_grammar",
            "no_claim_that_qw_2191_is_discharged"
        ]
    }

    out_path = GEN / "t9_route_admissibility_grammar_discharge_attempt_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True))
    print(out_path)


if __name__ == "__main__":
    main()
