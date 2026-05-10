#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P837 = GENERATED / "p837_current_strict_alpha_s_exact_required_form_statement_same_lane_exhaustion_boundary_audit_probe.json"
IN_F835 = GENERATED / "f835_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_statement_required_form_target_packet.json"
IN_F836 = GENERATED / "f836_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "f837_current_strict_alpha_s_exact_required_form_statement_continuation_boundary_packet.json"
OUT_SUMMARY = GENERATED / "f837_current_strict_alpha_s_exact_required_form_statement_continuation_boundary_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P837, IN_F835, IN_F836]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F837",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p837 = load_json(IN_P837)
    f835 = load_json(IN_F835)
    f836 = load_json(IN_F836)

    p837_theorem = p837.get("theorem_result") or {}
    f835_target = f835.get("target_object") or {}
    f836_target = f836.get("target_object") or {}

    if (
        p837.get("status")
        == "P837_CURRENT_STRICT_ALPHA_S_EXACT_REQUIRED_FORM_STATEMENT_SAME_LANE_EXHAUSTION_BOUNDARY_AUDIT_PROBE"
        and p837_theorem.get("exact_required_form_statement_exported_on_current_repo_state") is False
        and p837_theorem.get("same_lane_passive_groundwork_exhausted") is True
        and p837_theorem.get("next_honest_move_requires_continuation_boundary_export") is True
        and f835.get("status")
        == "F835_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET_NO_FALSE_PASS"
        and f836.get("status")
        == "F836_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F837_EXECUTED_CURRENT_STRICT_ALPHA_S_EXACT_REQUIRED_FORM_STATEMENT_CONTINUATION_BOUNDARY_PACKET_NO_FALSE_PASS"
    else:
        status = "F837_REQUIRES_REVIEW"

    artifact = {
        "stage": "F837",
        "packet_name": "CurrentStrictAlphaSLawfulExactRequiredFormStatementDomainAdmissionExactRequiredFormStatementContinuationBoundary_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p837_same_lane_exhaustion_boundary_probe": rel(IN_P837),
            "f835_exact_statement_required_form_target_packet": rel(IN_F835),
            "f836_exact_required_form_statement_target_packet": rel(IN_F836),
        },
        "exported_object": {
            "object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_continuation_boundary_v1",
            "goal": "Export the continuation boundary after the current same-lane passive groundwork under the missing exact_required_form_statement_ref is exhausted on the lawful exact-required-form-statement domain-admission lane.",
            "same_lane_exhaustion_boundary_ref": p837_theorem.get("boundary_name"),
            "current_exact_required_form_statement_target_ref": f836_target.get("object_id"),
            "upstream_exact_statement_required_form_target_ref": f835_target.get("object_id"),
            "admitted_next_move_classes": [
                "export_one_genuinely_new_same_lane_exact_required_form_statement_source",
                "shift_to_a_different_required_form_statement_provider_class_lane",
            ],
            "forbidden_repetition_clause": "no_further_same_level_passive_split_under_unchanged_exact_required_form_statement_ref_on_current_lawful_domain_admission_lane_as_primary_move",
            "scope": "strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_continuation_boundary_only",
        },
        "current_honest_reading": [
            "The repo now exports an explicit continuation boundary for the current lawful exact-required-form-statement domain-admission lane.",
            "The current passive chain has been exported as far as the present honest local split goes on that lane.",
            "The next honest move is now constrained to a genuinely new exact required-form statement source or a different provider class, not another fabricated passive split.",
        ],
        "recommended_next_move": "Attack one genuinely new same-lane exact required-form statement source candidate for the current lawful exact-required-form-statement domain-admission lane before attempting any further promotion on the current lane.",
        "hard_limits": [
            "Does not claim that the exact required-form statement already exists.",
            "Does not claim that the T213/T216 lane already enters the alpha_s exact-required-form-statement domain.",
            "Does not claim that any neighboring statement or form slot silently discharges the new lane.",
            "Does not claim that provider-class shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F837",
        "status": status,
        "as_of": AS_OF,
        "exported_object_id": artifact["exported_object"]["object_id"],
        "current_exact_required_form_statement_target_ref": artifact["exported_object"]["current_exact_required_form_statement_target_ref"],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
