#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P825 = GENERATED / "p825_current_strict_alpha_s_required_form_statement_class_candidate_supported_exact_required_form_statement_blocked.json"
IN_F824 = GENERATED / "f824_current_strict_alpha_s_exact_statement_required_form_target_packet.json"
IN_F823 = GENERATED / "f823_current_strict_alpha_s_exact_output_schema_statement_target_packet.json"
IN_F822 = GENERATED / "f822_current_strict_alpha_s_lawful_schema_domain_admission_output_schema_target_packet.json"
IN_F821 = GENERATED / "f821_current_strict_alpha_s_lawful_schema_domain_admission_target_packet.json"
IN_F820 = GENERATED / "f820_current_strict_alpha_s_schema_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F819 = GENERATED / "f819_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_schema_adapter_or_carrier_identification_rule_target_packet.json"
IN_F818 = GENERATED / "f818_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_schema_interface_target_packet.json"
IN_F814 = GENERATED / "f814_current_strict_alpha_s_strict_source_shannon_source_binding_selection_or_preference_rule_schema_target_packet.json"
IN_F812 = GENERATED / "f812_current_strict_alpha_s_strict_source_shannon_shift_to_domain_source_binding_selection_or_preference_rule_target_packet.json"
IN_F789 = GENERATED / "f789_current_strict_alpha_s_normalized_boundary_interface_target_packet.json"

OUT = GENERATED / "f825_current_strict_alpha_s_exact_required_form_statement_target_packet.json"
OUT_SUMMARY = GENERATED / "f825_current_strict_alpha_s_exact_required_form_statement_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [
        IN_P825,
        IN_F824,
        IN_F823,
        IN_F822,
        IN_F821,
        IN_F820,
        IN_F819,
        IN_F818,
        IN_F814,
        IN_F812,
        IN_F789,
    ]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F825",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p825 = load_json(IN_P825)
    f824 = load_json(IN_F824)
    f823 = load_json(IN_F823)
    f822 = load_json(IN_F822)
    f821 = load_json(IN_F821)
    f820 = load_json(IN_F820)
    f819 = load_json(IN_F819)
    f818 = load_json(IN_F818)
    f814 = load_json(IN_F814)
    f812 = load_json(IN_F812)
    f789 = load_json(IN_F789)

    p825_theorem = p825.get("theorem_result") or {}
    f824_target = f824.get("target_object") or {}
    f823_target = f823.get("target_object") or {}
    f822_target = f822.get("target_object") or {}
    f821_target = f821.get("target_object") or {}
    f820_target = f820.get("target_object") or {}
    f819_target = f819.get("target_object") or {}
    f818_target = f818.get("target_object") or {}
    f814_target = f814.get("target_object") or {}
    f812_target = f812.get("target_object") or {}
    f789_target = f789.get("target_interface") or {}

    if (
        p825.get("status")
        == "P825_CURRENT_STRICT_ALPHA_S_REQUIRED_FORM_STATEMENT_CLASS_CANDIDATE_SUPPORTED_EXACT_REQUIRED_FORM_STATEMENT_BLOCKED"
        and p825_theorem.get("required_form_statement_class_candidate_supported_now") is True
        and p825_theorem.get("exact_required_form_statement_exported_now") is False
        and p825_theorem.get("sharp_blocker_field") == "exact_required_form_statement_ref"
        and p825_theorem.get("next_honest_move_is_freeze_exact_required_form_statement_target") is True
        and f824.get("status")
        == "F824_EXECUTED_CURRENT_STRICT_ALPHA_S_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F825_EXECUTED_CURRENT_STRICT_ALPHA_S_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F825_REQUIRES_REVIEW"

    artifact = {
        "stage": "F825",
        "packet_name": "CurrentStrictAlphaSExactRequiredFormStatementTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p825_required_form_statement_probe": rel(IN_P825),
            "f824_required_form_target_packet": rel(IN_F824),
            "f823_exact_statement_target_packet": rel(IN_F823),
            "f822_output_schema_target_packet": rel(IN_F822),
            "f821_lawful_admission_target_packet": rel(IN_F821),
            "f820_boundary_target_packet": rel(IN_F820),
            "f819_rule_target_packet": rel(IN_F819),
            "f818_interface_target_packet": rel(IN_F818),
            "f814_downstream_schema_target_packet": rel(IN_F814),
            "f812_downstream_rule_target_packet": rel(IN_F812),
            "f789_generic_interface_target_packet": rel(IN_F789),
        },
        "why_this_packet_exists": [
            "F824 already freezes the exact required-form object and names one exact missing required-form statement field.",
            "P825 shows that neighboring statement/form slots exist, but the exact required-form statement needed by the new lane is still unexported.",
        ],
        "target_object": {
            "object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_exact_required_form_statement_target_v1",
            "goal": "Freeze the exact required-form statement object still missing for the new T213/T216 -> alpha_s schema lane.",
            "required_fields": [
                {
                    "name": "exact_statement_required_form_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F824 required-form target and not silently replace the problem.",
                },
                {
                    "name": "required_form_statement_class_candidate_support_refs",
                    "required": True,
                    "hard_limit": "Must preserve only candidate-grade neighboring statement/form-class support and must not promote it into exact discharge.",
                },
                {
                    "name": "neighboring_statement_or_form_slot_refs",
                    "required": True,
                    "hard_limit": "Must keep explicit which neighboring statement or form slots remain nonidentical to the new lane statement.",
                },
                {
                    "name": "exact_required_form_statement_ref",
                    "required": True,
                    "hard_limit": "Must state what exact required-form statement is needed for the new lane without claiming that the statement already exists.",
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny exact new-lane admission, silent reuse of neighboring statement or form slots, provider-class shift success, QCD closure, and ToE closure.",
                },
            ],
        },
        "target_refs": {
            "exact_statement_required_form_target_ref": f824_target.get("object_id"),
            "required_form_statement_class_candidate_support_refs": [
                f823_target.get("object_id"),
                f822_target.get("object_id"),
                f821_target.get("object_id"),
                f820_target.get("object_id"),
                f819_target.get("object_id"),
                f818_target.get("object_id"),
                f814_target.get("object_id"),
                f812_target.get("object_id"),
                f789_target.get("object_id"),
            ],
            "neighboring_statement_or_form_slot_refs": [
                "exact_statement_required_form_ref",
                "exact_output_schema_statement",
                "lawful_schema_domain_admission_output_schema",
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
                "selected_source_binding_output_schema",
                "normalization_rule_ref",
            ],
        },
        "current_honest_reading": [
            "The repo preserves statement and form-like structure around the new lane, but only through neighboring target fields and generic interface scaffolding.",
            "No current export yet names the exact required-form statement required by F824.",
            "F825 freezes that exact missing required-form statement object without pretending that lawful admission already exists.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current export can lawfully supply exact_required_form_statement_ref for alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_exact_required_form_statement_target_v1 without silent domain identification.",
        "hard_limits": [
            "Does not claim that exact required-form statement already exists.",
            "Does not claim that any neighboring statement or form slot silently discharges the new lane.",
            "Does not claim that the T213/T216 lane already enters the alpha_s schema domain.",
            "Does not claim that provider-class shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F825",
        "status": status,
        "as_of": AS_OF,
        "target_object_id": artifact["target_object"]["object_id"],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
