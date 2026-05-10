#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P822 = GENERATED / "p822_current_strict_alpha_s_output_schema_class_candidate_supported_lawful_schema_domain_admission_output_schema_blocked.json"
IN_F821 = GENERATED / "f821_current_strict_alpha_s_lawful_schema_domain_admission_target_packet.json"
IN_F820 = GENERATED / "f820_current_strict_alpha_s_schema_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F819 = GENERATED / "f819_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_schema_adapter_or_carrier_identification_rule_target_packet.json"
IN_F818 = GENERATED / "f818_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_schema_interface_target_packet.json"
IN_F814 = GENERATED / "f814_current_strict_alpha_s_strict_source_shannon_source_binding_selection_or_preference_rule_schema_target_packet.json"

OUT = GENERATED / "f822_current_strict_alpha_s_lawful_schema_domain_admission_output_schema_target_packet.json"
OUT_SUMMARY = GENERATED / "f822_current_strict_alpha_s_lawful_schema_domain_admission_output_schema_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P822, IN_F821, IN_F820, IN_F819, IN_F818, IN_F814]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F822",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p822 = load_json(IN_P822)
    f821 = load_json(IN_F821)
    f820 = load_json(IN_F820)
    f819 = load_json(IN_F819)
    f818 = load_json(IN_F818)
    f814 = load_json(IN_F814)

    p822_theorem = p822.get("theorem_result") or {}
    f821_target = f821.get("target_object") or {}
    f820_target = f820.get("target_object") or {}
    f819_target = f819.get("target_object") or {}
    f818_target = f818.get("target_object") or {}
    f814_target = f814.get("target_object") or {}

    if (
        p822.get("status")
        == "P822_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_CLASS_CANDIDATE_SUPPORTED_LAWFUL_SCHEMA_DOMAIN_ADMISSION_OUTPUT_SCHEMA_BLOCKED"
        and p822_theorem.get("output_schema_class_candidate_supported_now") is True
        and p822_theorem.get("lawful_schema_domain_admission_output_schema_exported_now") is False
        and p822_theorem.get("sharp_blocker_field") == "lawful_schema_domain_admission_output_schema"
        and p822_theorem.get("next_honest_move_is_freeze_exact_output_schema_target") is True
        and f821.get("status")
        == "F821_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_SCHEMA_DOMAIN_ADMISSION_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F822_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_SCHEMA_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F822_REQUIRES_REVIEW"

    artifact = {
        "stage": "F822",
        "packet_name": "CurrentStrictAlphaSLawfulSchemaDomainAdmissionOutputSchemaTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p822_output_schema_probe": rel(IN_P822),
            "f821_lawful_admission_target_packet": rel(IN_F821),
            "f820_combined_boundary_target_packet": rel(IN_F820),
            "f819_upstream_rule_target_packet": rel(IN_F819),
            "f818_upstream_interface_target_packet": rel(IN_F818),
            "f814_downstream_schema_target_packet": rel(IN_F814),
        },
        "why_this_packet_exists": [
            "F821 already freezes the lawful schema-domain admission object and names one exact missing output field.",
            "P822 shows that output-schema class support exists on nearby lanes, but the exact lawful schema-domain admission output schema for the new lane is still unexported.",
        ],
        "target_object": {
            "object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_output_schema_target_v1",
            "goal": "Freeze the exact lawful schema-domain admission output-schema object still missing for the new T213/T216 -> alpha_s schema lane.",
            "required_fields": [
                {
                    "name": "lawful_schema_domain_admission_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F821 lawful-admission target and not silently replace the problem.",
                },
                {
                    "name": "output_schema_class_candidate_support_refs",
                    "required": True,
                    "hard_limit": "Must preserve only candidate-grade output-schema class support from neighboring targets and must not promote it into exact discharge.",
                },
                {
                    "name": "upstream_rule_or_interface_output_refs",
                    "required": True,
                    "hard_limit": "Must keep explicit the upstream output-schema classes that remain nonidentical to the new lane output schema.",
                },
                {
                    "name": "downstream_schema_output_ref",
                    "required": True,
                    "hard_limit": "Must keep explicit the downstream schema-output class without silently reusing it as the new lawful-admission output schema.",
                },
                {
                    "name": "exact_output_schema_statement",
                    "required": True,
                    "hard_limit": "Must state what exact lawful-admission output would be required for the new lane without claiming that such output already exists.",
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny exact new-lane admission, silent reuse of neighboring output-schema classes, provider-class shift success, QCD closure, and ToE closure.",
                },
            ],
        },
        "target_refs": {
            "lawful_schema_domain_admission_target_ref": f821_target.get("object_id"),
            "output_schema_class_candidate_support_refs": [
                f820_target.get("object_id"),
                f819_target.get("object_id"),
                f818_target.get("object_id"),
                f814_target.get("object_id"),
            ],
            "upstream_rule_or_interface_output_refs": [
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
            ],
            "downstream_schema_output_ref": "selected_source_binding_output_schema",
        },
        "current_honest_reading": [
            "The repo preserves output-schema class around the new lane, but only at neighboring target level.",
            "No current export yet names the exact lawful schema-domain admission output schema required by F821.",
            "F822 freezes that exact missing output-schema object without pretending that lawful admission already exists.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current export can lawfully supply exact_output_schema_statement for alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_output_schema_target_v1 without silent domain identification.",
        "hard_limits": [
            "Does not claim that exact lawful schema-domain admission output schema already exists.",
            "Does not claim that any upstream or downstream output-schema class silently discharges the new lane.",
            "Does not claim that the T213/T216 lane already enters the alpha_s schema domain.",
            "Does not claim that provider-class shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F822",
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
