#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P823 = GENERATED / "p823_current_strict_alpha_s_output_schema_statement_class_candidate_supported_exact_output_schema_statement_blocked.json"
IN_P824 = GENERATED / "p824_current_strict_alpha_s_statement_form_class_candidate_supported_exact_statement_required_form_blocked.json"
IN_P825 = GENERATED / "p825_current_strict_alpha_s_required_form_statement_class_candidate_supported_exact_required_form_statement_blocked.json"
IN_F823 = GENERATED / "f823_current_strict_alpha_s_exact_output_schema_statement_target_packet.json"
IN_F824 = GENERATED / "f824_current_strict_alpha_s_exact_statement_required_form_target_packet.json"
IN_F825 = GENERATED / "f825_current_strict_alpha_s_exact_required_form_statement_target_packet.json"
IN_S2 = ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md"

OUT = GENERATED / "p826_current_strict_alpha_s_exact_required_form_statement_same_lane_exhaustion_boundary_audit_probe.json"
OUT_SUMMARY = GENERATED / "p826_current_strict_alpha_s_exact_required_form_statement_same_lane_exhaustion_boundary_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def has_required_field(obj: dict[str, Any], name: str) -> bool:
    return any(
        isinstance(item, dict) and item.get("name") == name
        for item in (obj.get("required_fields") or [])
    )


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P823, IN_P824, IN_P825, IN_F823, IN_F824, IN_F825, IN_S2]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P826",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p823 = load_json(IN_P823)
    p824 = load_json(IN_P824)
    p825 = load_json(IN_P825)
    f823 = load_json(IN_F823)
    f824 = load_json(IN_F824)
    f825 = load_json(IN_F825)
    s2_text = load_text(IN_S2)

    p823_theorem = p823.get("theorem_result") or {}
    p824_theorem = p824.get("theorem_result") or {}
    p825_theorem = p825.get("theorem_result") or {}
    f823_target = f823.get("target_object") or {}
    f824_target = f824.get("target_object") or {}
    f825_target = f825.get("target_object") or {}

    support_refs = (f825.get("target_refs") or {}).get("required_form_statement_class_candidate_support_refs") or []
    slot_refs = (f825.get("target_refs") or {}).get("neighboring_statement_or_form_slot_refs") or []

    checks = [
        {
            "id": "f823_exact_output_schema_statement_target_already_frozen",
            "pass": (
                f823.get("status")
                == "F823_EXECUTED_CURRENT_STRICT_ALPHA_S_EXACT_OUTPUT_SCHEMA_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f823_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_exact_output_schema_statement_target_v1"
            ),
            "details": "F823 already freezes the exact output-schema statement target for the new lane.",
        },
        {
            "id": "f824_exact_statement_required_form_target_already_frozen",
            "pass": (
                f824.get("status")
                == "F824_EXECUTED_CURRENT_STRICT_ALPHA_S_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET_NO_FALSE_PASS"
                and f824_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_exact_statement_required_form_target_v1"
                and has_required_field(f824_target, "exact_required_form_statement_ref")
            ),
            "details": "F824 already freezes the exact statement-required-form target and names exact_required_form_statement_ref as one exact missing field.",
        },
        {
            "id": "f825_exact_required_form_statement_target_already_frozen",
            "pass": (
                f825.get("status")
                == "F825_EXECUTED_CURRENT_STRICT_ALPHA_S_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f825_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_exact_required_form_statement_target_v1"
                and has_required_field(f825_target, "exact_required_form_statement_ref")
            ),
            "details": "F825 already freezes the exact required-form statement target itself.",
        },
        {
            "id": "p825_keeps_exact_required_form_statement_unexported",
            "pass": (
                p825.get("status")
                == "P825_CURRENT_STRICT_ALPHA_S_REQUIRED_FORM_STATEMENT_CLASS_CANDIDATE_SUPPORTED_EXACT_REQUIRED_FORM_STATEMENT_BLOCKED"
                and p825_theorem.get("exact_required_form_statement_exported_now") is False
                and p825_theorem.get("sharp_blocker_field") == "exact_required_form_statement_ref"
            ),
            "details": "P825 already keeps the exact required-form statement unexported on the current repo state.",
        },
        {
            "id": "f825_already_packs_neighboring_support_and_slot_context",
            "pass": (
                len(support_refs) >= 9
                and len(slot_refs) >= 8
                and has_required_field(f825_target, "required_form_statement_class_candidate_support_refs")
                and has_required_field(f825_target, "neighboring_statement_or_form_slot_refs")
            ),
            "details": "F825 already packs the neighboring statement/form-class support refs and neighboring slot refs into one exact target context.",
        },
        {
            "id": "blocker_chain_has_already_descended_to_current_floor",
            "pass": (
                p823.get("status")
                == "P823_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_STATEMENT_CLASS_CANDIDATE_SUPPORTED_EXACT_OUTPUT_SCHEMA_STATEMENT_BLOCKED"
                and p823_theorem.get("sharp_blocker_field") == "exact_output_schema_statement"
                and p824.get("status")
                == "P824_CURRENT_STRICT_ALPHA_S_STATEMENT_FORM_CLASS_CANDIDATE_SUPPORTED_EXACT_STATEMENT_REQUIRED_FORM_BLOCKED"
                and p824_theorem.get("sharp_blocker_field") == "exact_statement_required_form_ref"
                and p825.get("status")
                == "P825_CURRENT_STRICT_ALPHA_S_REQUIRED_FORM_STATEMENT_CLASS_CANDIDATE_SUPPORTED_EXACT_REQUIRED_FORM_STATEMENT_BLOCKED"
                and p825_theorem.get("sharp_blocker_field") == "exact_required_form_statement_ref"
            ),
            "details": "The blocker chain has already been pushed down from output-schema statement to statement-required-form to exact required-form statement, so the current lane is already at its present passive floor.",
        },
        {
            "id": "s2_noncyclic_continuation_discipline_applies",
            "pass": (
                "strict-core ToE closure using only strict-side sources" in s2_text
                and "new provider class and noncyclic anchor, not a repetition of L5/L12." in s2_text
            ),
            "details": "S2 still enforces noncyclic continuation discipline against repeated same-blocker passive splitting.",
        },
        {
            "id": "no_residual_passive_same_lane_loophole_below_current_blocker",
            "pass": (
                p825_theorem.get("required_form_statement_class_candidate_supported_now") is True
                and p825_theorem.get("exact_required_form_statement_exported_now") is False
                and len(support_refs) >= 9
                and len(slot_refs) >= 8
            ),
            "details": "Given that the support stack is already packed and the exact statement still remains unexported, no residual passive same-lane loophole remains below the current blocker.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "boundary_name": "CurrentStrictAlphaSExactRequiredFormStatementSameLaneExhaustionBoundary_strict_v1",
        "exact_required_form_statement_exported_on_current_repo_state": False if all_pass else None,
        "same_lane_passive_groundwork_exhausted": True if all_pass else None,
        "next_honest_move_requires_continuation_boundary_export": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P826_CURRENT_STRICT_ALPHA_S_EXACT_REQUIRED_FORM_STATEMENT_SAME_LANE_EXHAUSTION_BOUNDARY_AUDIT_PROBE"
        if all_pass
        else "P826_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P826",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p823_output_schema_statement_probe": rel(IN_P823),
            "p824_statement_required_form_probe": rel(IN_P824),
            "p825_required_form_statement_probe": rel(IN_P825),
            "f823_exact_output_schema_statement_target_packet": rel(IN_F823),
            "f824_exact_statement_required_form_target_packet": rel(IN_F824),
            "f825_exact_required_form_statement_target_packet": rel(IN_F825),
            "s2_priority_packet": rel(IN_S2),
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The exact target chain is already frozen down to the exact required-form statement object itself.",
            "The current repo state still keeps that exact required-form statement unexported, while the neighboring support and slot context is already packed inside F825.",
            "So no residual passive same-lane loophole remains below the current blocker on this lane.",
        ],
        "recommended_next_packet": {
            "id": "F826_CURRENT_STRICT_ALPHA_S_EXACT_REQUIRED_FORM_STATEMENT_CONTINUATION_BOUNDARY_PACKET",
            "goal": "Export the continuation boundary after the current same-lane passive groundwork is exhausted but the exact required-form statement still remains missing.",
            "export_object_id": "alpha_s_exact_required_form_statement_continuation_boundary_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P826",
        "status": status,
        "as_of": AS_OF,
        "exact_required_form_statement_exported_on_current_repo_state": theorem_result[
            "exact_required_form_statement_exported_on_current_repo_state"
        ],
        "same_lane_passive_groundwork_exhausted": theorem_result["same_lane_passive_groundwork_exhausted"],
        "next_honest_move_requires_continuation_boundary_export": theorem_result[
            "next_honest_move_requires_continuation_boundary_export"
        ],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
