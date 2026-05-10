#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F817 = GENERATED / "f817_current_strict_alpha_s_different_selection_provider_class_shift_candidate_reference_packet.json"
IN_F789_SUMMARY = GENERATED / "f789_current_strict_alpha_s_normalized_boundary_interface_target_packet_summary.json"
IN_P764_SUMMARY = GENERATED / "p764_current_strict_t218_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_probe_summary.json"
IN_P765_SUMMARY = GENERATED / "p765_current_strict_t219_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_nonexport_audit_probe_summary.json"
IN_P766_SUMMARY = GENERATED / "p766_current_strict_t220_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_probe_summary.json"
IN_P809_SUMMARY = GENERATED / "p809_current_strict_alpha_s_strict_source_shannon_shift_to_alpha_s_domain_interface_target_audit_probe_summary.json"
IN_P810_SUMMARY = GENERATED / "p810_current_strict_alpha_s_strict_source_shannon_shift_to_domain_carrier_identification_or_adapter_rule_audit_probe_summary.json"

OUT = GENERATED / "p818_current_strict_alpha_s_different_provider_class_shift_interface_target_first_continuation_class_audit_probe.json"
OUT_SUMMARY = GENERATED / "p818_current_strict_alpha_s_different_provider_class_shift_interface_target_first_continuation_class_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [
        IN_F817,
        IN_F789_SUMMARY,
        IN_P764_SUMMARY,
        IN_P765_SUMMARY,
        IN_P766_SUMMARY,
        IN_P809_SUMMARY,
        IN_P810_SUMMARY,
    ]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P818",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    f817 = load_json(IN_F817)
    f789_summary = load_json(IN_F789_SUMMARY)
    p764 = load_json(IN_P764_SUMMARY)
    p765 = load_json(IN_P765_SUMMARY)
    p766 = load_json(IN_P766_SUMMARY)
    p809 = load_json(IN_P809_SUMMARY)
    p810 = load_json(IN_P810_SUMMARY)

    f817_export = f817.get("exported_object") or {}

    checks = [
        {
            "id": "f817_leaves_the_blocker_exactly_at_alpha_s_shift_interface_level",
            "pass": (
                f817.get("status")
                == "F817_EXECUTED_CURRENT_STRICT_ALPHA_S_DIFFERENT_SELECTION_PROVIDER_CLASS_SHIFT_CANDIDATE_REFERENCE_PACKET_NO_FALSE_PASS"
                and f817_export.get("candidate_reference_lane_grade") == "reference_context_candidate_only"
                and f817_export.get("alpha_s_shift_interface_status") == "blocked_nonexport"
                and f817_export.get("provider_class_shift_realization_status") == "not_realized"
            ),
            "details": "F817 leaves the new provider-class lane admitted only as a candidate, with the blocker localized at the alpha_s shift-interface level.",
        },
        {
            "id": "generic_alpha_s_bridge_history_uses_interface_target_first",
            "pass": (
                f789_summary.get("status")
                == "F789_EXECUTED_CURRENT_STRICT_ALPHA_S_NORMALIZED_BOUNDARY_INTERFACE_TARGET_PACKET_NO_FALSE_PASS"
            ),
            "details": "The generic alpha_s bridge lane already uses an interface-target-first move when re-entry is blocked.",
        },
        {
            "id": "earlier_provider_shift_precedent_uses_interface_target_before_adapter_rule",
            "pass": (
                p809.get("status")
                == "P809_CURRENT_STRICT_ALPHA_S_NO_EXACT_STRICT_SOURCE_SHANNON_SHIFT_TO_ALPHA_S_DOMAIN_INTERFACE_TARGET_EXPORTED_TARGET_FREEZE_REQUIRED"
                and p809.get("next_honest_move_requires_freezing_exact_shift_to_domain_interface_target")
                is True
                and p810.get("status")
                == "P810_CURRENT_STRICT_ALPHA_S_NO_CARRIER_IDENTIFICATION_OR_ADAPTER_RULE_FOR_STRICT_SOURCE_SHANNON_SHIFT_TO_DOMAIN_INTERFACE_TARGET_EXPORTED_TARGET_FREEZE_REQUIRED"
                and p810.get("exact_shift_to_domain_interface_target_frozen") is True
                and p810.get("next_honest_move_requires_freezing_exact_carrier_identification_or_adapter_rule_target")
                is True
            ),
            "details": "The earlier alpha_s provider-shift precedent clearly orders interface-target freeze before adapter-rule targeting.",
        },
        {
            "id": "t213_t216_own_lane_also_uses_interface_target_before_actual_interface_attempt",
            "pass": (
                p764.get("status")
                == "PASS_STRICT_T218_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_TARGET_EXPORTED"
                and p764.get("t218_target_exported_on_current_repo_state") is True
                and p765.get("status")
                == "PASS_STRICT_T219_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
                and p765.get("next_honest_move_is_actual_t218_interface_realization_attempt_or_attempt_level_failure_boundary")
                is True
                and p766.get("status")
                == "PASS_STRICT_T220_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
                and p766.get("t220_attempt_exported_on_current_repo_state") is True
            ),
            "details": "The T213/T216 lane itself progresses through exact interface-target freeze before any actual interface attempt, so its own discipline is interface-target-first as well.",
        },
        {
            "id": "adapter_or_carrier_target_first_would_skip_the_localized_blocker_layer",
            "pass": (
                f817_export.get("alpha_s_shift_interface_status") == "blocked_nonexport"
                and p809.get("next_honest_move_requires_freezing_exact_shift_to_domain_interface_target")
                is True
            ),
            "details": "Because F817 localizes the blocker at the shift-interface layer, jumping directly to adapter/carrier targeting would skip the still-unfrozen interface object class.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "interface_target_first_continuation_class_admitted": all_pass,
        "immediate_next_move_class": "freeze_exact_alpha_s_side_shift_interface_target" if all_pass else None,
        "adapter_rule_target_first": False if all_pass else None,
        "carrier_identification_target_first": False if all_pass else None,
        "next_honest_move_requires_interface_target_oriented_probe_or_freeze": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P818_CURRENT_STRICT_ALPHA_S_DIFFERENT_PROVIDER_CLASS_SHIFT_INTERFACE_TARGET_FIRST_CONTINUATION_CLASS_AUDITED"
        if all_pass
        else "P818_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P818",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f817_different_provider_class_candidate_reference_packet": rel(IN_F817),
            "f789_generic_alpha_s_interface_target_summary": rel(IN_F789_SUMMARY),
            "p764_t218_interface_target_summary": rel(IN_P764_SUMMARY),
            "p765_t219_interface_nonexport_boundary_summary": rel(IN_P765_SUMMARY),
            "p766_t220_interface_attempt_summary": rel(IN_P766_SUMMARY),
            "p809_shannon_shift_interface_target_precedent_summary": rel(IN_P809_SUMMARY),
            "p810_shannon_adapter_rule_after_interface_precedent_summary": rel(IN_P810_SUMMARY),
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "F817 leaves the new provider-class route blocked exactly at the alpha_s shift-interface layer.",
            "Both the generic alpha_s bridge history and the earlier provider-shift precedent use interface-target-first continuation before any adapter-rule targeting.",
            "The T213/T216 lane itself also uses exact interface-target freeze before actual interface attempts.",
            "So the immediate honest move class after F817 is interface-target-first, not adapter-rule-first.",
        ],
        "recommended_next_packet": {
            "id": "F818_CURRENT_STRICT_ALPHA_S_DIFFERENT_PROVIDER_CLASS_SHIFT_INTERFACE_TARGET_FIRST_CONTINUATION_BOUNDARY_PACKET",
            "goal": "Export the continuation boundary that freezes interface-target-first ordering after F817.",
            "export_object_id": "alpha_s_different_provider_class_shift_interface_target_first_continuation_boundary_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P818",
        "status": status,
        "as_of": AS_OF,
        "interface_target_first_continuation_class_admitted": theorem_result[
            "interface_target_first_continuation_class_admitted"
        ],
        "immediate_next_move_class": theorem_result["immediate_next_move_class"],
        "adapter_rule_target_first": theorem_result["adapter_rule_target_first"],
        "carrier_identification_target_first": theorem_result["carrier_identification_target_first"],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
