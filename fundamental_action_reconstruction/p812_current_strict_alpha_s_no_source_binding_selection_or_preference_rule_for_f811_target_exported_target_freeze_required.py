#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F811 = GENERATED / "f811_current_strict_alpha_s_strict_source_shannon_shift_to_domain_source_binding_target_packet.json"
IN_P811 = GENERATED / "p811_current_strict_alpha_s_source_support_present_exact_source_binding_unexported_adapter_action_schema_blocked_source_binding_target_freeze_required.json"
IN_P791 = GENERATED / "p791_current_strict_alpha_s_selection_principle_reuse_audit_probe.json"
IN_P792 = GENERATED / "p792_current_strict_alpha_s_family_selection_order_rule_probe.json"
IN_F792 = GENERATED / "f792_current_strict_alpha_s_family_selection_order_rule_target_packet.json"
IN_F447 = ROOT / "F447_CURRENT_STRICT_T169_QW2122_SCALAR_TO_T168_PER_SITE_VALUE_PROVIDER_POWERLAW_ELEMENT_ORDER_REFERENCE_PACKET.md"
IN_T169 = ROOT / "T169_CURRENT_STRICT_CONSTRAINED_LIFT_FROM_QW2122_SCALAR_TO_T168_PER_SITE_PROVIDER_TARGET_SPEC.md"

OUT = GENERATED / "p812_current_strict_alpha_s_no_source_binding_selection_or_preference_rule_for_f811_target_exported_target_freeze_required.json"
OUT_SUMMARY = GENERATED / "p812_current_strict_alpha_s_no_source_binding_selection_or_preference_rule_for_f811_target_exported_target_freeze_required_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def normalize(text: str) -> str:
    return " ".join(
        text.lower()
        .replace("“", '"')
        .replace("”", '"')
        .replace("’", "'")
        .replace("‑", "-")
        .replace("–", "-")
        .replace("—", "-")
        .replace("->", " ")
        .replace("→", " ")
        .replace("/", " ")
        .replace("-", " ")
        .replace("_", " ")
        .split()
    )


def find_exact_rule_hits(target_id: str, rule_terms: list[str]) -> list[str]:
    hits: list[str] = []
    excluded_suffixes = [
        "generated/f811_current_strict_alpha_s_strict_source_shannon_shift_to_domain_source_binding_target_packet.json",
        "generated/f811_current_strict_alpha_s_strict_source_shannon_shift_to_domain_source_binding_target_packet_summary.json",
        "generated/p811_current_strict_alpha_s_source_support_present_exact_source_binding_unexported_adapter_action_schema_blocked_source_binding_target_freeze_required.json",
        "generated/p811_current_strict_alpha_s_source_support_present_exact_source_binding_unexported_adapter_action_schema_blocked_source_binding_target_freeze_required_summary.json",
        "F811_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_DOMAIN_SOURCE_BINDING_TARGET_PACKET.md",
        "P811_CURRENT_STRICT_ALPHA_S_SOURCE_SUPPORT_PRESENT_EXACT_SOURCE_BINDING_UNEXPORTED_ADAPTER_ACTION_SCHEMA_BLOCKED_SOURCE_BINDING_TARGET_FREEZE_REQUIRED.md",
        "f811_current_strict_alpha_s_strict_source_shannon_shift_to_domain_source_binding_target_packet.py",
        "p811_current_strict_alpha_s_source_support_present_exact_source_binding_unexported_adapter_action_schema_blocked_source_binding_target_freeze_required.py",
        "generated/p812_current_strict_alpha_s_no_source_binding_selection_or_preference_rule_for_f811_target_exported_target_freeze_required.json",
        "generated/p812_current_strict_alpha_s_no_source_binding_selection_or_preference_rule_for_f811_target_exported_target_freeze_required_summary.json",
        "generated/f812_current_strict_alpha_s_strict_source_shannon_shift_to_domain_source_binding_selection_or_preference_rule_target_packet.json",
        "generated/f812_current_strict_alpha_s_strict_source_shannon_shift_to_domain_source_binding_selection_or_preference_rule_target_packet_summary.json",
        "P812_CURRENT_STRICT_ALPHA_S_NO_SOURCE_BINDING_SELECTION_OR_PREFERENCE_RULE_FOR_F811_TARGET_EXPORTED_TARGET_FREEZE_REQUIRED.md",
        "F812_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_DOMAIN_SOURCE_BINDING_SELECTION_OR_PREFERENCE_RULE_TARGET_PACKET.md",
        "p812_current_strict_alpha_s_no_source_binding_selection_or_preference_rule_for_f811_target_exported_target_freeze_required.py",
        "f812_current_strict_alpha_s_strict_source_shannon_shift_to_domain_source_binding_selection_or_preference_rule_target_packet.py",
    ]
    norm_target = normalize(target_id)
    norm_terms = [normalize(term) for term in rule_terms]
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        if path.suffix not in {".md", ".json", ".py"}:
            continue
        relpath = rel(path)
        if any(relpath.endswith(excluded) for excluded in excluded_suffixes):
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        hay = normalize(text)
        if norm_target in hay and any(term in hay for term in norm_terms):
            hits.append(relpath)
    return sorted(hits)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F811, IN_P811, IN_P791, IN_P792, IN_F792, IN_F447, IN_T169]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P812",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    f811 = load_json(IN_F811)
    p811 = load_json(IN_P811)
    p791 = load_json(IN_P791)
    p792 = load_json(IN_P792)
    f792 = load_json(IN_F792)
    f447_text = load_text(IN_F447)
    t169_text = load_text(IN_T169)

    f811_target = f811.get("target_object") or {}
    f811_refs = f811.get("target_refs") or {}
    p811_support = p811.get("support_objects") or {}

    exact_rule_hits = find_exact_rule_hits(
        "alpha_s_strict_source_shannon_shift_to_domain_source_binding_target_v1",
        [
            "selection rule",
            "preference rule",
            "source authority rule",
            "selection objective",
            "order rule",
        ],
    )

    f447_norm = normalize(f447_text)
    t169_norm = normalize(t169_text)
    p792_rule = (p792.get("probe_local_order_rule_tuple") or {})

    checks = [
        {
            "id": "f811_freezes_exact_source_binding_target_and_requires_selection_rule",
            "pass": (
                f811.get("status")
                == "F811_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_DOMAIN_SOURCE_BINDING_TARGET_PACKET_NO_FALSE_PASS"
                and f811_target.get("object_id")
                == "alpha_s_strict_source_shannon_shift_to_domain_source_binding_target_v1"
                and any(
                    isinstance(item, dict) and item.get("name") == "source_binding_selection_or_preference_rule_ref"
                    for item in (f811_target.get("required_fields") or [])
                )
            ),
            "details": "F811 already freezes the exact source-binding target and names the still-missing source_binding_selection_or_preference_rule_ref slot.",
        },
        {
            "id": "source_side_support_objects_are_real_but_unbound",
            "pass": (
                p811.get("status")
                == "P811_CURRENT_STRICT_ALPHA_S_SOURCE_SUPPORT_PRESENT_EXACT_SOURCE_BINDING_UNEXPORTED_ADAPTER_ACTION_SCHEMA_BLOCKED_SOURCE_BINDING_TARGET_FREEZE_REQUIRED"
                and (p811.get("theorem_result") or {}).get("current_source_candidate_support_present") is True
                and (p811.get("theorem_result") or {}).get("lawful_future_entry_object_support_present") is True
                and (p811.get("theorem_result") or {}).get("exact_source_candidate_lane_or_entry_ref_exported_for_f810_target")
                is False
            ),
            "details": "Current source-side support objects are real, but no exact source object is yet bound to the F811 route.",
        },
        {
            "id": "f447_t169_selection_template_is_foreign_domain_only",
            "pass": (
                "selection principle" in t169_norm
                and "per site provider" in t169_norm
                and "qw 2122 scalar" in f447_norm
                and p791.get("reusable_selection_principle_found") is False
            ),
            "details": "The repo knows a real strict selection template elsewhere, but only on the foreign QW-2122 -> T168 provider domain.",
        },
        {
            "id": "p792_f792_order_rule_is_alpha_s_family_specific_not_source_binding_rule",
            "pass": (
                p792.get("status")
                == "P792_CURRENT_STRICT_ALPHA_S_PROBE_LOCAL_ORDER_RULE_UNIQUE_WINNER_NONEXPORT"
                and p792_rule.get("export_status") == "probe_local_only"
                and f792.get("status")
                == "F792_EXECUTED_CURRENT_STRICT_ALPHA_S_FAMILY_SELECTION_ORDER_RULE_TARGET_PACKET_NO_FALSE_PASS"
                and "candidate_family_domain_ref"
                in {item.get("name") for item in (f792.get("target_object") or {}).get("required_fields", []) if isinstance(item, dict)}
            ),
            "details": "The only close alpha_s-side order rule today is the family-selection rule on the P789 domain; it is probe-local and acts on a different domain than F811 source binding.",
        },
        {
            "id": "no_current_export_names_exact_source_binding_selection_or_preference_rule_for_f811",
            "pass": exact_rule_hits == [],
            "details": "No current export names one exact source-binding selection or preference rule for the frozen F811 target.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "exact_source_binding_target_frozen": checks[0]["pass"],
        "source_side_support_present_but_unbound": checks[1]["pass"],
        "foreign_domain_selection_template_reusable": False if checks[2]["pass"] else None,
        "alpha_s_family_order_rule_reusable_for_source_binding": False if checks[3]["pass"] else None,
        "exact_source_binding_selection_or_preference_rule_exported_for_f811_target": False if all_pass else None,
        "next_honest_move_requires_freezing_exact_source_binding_selection_or_preference_rule_target": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P812_CURRENT_STRICT_ALPHA_S_NO_SOURCE_BINDING_SELECTION_OR_PREFERENCE_RULE_FOR_F811_TARGET_EXPORTED_TARGET_FREEZE_REQUIRED"
        if all_pass
        else "P812_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P812",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f811_source_binding_target_packet": rel(IN_F811),
            "p811_source_binding_audit_probe": rel(IN_P811),
            "p791_selection_reuse_audit": rel(IN_P791),
            "p792_family_order_rule_probe": rel(IN_P792),
            "f792_family_order_rule_target_packet": rel(IN_F792),
            "f447_foreign_selection_template": rel(IN_F447),
            "t169_foreign_selection_template_spec": rel(IN_T169),
        },
        "exact_missing_rule_target_candidate": {
            "candidate_id": "alpha_s_strict_source_shannon_shift_to_domain_source_binding_selection_or_preference_rule_target_missing_v1",
            "source_binding_target_ref": f811_target.get("object_id"),
            "current_source_candidate_support_ref": p811_support.get("current_source_candidate_support_ref"),
            "lawful_future_entry_object_support_ref": p811_support.get("lawful_future_entry_object_support_ref"),
            "repo_scan_hits_for_exact_rule": exact_rule_hits,
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The repo now freezes the exact source-binding target and already has source-side supports on hand.",
            "But its real selection templates are foreign-domain, and its only close alpha_s-side order rule acts on the family domain rather than the source-binding domain.",
            "So the next honest move is to freeze the exact missing source-binding selection/preference-rule target, not to reuse an existing rule by analogy.",
        ],
        "recommended_next_packet": {
            "id": "F812_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_DOMAIN_SOURCE_BINDING_SELECTION_OR_PREFERENCE_RULE_TARGET_PACKET",
            "goal": "Freeze the exact missing source-binding selection/preference-rule target required before one source can be lawfully chosen for the F811 route.",
            "export_object_id": "alpha_s_strict_source_shannon_shift_to_domain_source_binding_selection_or_preference_rule_target_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P812",
        "status": status,
        "as_of": AS_OF,
        "exact_source_binding_target_frozen": theorem_result["exact_source_binding_target_frozen"],
        "exact_source_binding_selection_or_preference_rule_exported_for_f811_target": theorem_result[
            "exact_source_binding_selection_or_preference_rule_exported_for_f811_target"
        ],
        "next_honest_move_requires_freezing_exact_source_binding_selection_or_preference_rule_target": theorem_result[
            "next_honest_move_requires_freezing_exact_source_binding_selection_or_preference_rule_target"
        ],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
