#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P802 = GENERATED / "p802_current_strict_alpha_s_provider_class_reorganization_audit_probe.json"
IN_F802 = GENERATED / "f802_current_strict_alpha_s_provider_action_rule_target_packet.json"
IN_F801 = GENERATED / "f801_current_strict_alpha_s_same_domain_provider_skeleton_packet.json"
IN_P792 = GENERATED / "p792_current_strict_alpha_s_family_selection_order_rule_probe.json"
IN_P797 = GENERATED / "p797_current_strict_alpha_s_invariant_extremum_vs_reference_scale_audit_probe.json"
IN_P794 = GENERATED / "p794_current_strict_alpha_s_normalization_boundary_subclause_audit_probe.json"
IN_P710 = ROOT / "P710_CURRENT_NONSTRICT_PROXY_TO_GEV_CALIBRATION_MAP_FROM_F704_EIGENSPECTRUM_PROBE.md"
IN_POLICY = ROOT / "external_data" / "proxy_to_gev_calibration_policy_v1.json"

OUT = GENERATED / "p803_current_strict_alpha_s_same_domain_relation_bundle_admission_probe.json"
OUT_SUMMARY = GENERATED / "p803_current_strict_alpha_s_same_domain_relation_bundle_admission_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


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


def contains_all(text: str, needles: list[str]) -> bool:
    hay = normalize(text)
    return all(normalize(needle) in hay for needle in needles)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [
        IN_P802,
        IN_F802,
        IN_F801,
        IN_P792,
        IN_P797,
        IN_P794,
        IN_P710,
        IN_POLICY,
    ]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P803",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p802 = load_json(IN_P802)
    f802 = load_json(IN_F802)
    f801 = load_json(IN_F801)
    p792 = load_json(IN_P792)
    p797 = load_json(IN_P797)
    p794 = load_json(IN_P794)
    policy = load_json(IN_POLICY)
    p710_text = IN_P710.read_text(encoding="utf-8")

    skeleton = f801.get("exported_object") or {}
    action_target = f802.get("target_object") or {}
    p792_rule = p792.get("probe_local_order_rule_tuple") or {}
    p797_split = p797.get("clause_split") or {}
    p794_split = (p794.get("subclause_split_audit") or {})

    checks = [
        {
            "id": "f802_already_isolates_missing_provider_action_layer",
            "pass": (
                p802.get("status")
                == "P802_CURRENT_STRICT_ALPHA_S_PASSIVE_PROVIDER_SKELETON_SUPPORTED_ACTIVE_PROVIDER_ACTION_RULE_BLOCKED"
                and f802.get("status")
                == "F802_EXECUTED_CURRENT_STRICT_ALPHA_S_PROVIDER_ACTION_RULE_TARGET_PACKET_NO_FALSE_PASS"
                and skeleton.get("object_id") == "alpha_s_reference_scale_same_domain_provider_skeleton_v1"
                and action_target.get("object_id") == "alpha_s_reference_scale_provider_action_rule_target_v1"
            ),
            "details": "P802/F802 already isolate the missing provider-action layer above the passive same-domain skeleton.",
        },
        {
            "id": "family_winner_relation_present_on_same_domain",
            "pass": (
                p792.get("status")
                == "P792_CURRENT_STRICT_ALPHA_S_PROBE_LOCAL_ORDER_RULE_UNIQUE_WINNER_NONEXPORT"
                and p792.get("unique_winner_exists") is True
                and p792.get("unique_winner_family_id") == "f704_max_mode_anchor_family"
                and p792_rule.get("export_status") == "probe_local_only"
            ),
            "details": "P792 exports one same-domain winner relation on the finite alpha_s family domain, but only at probe-local grade.",
        },
        {
            "id": "extremum_relation_present_on_same_domain",
            "pass": (
                p797.get("status")
                == "P797_CURRENT_STRICT_ALPHA_S_INVARIANT_EXTREMUM_CANDIDATE_SUPPORTED_REFERENCE_SCALE_RULE_BLOCKED"
                and p797_split.get("invariant_extremum_clause_status")
                == "candidate_supported_not_yet_exported"
                and p797_split.get("sharp_blocker_clause") == "reference_scale_rule"
            ),
            "details": "P797 exports one same-domain extremum relation: unique stable maximum support below any reference-scale rule.",
        },
        {
            "id": "bounded_boundary_relation_present_on_same_domain",
            "pass": (
                p794.get("status")
                == "P794_CURRENT_STRICT_ALPHA_S_BOUNDED_GRID_CANDIDATE_SUPPORTED_TOP_BOUNDARY_ANCHOR_BLOCKED"
                and p794_split.get("bounded_grid_clause_status")
                == "candidate_supported_not_yet_exported"
                and p794.get("derived_arithmetic_snapshot", {}).get("top_boundary_one_present") is True
            ),
            "details": "P794 exports one same-domain bounded-boundary relation: max-normalized arithmetic yields a bounded grid with point 1 present.",
        },
        {
            "id": "relation_bundle_still_stops_before_active_provider_action_rule",
            "pass": any(
                isinstance(item, dict) and item.get("name") == "provider_action_rule_ref"
                for item in (action_target.get("required_fields") or [])
            ),
            "details": "F802 itself records the missing active rule above these relations, so the relation bundle does not yet act as a provider.",
        },
        {
            "id": "nonstrict_calibration_lane_explicitly_excluded",
            "pass": (
                policy.get("scope") == "nonstrict_physical_units_calibration_only"
                and contains_all(
                    p710_text,
                    [
                        "non-strict calibration map",
                        "no strict physical-unit map",
                    ],
                )
            ),
            "details": "The proxy-to-GeV calibration lane remains explicitly excluded from the same-domain relation-bundle question.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]

    clause_split = {
        "relation_bundle_clause_status": "export_admitted"
        if checks[1]["pass"] and checks[2]["pass"] and checks[3]["pass"] and checks[5]["pass"]
        else "requires_review",
        "provider_action_rule_clause_status": "blocked_nonexport"
        if checks[0]["pass"] and checks[4]["pass"]
        else "requires_review",
        "sharp_blocker_clause": "provider_action_rule_ref"
        if all(item["pass"] for item in checks)
        else None,
    }

    if all(item["pass"] for item in checks):
        status = "P803_CURRENT_STRICT_ALPHA_S_SAME_DOMAIN_RELATION_BUNDLE_EXPORT_ADMITTED_PROVIDER_ACTION_RULE_STILL_BLOCKED"
    else:
        status = "P803_REQUIRES_REVIEW"

    artifact = {
        "stage": "P803",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p802_provider_action_probe": rel(IN_P802),
            "f802_provider_action_target": rel(IN_F802),
            "f801_provider_skeleton_packet": rel(IN_F801),
            "p792_family_selection_relation": rel(IN_P792),
            "p797_extremum_relation": rel(IN_P797),
            "p794_bounded_boundary_relation": rel(IN_P794),
            "p710_nonstrict_proxy_to_gev_probe": rel(IN_P710),
            "proxy_to_gev_policy": rel(IN_POLICY),
        },
        "same_domain_relation_bundle_candidate": {
            "provider_skeleton_ref": skeleton.get("object_id"),
            "provider_action_target_ref": action_target.get("object_id"),
            "relation_refs": [
                rel(IN_P792),
                rel(IN_P797),
                rel(IN_P794),
            ],
            "relation_snapshot": {
                "family_winner_id": p792.get("unique_winner_family_id"),
                "extremum_max_m_proxy": p797.get("extremum_snapshot", {}).get("max_m_proxy"),
                "bounded_grid_top_boundary_one_present": p794.get("derived_arithmetic_snapshot", {}).get("top_boundary_one_present"),
            },
        },
        "clause_split": clause_split,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The current same-domain alpha_s lane already exports a relation bundle beneath the missing provider action rule.",
            "That bundle includes one family-winner relation, one unique-extremum relation, and one bounded-boundary relation on the same lane.",
            "What remains missing is the active rule that would make this relation bundle act as a provider.",
        ],
        "recommended_next_packet": {
            "id": "F803_CURRENT_STRICT_ALPHA_S_SAME_DOMAIN_RELATION_BUNDLE_PACKET",
            "goal": "Export the admitted same-domain relation bundle explicitly while keeping the provider action rule unresolved.",
            "export_object_id": "alpha_s_reference_scale_same_domain_relation_bundle_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P803",
        "status": status,
        "as_of": AS_OF,
        "relation_bundle_clause_status": clause_split["relation_bundle_clause_status"],
        "provider_action_rule_clause_status": clause_split["provider_action_rule_clause_status"],
        "sharp_blocker_clause": clause_split["sharp_blocker_clause"],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
