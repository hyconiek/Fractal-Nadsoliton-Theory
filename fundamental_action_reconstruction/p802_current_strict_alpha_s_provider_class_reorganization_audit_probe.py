#!/usr/bin/env python3
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P801 = GENERATED / "p801_current_strict_alpha_s_same_domain_provider_skeleton_admission_probe.json"
IN_F801 = GENERATED / "f801_current_strict_alpha_s_same_domain_provider_skeleton_packet.json"
IN_F800 = GENERATED / "f800_current_strict_alpha_s_reference_scale_provider_class_target_packet.json"
IN_F704 = ROOT / "F704_CURRENT_STRICT_INVARIANT_MASS_OBSERVABLE_FROM_DIAGONAL_LOCAL_PSI_HESSIAN_EIGENSYSTEM_EXPORT_PACKET.md"
IN_N705 = ROOT / "N705_CURRENT_FIRST_STRICT_PROJECTIVE_OPERATIONAL_TOE_OS_SUPPORT_WITH_INVARIANT_MASS_OBSERVABLE_THEOREM.md"
IN_N706 = ROOT / "N706_CURRENT_STRICT_RELEASE_7_OS_RESIDUAL_SIGN_GAUGE_IRRELEVANCE_PACKAGE_THEOREM.md"
IN_N703 = ROOT / "N703_CURRENT_STRICT_QUADRATIC_MASS_PROXY_MEANING_DEFINITION_THEOREM.md"
IN_P696 = ROOT / "P696_CURRENT_STRICT_PHYSICAL_COMPUTABILITY_SELECTOR_ALIGNED_CHANNEL_SPECTRUM_PROXY_FROM_PROJECTIVE_SELECTOR_CLOSURE_PROBE.md"
IN_P709 = ROOT / "P709_CURRENT_STRICT_RELEASE_7_OS_RESIDUAL_SIGN_GAUGE_IRRELEVANCE_AUDIT_PROBE.md"
IN_P710 = ROOT / "P710_CURRENT_NONSTRICT_PROXY_TO_GEV_CALIBRATION_MAP_FROM_F704_EIGENSPECTRUM_PROBE.md"
IN_POLICY = ROOT / "external_data" / "proxy_to_gev_calibration_policy_v1.json"

OUT = GENERATED / "p802_current_strict_alpha_s_provider_class_reorganization_audit_probe.json"
OUT_SUMMARY = GENERATED / "p802_current_strict_alpha_s_provider_class_reorganization_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def normalize(text: str) -> str:
    lowered = (
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
    )
    return " ".join(re.sub(r"[^a-z0-9]+", " ", lowered).split())


def contains_all(text: str, needles: list[str]) -> bool:
    hay = normalize(text)
    return all(normalize(needle) in hay for needle in needles)


def contains_any(text: str, needles: list[str]) -> bool:
    hay = normalize(text)
    return any(normalize(needle) in hay for needle in needles)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [
        IN_P801,
        IN_F801,
        IN_F800,
        IN_F704,
        IN_N705,
        IN_N706,
        IN_N703,
        IN_P696,
        IN_P709,
        IN_P710,
        IN_POLICY,
    ]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P802",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p801 = load_json(IN_P801)
    f801 = load_json(IN_F801)
    f800 = load_json(IN_F800)
    policy = load_json(IN_POLICY)

    f704_text = IN_F704.read_text(encoding="utf-8")
    n705_text = IN_N705.read_text(encoding="utf-8")
    n706_text = IN_N706.read_text(encoding="utf-8")
    n703_text = IN_N703.read_text(encoding="utf-8")
    p696_text = IN_P696.read_text(encoding="utf-8")
    p709_text = IN_P709.read_text(encoding="utf-8")
    p710_text = IN_P710.read_text(encoding="utf-8")

    same_domain_text = "\n".join([f704_text, n705_text, n706_text, n703_text, p696_text, p709_text])
    skeleton = f801.get("exported_object") or {}
    provider_target = f800.get("target_object") or {}

    checks = [
        {
            "id": "f801_already_exports_passive_same_domain_skeleton",
            "pass": (
                p801.get("status")
                == "P801_CURRENT_STRICT_ALPHA_S_SAME_DOMAIN_PROVIDER_SKELETON_EXPORT_ADMITTED_PROVIDER_CLASS_STILL_BLOCKED"
                and f801.get("status")
                == "F801_EXECUTED_CURRENT_STRICT_ALPHA_S_SAME_DOMAIN_PROVIDER_SKELETON_PACKET_NO_FALSE_PASS"
                and skeleton.get("object_id") == "alpha_s_reference_scale_same_domain_provider_skeleton_v1"
                and skeleton.get("provider_class_target_ref") == "alpha_s_reference_scale_provider_class_target_v1"
            ),
            "details": "F801 already exports the passive same-domain provider skeleton beneath the missing provider class.",
        },
        {
            "id": "same_domain_lane_supplies_passive_roles_under_future_provider_class",
            "pass": (
                contains_all(f704_text, ["mass like", "basis invariant", "dimensionless"])
                and contains_all(n705_text, ["tracked component", "basis invariant mass observable object"])
                and contains_all(n706_text, ["gauge convention layer", "exported os observable values"])
                and contains_all(n703_text, ["meaning definition", "dimensionless quadratic coefficients"])
                and contains_all(p696_text, ["physical computability", "basis level operational proxy"])
                and contains_all(p709_text, ["audit only", "invariants are unchanged"])
            ),
            "details": "The same-domain lane already supplies passive roles: observable carrier, tracked packaging, gauge safety, meaning discipline, and computability/audit support.",
        },
        {
            "id": "same_domain_lane_exports_no_active_provider_action_rule_or_supply_schema",
            "pass": not contains_any(
                same_domain_text,
                [
                    "provider class",
                    "provider action rule",
                    "supply schema",
                    "semantic principle",
                    "reference scale rule",
                    "reference datum",
                ],
            ),
            "details": "The same-domain lane still exports no active provider action rule or semantic-principle supply schema.",
        },
        {
            "id": "provider_target_itself_demands_active_rule_beyond_passive_skeleton",
            "pass": (
                provider_target.get("object_id") == "alpha_s_reference_scale_provider_class_target_v1"
                and any(
                    isinstance(item, dict) and item.get("name") == "provider_class_ref"
                    for item in (provider_target.get("required_fields") or [])
                )
                and any(
                    isinstance(item, dict) and item.get("name") == "selected_semantic_principle_supply_schema"
                    for item in (provider_target.get("required_fields") or [])
                )
            ),
            "details": "F800 itself makes clear that the missing object is active: it requires a provider-class ref and a semantic-principle supply schema.",
        },
        {
            "id": "foreign_domain_and_nonstrict_imports_remain_excluded",
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
            "details": "Foreign-domain reference imports and nonstrict calibration remain excluded from the provider-class reorganization question.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]

    clause_split = {
        "passive_skeleton_clause_status": "supported_exported"
        if checks[0]["pass"] and checks[1]["pass"]
        else "requires_review",
        "active_provider_action_rule_clause_status": "blocked_nonexport"
        if checks[2]["pass"] and checks[3]["pass"] and checks[4]["pass"]
        else "requires_review",
        "sharp_blocker_clause": "provider_action_rule_ref"
        if all(item["pass"] for item in checks)
        else None,
    }

    if all(item["pass"] for item in checks):
        status = "P802_CURRENT_STRICT_ALPHA_S_PASSIVE_PROVIDER_SKELETON_SUPPORTED_ACTIVE_PROVIDER_ACTION_RULE_BLOCKED"
    else:
        status = "P802_REQUIRES_REVIEW"

    artifact = {
        "stage": "P802",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p801_provider_skeleton_probe": rel(IN_P801),
            "f801_provider_skeleton_packet": rel(IN_F801),
            "f800_provider_class_target": rel(IN_F800),
            "f704_packet": rel(IN_F704),
            "n705_os_support_theorem": rel(IN_N705),
            "n706_sign_gauge_theorem": rel(IN_N706),
            "n703_meaning_theorem": rel(IN_N703),
            "p696_computability_probe": rel(IN_P696),
            "p709_sign_gauge_audit": rel(IN_P709),
            "p710_nonstrict_proxy_to_gev_probe": rel(IN_P710),
            "proxy_to_gev_policy": rel(IN_POLICY),
        },
        "provider_skeleton_ref": skeleton.get("object_id"),
        "provider_class_target_ref": skeleton.get("provider_class_target_ref"),
        "clause_split": clause_split,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The current same-domain alpha_s lane already exports the passive provider skeleton beneath the missing provider class.",
            "That skeleton supplies passive roles only: observable carrier, packaging, gauge safety, computability, and meaning discipline.",
            "What remains missing is the active provider action rule that would make this passive base act as a provider class.",
        ],
        "recommended_next_packet": {
            "id": "F802_CURRENT_STRICT_ALPHA_S_PROVIDER_ACTION_RULE_TARGET_PACKET",
            "goal": "Freeze the exact active-rule object still missing before the passive same-domain skeleton can count as a real provider class.",
            "minimum_fields": [
                "provider_skeleton_ref",
                "support_bundle_ref",
                "acting_same_domain_input_ref",
                "provider_action_rule_ref",
                "semantic_principle_supply_ref",
                "passive_to_active_upgrade_block_ref",
                "foreign_domain_exclusion_ref",
                "selected_provider_class_output_schema",
                "hard_limits",
            ],
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P802",
        "status": status,
        "as_of": AS_OF,
        "passive_skeleton_clause_status": clause_split["passive_skeleton_clause_status"],
        "active_provider_action_rule_clause_status": clause_split["active_provider_action_rule_clause_status"],
        "sharp_blocker_clause": clause_split["sharp_blocker_clause"],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
