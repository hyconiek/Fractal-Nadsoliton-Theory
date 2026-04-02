#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F800 = GENERATED / "f800_current_strict_alpha_s_reference_scale_provider_class_target_packet.json"
IN_F805 = GENERATED / "f805_current_strict_alpha_s_acting_input_bundle_packet.json"
IN_F806 = GENERATED / "f806_current_strict_alpha_s_provider_action_continuation_boundary_packet.json"
IN_ALPHA = GENERATED / "alpha_geo_strict_derived_v1.json"
IN_P748 = GENERATED / "p748_current_strict_t202_strict_source_shannon_pair_population_refinement_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe.json"
IN_P749 = GENERATED / "p749_current_strict_t203_strict_source_shannon_pair_population_support_refinement_to_pair_indexed_population_anchor_entry_nonexport_audit_probe_summary.json"
IN_P754 = GENERATED / "p754_current_strict_t208_strict_source_shannon_minimal_designated_pair12_entry_lane_provider_shift_requirement_boundary_audit_probe_summary.json"
IN_F321 = ROOT / "F321_FIRST_ACTUAL_STRICT_NAD12_SIGMA_STRICT_SOURCE_SHANNON_WEIGHTED_PAIR_POPULATION_REFINEMENT_CANDIDATE_PACKET.md"
IN_P422 = ROOT / "P422_CURRENT_STRICT_SHANNON_SYMMETRY_BREAKING_SELECTOR_INGREDIENT_QW2191_ATTACK_AUDIT_PROBE.md"

OUT = GENERATED / "p807_current_strict_alpha_s_strict_source_shannon_reuse_as_same_domain_provider_action_source_audit_probe.json"
OUT_SUMMARY = GENERATED / "p807_current_strict_alpha_s_strict_source_shannon_reuse_as_same_domain_provider_action_source_audit_probe_summary.json"


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
        IN_F800,
        IN_F805,
        IN_F806,
        IN_ALPHA,
        IN_P748,
        IN_P749,
        IN_P754,
        IN_F321,
        IN_P422,
    ]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P807",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    f800 = load_json(IN_F800)
    f805 = load_json(IN_F805)
    f806 = load_json(IN_F806)
    alpha = load_json(IN_ALPHA)
    p748 = load_json(IN_P748)
    p749 = load_json(IN_P749)
    p754 = load_json(IN_P754)
    f321_text = IN_F321.read_text(encoding="utf-8")
    p422_text = IN_P422.read_text(encoding="utf-8")

    provider_target = f800.get("target_object") or {}
    acting_input = f805.get("exported_object") or {}
    continuation_boundary = f806.get("exported_object") or {}

    checks = [
        {
            "id": "f806_already_requires_new_same_domain_source_or_provider_shift",
            "pass": (
                f806.get("status")
                == "F806_EXECUTED_CURRENT_STRICT_ALPHA_S_PROVIDER_ACTION_CONTINUATION_BOUNDARY_PACKET_NO_FALSE_PASS"
                and continuation_boundary.get("object_id") == "alpha_s_provider_action_continuation_boundary_v1"
                and continuation_boundary.get("admitted_next_move_classes")
                == [
                    "export_one_genuinely_new_same_domain_provider_action_source",
                    "shift_to_a_different_provider_class_lane",
                ]
            ),
            "details": "F806 already narrows the continuation to one genuinely new same-domain source or a provider-class shift.",
        },
        {
            "id": "alpha_geo_strict_source_object_is_real",
            "pass": (
                alpha.get("object") == "alpha_geo_strict_derived_v1"
                and alpha.get("status") == "actual_exported_strict_derived_source_upgrade_value"
                and alpha.get("value") == "4 ln(2)"
            ),
            "details": "alpha_geo_strict_derived_v1 is a real strict-derived source object, not a legacy import.",
        },
        {
            "id": "strict_source_shannon_candidate_lane_exists",
            "pass": (
                p749.get("current_strict_source_shannon_source_upgrades_exported") is True
                and p749.get("current_strict_source_shannon_pair_population_support_refinement_candidate_exported")
                is True
                and contains_all(
                    f321_text,
                    [
                        "strict-source Shannon-weighted pair-population refinement candidate",
                        "still below actual pair population",
                        "still below actual theta export",
                    ],
                )
            ),
            "details": "The repo does export strict-source alpha_geo/Shannon candidate lanes.",
        },
        {
            "id": "strict_source_shannon_lane_remains_candidate_only_nonentering_and_unbridged",
            "pass": (
                p749.get("current_strict_source_shannon_pair_population_support_refinement_candidate_remains_below_actual_pair_population_and_theta_export")
                is True
                and p749.get("current_strict_source_shannon_pair_population_support_refinement_route_remains_nonentering_for_pair_indexed_population_anchor")
                is True
                and p748.get("current_strict_source_shannon_pair_population_refinement_route_remains_unbridged_to_pair12_typed_carrier")
                is True
                and p748.get("current_strict_source_shannon_pair_population_refinement_route_has_exported_pair12_typed_carrier_bridge")
                is False
                and p754.get("same_level_shannon_entry_lane_continuation_no_longer_admitted_primary_move")
                is True
            ),
            "details": "The strict-source Shannon route remains candidate-only, nonentering, and unbridged even on its own declared future carrier lanes.",
        },
        {
            "id": "current_alpha_s_lane_still_requires_same_domain_carrier_discipline",
            "pass": (
                acting_input.get("object_id") == "alpha_s_reference_scale_acting_input_bundle_v1"
                and acting_input.get("scope") == "strict_same_domain_acting_input_only"
                and provider_target.get("object_id") == "alpha_s_reference_scale_provider_class_target_v1"
                and any(
                    isinstance(item, dict) and item.get("name") == "same_domain_carrier_ref"
                    for item in (provider_target.get("required_fields") or [])
                )
                and any(
                    isinstance(item, dict) and item.get("name") == "foreign_domain_reuse_block_ref"
                    for item in (provider_target.get("required_fields") or [])
                )
            ),
            "details": "The current alpha_s lane still explicitly requires same-domain carrier discipline and fences off foreign-domain reuse.",
        },
        {
            "id": "no_current_strict_source_shannon_same_domain_provider_action_source_is_exported_for_alpha_s",
            "pass": (
                contains_all(
                    p422_text,
                    [
                        "No strict-core Shannon symmetry-breaking selector ingredient is exported today",
                        "not admissible as a strict-core claim",
                    ],
                )
                and p748.get("current_strict_source_shannon_pair_population_refinement_route_has_exported_pair12_typed_carrier_bridge")
                is False
                and p749.get("current_strict_source_shannon_pair_population_support_refinement_route_has_exported_pair_indexed_population_anchor_entry")
                is False
            ),
            "details": "Nothing currently exported upgrades strict-source alpha_geo/Shannon rhetoric into a same-domain provider-action source for the alpha_s lane.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]

    theorem_result = {
        "strict_source_alpha_geo_exists": checks[1]["pass"],
        "strict_source_shannon_candidate_lane_exists": checks[2]["pass"],
        "strict_source_shannon_lane_is_same_domain_provider_action_source_for_alpha_s": False
        if checks[3]["pass"] and checks[4]["pass"] and checks[5]["pass"]
        else None,
        "current_repo_exports_no_genuinely_new_same_domain_provider_action_source_for_alpha_s": all(
            item["pass"] for item in checks
        ),
        "next_honest_move_requires_provider_class_shift": all(item["pass"] for item in checks),
        "no_false_pass": True,
    }

    if all(item["pass"] for item in checks):
        status = "P807_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_REUSE_AS_SAME_DOMAIN_PROVIDER_ACTION_SOURCE_NEGATIVE_PROVIDER_SHIFT_REQUIRED"
    else:
        status = "P807_REQUIRES_REVIEW"

    artifact = {
        "stage": "P807",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f800_provider_class_target": rel(IN_F800),
            "f805_acting_input_bundle_packet": rel(IN_F805),
            "f806_provider_action_continuation_boundary_packet": rel(IN_F806),
            "alpha_geo_strict_source": rel(IN_ALPHA),
            "p748_strict_source_shannon_pair12_bridge_boundary": rel(IN_P748),
            "p749_strict_source_shannon_pair_indexed_entry_boundary": rel(IN_P749),
            "p754_strict_source_shannon_provider_shift_boundary": rel(IN_P754),
            "f321_strict_source_shannon_pair_population_candidate_packet": rel(IN_F321),
            "p422_strict_shannon_selector_attack_audit": rel(IN_P422),
        },
        "rejected_same_domain_source_candidate_class": {
            "candidate_id": "strict_source_alpha_geo_shannon_reuse_same_domain_provider_action_source_candidate_v1",
            "alpha_source_ref": alpha.get("object"),
            "current_alpha_s_acting_input_bundle_ref": acting_input.get("object_id"),
            "continuation_boundary_ref": continuation_boundary.get("object_id"),
            "rejection_basis_refs": [
                rel(IN_F800),
                rel(IN_F805),
                rel(IN_F806),
                rel(IN_ALPHA),
                rel(IN_P748),
                rel(IN_P749),
                rel(IN_P754),
                rel(IN_F321),
                rel(IN_P422),
            ],
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "alpha_geo_strict_derived_v1 and strict-source Shannon candidate lanes do exist on the current repo state.",
            "But those Shannon lanes remain candidate-only, nonentering, and unbridged even on their own declared routes.",
            "So they are not currently exported as a genuinely new same-domain provider-action source for the alpha_s lane.",
        ],
        "recommended_next_packet": {
            "id": "F807_CURRENT_STRICT_ALPHA_S_PROVIDER_CLASS_SHIFT_REQUIREMENT_PACKET",
            "goal": "Freeze the current-repo-state requirement that continuation now proceed by provider-class shift rather than same-domain verbal promotion of strict-source alpha_geo/Shannon language.",
            "export_object_id": "alpha_s_provider_class_shift_requirement_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P807",
        "status": status,
        "as_of": AS_OF,
        "strict_source_alpha_geo_exists": theorem_result["strict_source_alpha_geo_exists"],
        "strict_source_shannon_candidate_lane_exists": theorem_result["strict_source_shannon_candidate_lane_exists"],
        "current_repo_exports_no_genuinely_new_same_domain_provider_action_source_for_alpha_s": theorem_result[
            "current_repo_exports_no_genuinely_new_same_domain_provider_action_source_for_alpha_s"
        ],
        "next_honest_move_requires_provider_class_shift": theorem_result[
            "next_honest_move_requires_provider_class_shift"
        ],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
