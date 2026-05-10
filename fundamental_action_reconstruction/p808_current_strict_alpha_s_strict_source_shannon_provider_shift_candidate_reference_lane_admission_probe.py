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
IN_F807 = GENERATED / "f807_current_strict_alpha_s_provider_class_shift_requirement_packet.json"
IN_P807 = GENERATED / "p807_current_strict_alpha_s_strict_source_shannon_reuse_as_same_domain_provider_action_source_audit_probe.json"
IN_ALPHA = GENERATED / "alpha_geo_strict_derived_v1.json"
IN_P748 = GENERATED / "p748_current_strict_t202_strict_source_shannon_pair_population_refinement_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe.json"
IN_P749 = GENERATED / "p749_current_strict_t203_strict_source_shannon_pair_population_support_refinement_to_pair_indexed_population_anchor_entry_nonexport_audit_probe_summary.json"
IN_P754 = GENERATED / "p754_current_strict_t208_strict_source_shannon_minimal_designated_pair12_entry_lane_provider_shift_requirement_boundary_audit_probe_summary.json"
IN_F321 = ROOT / "F321_FIRST_ACTUAL_STRICT_NAD12_SIGMA_STRICT_SOURCE_SHANNON_WEIGHTED_PAIR_POPULATION_REFINEMENT_CANDIDATE_PACKET.md"
IN_N750 = ROOT / "N750_CURRENT_STRICT_T208_STRICT_SOURCE_SHANNON_MINIMAL_DESIGNATED_PAIR12_ENTRY_LANE_PROVIDER_SHIFT_REQUIREMENT_BOUNDARY_THEOREM.md"

OUT = GENERATED / "p808_current_strict_alpha_s_strict_source_shannon_provider_shift_candidate_reference_lane_admission_probe.json"
OUT_SUMMARY = GENERATED / "p808_current_strict_alpha_s_strict_source_shannon_provider_shift_candidate_reference_lane_admission_probe_summary.json"


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
        IN_F807,
        IN_P807,
        IN_ALPHA,
        IN_P748,
        IN_P749,
        IN_P754,
        IN_F321,
        IN_N750,
    ]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P808",
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
    f807 = load_json(IN_F807)
    p807 = load_json(IN_P807)
    alpha = load_json(IN_ALPHA)
    p748 = load_json(IN_P748)
    p749 = load_json(IN_P749)
    p754 = load_json(IN_P754)
    f321_text = IN_F321.read_text(encoding="utf-8")
    n750_text = IN_N750.read_text(encoding="utf-8")

    provider_target = f800.get("target_object") or {}
    shift_requirement = f807.get("exported_object") or {}
    p807_theorem = p807.get("theorem_result") or {}

    checks = [
        {
            "id": "f807_already_requires_provider_class_shift",
            "pass": (
                f807.get("status")
                == "F807_EXECUTED_CURRENT_STRICT_ALPHA_S_PROVIDER_CLASS_SHIFT_REQUIREMENT_PACKET_NO_FALSE_PASS"
                and shift_requirement.get("object_id") == "alpha_s_provider_class_shift_requirement_v1"
                and shift_requirement.get("remaining_admitted_move_class")
                == "shift_to_a_different_provider_class_lane"
                and shift_requirement.get("candidate_shift_lane_hint")
                == "strict_source_shannon_provider_shift_lane_candidate_not_yet_alpha_s_domain_audited"
            ),
            "details": "F807 already freezes provider-class shift as the remaining admitted continuation class for alpha_s.",
        },
        {
            "id": "strict_source_shannon_infrastructure_is_real",
            "pass": (
                alpha.get("object") == "alpha_geo_strict_derived_v1"
                and alpha.get("status") == "actual_exported_strict_derived_source_upgrade_value"
                and p748.get("current_strict_source_shannon_source_upgrades_exported") is True
                and p749.get("current_strict_source_shannon_pair_population_support_refinement_candidate_exported")
                is True
                and contains_all(
                    f321_text,
                    [
                        "strict-source Shannon-weighted pair-population refinement candidate",
                        "strict_weight: alpha_geo_strict_derived_v1",
                    ],
                )
            ),
            "details": "The strict-source alpha_geo/Shannon lane is real and exported at source/candidate grade.",
        },
        {
            "id": "strict_source_shannon_lane_remains_candidate_grade_reference_only",
            "pass": (
                p748.get("current_strict_source_shannon_pair_population_refinement_route_remains_unbridged_to_pair12_typed_carrier")
                is True
                and p749.get("current_strict_source_shannon_pair_population_support_refinement_route_remains_nonentering_for_pair_indexed_population_anchor")
                is True
                and p754.get("same_level_shannon_entry_lane_continuation_no_longer_admitted_primary_move")
                is True
                and contains_all(
                    n750_text,
                    [
                        "no longer an admitted primary move",
                        "or shift to a different provider class",
                    ],
                )
            ),
            "details": "The strict-source Shannon lane remains candidate-grade and reference-only; it is not yet an entered or realized provider lane.",
        },
        {
            "id": "alpha_s_domain_admission_still_blocked",
            "pass": (
                p807_theorem.get("strict_source_shannon_lane_is_same_domain_provider_action_source_for_alpha_s")
                is False
                and any(
                    isinstance(item, dict) and item.get("name") == "same_domain_carrier_ref"
                    for item in (provider_target.get("required_fields") or [])
                )
                and any(
                    isinstance(item, dict) and item.get("name") == "foreign_domain_reuse_block_ref"
                    for item in (provider_target.get("required_fields") or [])
                )
            ),
            "details": "Alpha_s-domain admission is still blocked: same-domain carrier discipline remains mandatory and foreign-domain reuse is fenced off.",
        },
        {
            "id": "strict_source_shannon_lane_can_be_admitted_as_provider_shift_candidate_reference_lane_only",
            "pass": (
                p807_theorem.get("next_honest_move_requires_provider_class_shift") is True
                and p807_theorem.get(
                    "current_repo_exports_no_genuinely_new_same_domain_provider_action_source_for_alpha_s"
                )
                is True
            ),
            "details": "Given the failed same-domain reuse and the real existence of a separate strict-source Shannon candidate lane, the honest admission grade is reference-lane candidate only.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]

    theorem_result = {
        "strict_source_shannon_provider_shift_candidate_reference_lane_admitted": all(
            item["pass"] for item in checks
        ),
        "strict_source_shannon_provider_shift_candidate_reference_lane_grade": "reference_context_candidate_only"
        if all(item["pass"] for item in checks)
        else None,
        "alpha_s_domain_interface_exported": False if all(item["pass"] for item in checks) else None,
        "provider_shift_realized": False if all(item["pass"] for item in checks) else None,
        "next_honest_move_requires_alpha_s_domain_interface_audit_for_shift_candidate": all(
            item["pass"] for item in checks
        ),
        "no_false_pass": True,
    }

    if all(item["pass"] for item in checks):
        status = "P808_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_PROVIDER_SHIFT_CANDIDATE_REFERENCE_LANE_ADMITTED_ALPHA_S_DOMAIN_INTERFACE_BLOCKED"
    else:
        status = "P808_REQUIRES_REVIEW"

    artifact = {
        "stage": "P808",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f800_provider_class_target": rel(IN_F800),
            "f807_provider_class_shift_requirement_packet": rel(IN_F807),
            "p807_same_domain_source_negative_audit": rel(IN_P807),
            "alpha_geo_strict_source": rel(IN_ALPHA),
            "p748_strict_source_shannon_pair12_bridge_boundary": rel(IN_P748),
            "p749_strict_source_shannon_pair_indexed_entry_boundary": rel(IN_P749),
            "p754_strict_source_shannon_provider_shift_boundary": rel(IN_P754),
            "f321_strict_source_shannon_candidate_packet": rel(IN_F321),
            "n750_strict_source_shannon_provider_shift_requirement_theorem": rel(IN_N750),
        },
        "provider_shift_candidate_reference_lane": {
            "candidate_id": "strict_source_shannon_provider_shift_candidate_reference_lane_v1",
            "shift_requirement_ref": shift_requirement.get("object_id"),
            "source_upgrade_ref": alpha.get("object"),
            "supporting_candidate_lane_refs": [
                rel(IN_P748),
                rel(IN_P749),
                rel(IN_P754),
                rel(IN_F321),
                rel(IN_N750),
            ],
            "candidate_grade": theorem_result["strict_source_shannon_provider_shift_candidate_reference_lane_grade"],
            "alpha_s_domain_interface_status": "blocked_nonexport",
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The strict-source Shannon lane is real enough to be named as a provider-shift candidate reference lane for alpha_s.",
            "It is not admitted as a same-domain alpha_s source and it does not yet enter the alpha_s domain.",
            "So the strongest honest export is candidate-reference status only, with alpha_s-domain interface still blocked.",
        ],
        "recommended_next_packet": {
            "id": "F808_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_PROVIDER_SHIFT_CANDIDATE_REFERENCE_PACKET",
            "goal": "Export the admitted strict-source Shannon provider-shift candidate reference lane for alpha_s while keeping alpha_s-domain interface blocked.",
            "export_object_id": "alpha_s_strict_source_shannon_provider_shift_candidate_reference_lane_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P808",
        "status": status,
        "as_of": AS_OF,
        "strict_source_shannon_provider_shift_candidate_reference_lane_admitted": theorem_result[
            "strict_source_shannon_provider_shift_candidate_reference_lane_admitted"
        ],
        "alpha_s_domain_interface_exported": theorem_result["alpha_s_domain_interface_exported"],
        "provider_shift_realized": theorem_result["provider_shift_realized"],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
