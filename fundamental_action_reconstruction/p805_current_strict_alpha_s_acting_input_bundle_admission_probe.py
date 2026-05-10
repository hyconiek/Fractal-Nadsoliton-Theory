#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"
EPS = 1e-12

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P804 = GENERATED / "p804_current_strict_alpha_s_same_domain_compositional_alignment_bundle_admission_probe.json"
IN_F804 = GENERATED / "f804_current_strict_alpha_s_same_domain_compositional_alignment_bundle_packet.json"
IN_F802 = GENERATED / "f802_current_strict_alpha_s_provider_action_rule_target_packet.json"
IN_P789 = GENERATED / "p789_current_strict_alpha_s_normalized_boundary_interface_candidate_probe.json"
IN_P792 = GENERATED / "p792_current_strict_alpha_s_family_selection_order_rule_probe.json"
IN_P797 = GENERATED / "p797_current_strict_alpha_s_invariant_extremum_vs_reference_scale_audit_probe.json"
IN_P794 = GENERATED / "p794_current_strict_alpha_s_normalization_boundary_subclause_audit_probe.json"
IN_P710 = ROOT / "P710_CURRENT_NONSTRICT_PROXY_TO_GEV_CALIBRATION_MAP_FROM_F704_EIGENSPECTRUM_PROBE.md"
IN_POLICY = ROOT / "external_data" / "proxy_to_gev_calibration_policy_v1.json"

OUT = GENERATED / "p805_current_strict_alpha_s_acting_input_bundle_admission_probe.json"
OUT_SUMMARY = GENERATED / "p805_current_strict_alpha_s_acting_input_bundle_admission_probe_summary.json"


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
        IN_P804,
        IN_F804,
        IN_F802,
        IN_P789,
        IN_P792,
        IN_P797,
        IN_P794,
        IN_P710,
        IN_POLICY,
    ]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P805",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p804 = load_json(IN_P804)
    f804 = load_json(IN_F804)
    f802 = load_json(IN_F802)
    p789 = load_json(IN_P789)
    p792 = load_json(IN_P792)
    p797 = load_json(IN_P797)
    p794 = load_json(IN_P794)
    policy = load_json(IN_POLICY)
    p710_text = IN_P710.read_text(encoding="utf-8")

    alignment = f804.get("exported_object") or {}
    action_target = f802.get("target_object") or {}
    strongest_family_id = p789.get("strongest_family_id")
    strongest_family = ((p789.get("candidate_families") or {}).get(strongest_family_id)) or {}

    winner_family_id = p792.get("unique_winner_family_id")
    acting_mu0_tilde = strongest_family.get("mu0_tilde_candidate")
    acting_normalization_rule = strongest_family.get("normalization_rule")
    acting_points = strongest_family.get("normalized_validation_points_candidate") or []
    shared_max = float(p797.get("extremum_snapshot", {}).get("max_m_proxy"))
    arithmetic_max = float(p794.get("derived_arithmetic_snapshot", {}).get("m_proxy_max"))
    top_boundary_one_present = bool(p794.get("derived_arithmetic_snapshot", {}).get("top_boundary_one_present"))

    checks = [
        {
            "id": "f804_already_exports_same_domain_alignment_below_missing_action_rule",
            "pass": (
                p804.get("status")
                == "P804_CURRENT_STRICT_ALPHA_S_COMPOSITIONAL_ALIGNMENT_BUNDLE_EXPORT_ADMITTED_PROVIDER_ACTION_RULE_STILL_BLOCKED"
                and f804.get("status")
                == "F804_EXECUTED_CURRENT_STRICT_ALPHA_S_SAME_DOMAIN_COMPOSITIONAL_ALIGNMENT_BUNDLE_PACKET_NO_FALSE_PASS"
                and alignment.get("object_id") == "alpha_s_reference_scale_same_domain_compositional_alignment_bundle_v1"
                and alignment.get("provider_action_target_ref") == "alpha_s_reference_scale_provider_action_rule_target_v1"
            ),
            "details": "F804 already exports the aligned same-domain chain beneath the missing provider action rule.",
        },
        {
            "id": "winner_normalization_and_boundary_determine_one_acting_input_tuple",
            "pass": (
                winner_family_id == "f704_max_mode_anchor_family"
                and strongest_family_id == "f704_max_mode_anchor_family"
                and acting_normalization_rule == "mu_tilde_i := m_proxy_i / max_j(m_proxy_j)"
                and acting_mu0_tilde is not None
                and abs(float(acting_mu0_tilde) - 1.0) <= EPS
                and top_boundary_one_present
                and acting_points
                and abs(float(max(acting_points)) - 1.0) <= EPS
            ),
            "details": "The current same-domain chain determines one explicit acting tuple: winning family, max-normalization rule, and mu0_tilde/top-boundary point 1.",
        },
        {
            "id": "acting_input_tuple_stays_fully_dimensionless_and_same_domain",
            "pass": (
                strongest_family.get("checks", {}).get("dimensionless_points_present") is True
                and strongest_family.get("checks", {}).get("derived_from_exported_strict_objects") is True
                and strongest_family.get("checks", {}).get("no_gev_fields_used") is True
                and abs(shared_max - arithmetic_max) <= EPS
            ),
            "details": "The acting input tuple remains fully dimensionless, uses only exported strict objects, and stays on the same F704-based domain.",
        },
        {
            "id": "acting_input_bundle_still_stops_before_active_provider_action_rule",
            "pass": any(
                isinstance(item, dict) and item.get("name") == "acting_same_domain_input_ref"
                for item in (action_target.get("required_fields") or [])
            ) and any(
                isinstance(item, dict) and item.get("name") == "provider_action_rule_ref"
                for item in (action_target.get("required_fields") or [])
            ),
            "details": "F802 still records both the acting input ref and the active provider action rule as distinct missing layers.",
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
            "details": "The proxy-to-GeV calibration lane remains explicitly excluded from the acting-input question.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]

    clause_split = {
        "acting_input_bundle_clause_status": "export_admitted"
        if checks[1]["pass"] and checks[2]["pass"] and checks[4]["pass"]
        else "requires_review",
        "provider_action_rule_clause_status": "blocked_nonexport"
        if checks[0]["pass"] and checks[3]["pass"]
        else "requires_review",
        "sharp_blocker_clause": "provider_action_rule_ref"
        if all(item["pass"] for item in checks)
        else None,
    }

    if all(item["pass"] for item in checks):
        status = "P805_CURRENT_STRICT_ALPHA_S_ACTING_INPUT_BUNDLE_EXPORT_ADMITTED_PROVIDER_ACTION_RULE_STILL_BLOCKED"
    else:
        status = "P805_REQUIRES_REVIEW"

    artifact = {
        "stage": "P805",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p804_alignment_probe": rel(IN_P804),
            "f804_alignment_packet": rel(IN_F804),
            "f802_provider_action_target": rel(IN_F802),
            "p789_candidate_family_probe": rel(IN_P789),
            "p792_family_winner_relation": rel(IN_P792),
            "p797_extremum_relation": rel(IN_P797),
            "p794_bounded_boundary_relation": rel(IN_P794),
            "p710_nonstrict_proxy_to_gev_probe": rel(IN_P710),
            "proxy_to_gev_policy": rel(IN_POLICY),
        },
        "same_domain_acting_input_bundle_candidate": {
            "alignment_bundle_ref": alignment.get("object_id"),
            "provider_action_target_ref": alignment.get("provider_action_target_ref"),
            "acting_family_id": winner_family_id,
            "acting_normalization_rule": acting_normalization_rule,
            "acting_mu0_tilde": acting_mu0_tilde,
            "acting_top_boundary_point": 1.0,
            "acting_support_refs": [
                rel(IN_P789),
                rel(IN_P792),
                rel(IN_P797),
                rel(IN_P794),
            ],
        },
        "clause_split": clause_split,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The current same-domain alpha_s lane already determines one explicit acting input bundle beneath the missing provider action rule.",
            "That bundle selects the winning family, the winning max-normalization rule, and the acting point mu0_tilde = 1 on the same lane.",
            "What remains missing is the active rule that would act on this already-determined input bundle.",
        ],
        "recommended_next_packet": {
            "id": "F805_CURRENT_STRICT_ALPHA_S_ACTING_INPUT_BUNDLE_PACKET",
            "goal": "Export the admitted same-domain acting input bundle explicitly while keeping the provider action rule unresolved.",
            "export_object_id": "alpha_s_reference_scale_acting_input_bundle_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P805",
        "status": status,
        "as_of": AS_OF,
        "acting_input_bundle_clause_status": clause_split["acting_input_bundle_clause_status"],
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
