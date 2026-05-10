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

IN_P803 = GENERATED / "p803_current_strict_alpha_s_same_domain_relation_bundle_admission_probe.json"
IN_F803 = GENERATED / "f803_current_strict_alpha_s_same_domain_relation_bundle_packet.json"
IN_F802 = GENERATED / "f802_current_strict_alpha_s_provider_action_rule_target_packet.json"
IN_P792 = GENERATED / "p792_current_strict_alpha_s_family_selection_order_rule_probe.json"
IN_P797 = GENERATED / "p797_current_strict_alpha_s_invariant_extremum_vs_reference_scale_audit_probe.json"
IN_P794 = GENERATED / "p794_current_strict_alpha_s_normalization_boundary_subclause_audit_probe.json"
IN_P789 = GENERATED / "p789_current_strict_alpha_s_normalized_boundary_interface_candidate_probe.json"
IN_P710 = ROOT / "P710_CURRENT_NONSTRICT_PROXY_TO_GEV_CALIBRATION_MAP_FROM_F704_EIGENSPECTRUM_PROBE.md"
IN_POLICY = ROOT / "external_data" / "proxy_to_gev_calibration_policy_v1.json"

OUT = GENERATED / "p804_current_strict_alpha_s_same_domain_compositional_alignment_bundle_admission_probe.json"
OUT_SUMMARY = GENERATED / "p804_current_strict_alpha_s_same_domain_compositional_alignment_bundle_admission_probe_summary.json"


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
        IN_P803,
        IN_F803,
        IN_F802,
        IN_P792,
        IN_P797,
        IN_P794,
        IN_P789,
        IN_P710,
        IN_POLICY,
    ]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P804",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p803 = load_json(IN_P803)
    f803 = load_json(IN_F803)
    f802 = load_json(IN_F802)
    p792 = load_json(IN_P792)
    p797 = load_json(IN_P797)
    p794 = load_json(IN_P794)
    p789 = load_json(IN_P789)
    policy = load_json(IN_POLICY)
    p710_text = IN_P710.read_text(encoding="utf-8")

    relation_bundle = f803.get("exported_object") or {}
    action_target = f802.get("target_object") or {}
    strongest_family_id = p789.get("strongest_family_id")
    strongest_family = ((p789.get("candidate_families") or {}).get(strongest_family_id)) or {}
    normalized_points = strongest_family.get("normalized_validation_points_candidate") or []

    winner_family_id = p792.get("unique_winner_family_id")
    extremum_max = float(p797.get("extremum_snapshot", {}).get("max_m_proxy"))
    arithmetic_max = float(p794.get("derived_arithmetic_snapshot", {}).get("m_proxy_max"))
    top_boundary_one_present = bool(p794.get("derived_arithmetic_snapshot", {}).get("top_boundary_one_present"))
    normalized_max = max(normalized_points) if normalized_points else None

    checks = [
        {
            "id": "f803_already_exports_relation_bundle_below_missing_action_rule",
            "pass": (
                p803.get("status")
                == "P803_CURRENT_STRICT_ALPHA_S_SAME_DOMAIN_RELATION_BUNDLE_EXPORT_ADMITTED_PROVIDER_ACTION_RULE_STILL_BLOCKED"
                and f803.get("status")
                == "F803_EXECUTED_CURRENT_STRICT_ALPHA_S_SAME_DOMAIN_RELATION_BUNDLE_PACKET_NO_FALSE_PASS"
                and relation_bundle.get("object_id") == "alpha_s_reference_scale_same_domain_relation_bundle_v1"
                and relation_bundle.get("provider_action_target_ref") == "alpha_s_reference_scale_provider_action_rule_target_v1"
            ),
            "details": "F803 already exports the same-domain relation bundle beneath the missing provider action rule.",
        },
        {
            "id": "winner_relation_and_candidate_family_align_on_same_object",
            "pass": (
                winner_family_id == "f704_max_mode_anchor_family"
                and strongest_family_id == "f704_max_mode_anchor_family"
                and strongest_family.get("normalization_rule") == "mu_tilde_i := m_proxy_i / max_j(m_proxy_j)"
            ),
            "details": "P792 and P789 align on the same strongest family object: f704_max_mode_anchor_family.",
        },
        {
            "id": "extremum_and_boundary_relations_align_on_same_maximum",
            "pass": (
                abs(extremum_max - arithmetic_max) <= EPS
                and top_boundary_one_present
                and normalized_max is not None
                and abs(float(normalized_max) - 1.0) <= EPS
                and abs(float(strongest_family.get("mu0_tilde_candidate")) - 1.0) <= EPS
            ),
            "details": "P797, P794, and P789 align on one shared maximum: the same F704 maximum generates the normalized top boundary point 1.",
        },
        {
            "id": "relation_bundle_snapshot_is_compositionally_consistent",
            "pass": (
                relation_bundle.get("relation_snapshot", {}).get("family_winner_id") == winner_family_id
                and abs(float(relation_bundle.get("relation_snapshot", {}).get("extremum_max_m_proxy")) - extremum_max) <= EPS
                and relation_bundle.get("relation_snapshot", {}).get("bounded_grid_top_boundary_one_present") is True
            ),
            "details": "The already-exported relation bundle snapshot is compositionally consistent with the underlying probes.",
        },
        {
            "id": "alignment_bundle_still_stops_before_active_provider_action_rule",
            "pass": any(
                isinstance(item, dict) and item.get("name") == "provider_action_rule_ref"
                for item in (action_target.get("required_fields") or [])
            ),
            "details": "F802 still records provider_action_rule_ref as the missing active rule above this aligned chain.",
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
            "details": "The proxy-to-GeV calibration lane remains explicitly excluded from the compositional-alignment question.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]

    clause_split = {
        "alignment_bundle_clause_status": "export_admitted"
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
        status = "P804_CURRENT_STRICT_ALPHA_S_COMPOSITIONAL_ALIGNMENT_BUNDLE_EXPORT_ADMITTED_PROVIDER_ACTION_RULE_STILL_BLOCKED"
    else:
        status = "P804_REQUIRES_REVIEW"

    artifact = {
        "stage": "P804",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p803_relation_bundle_probe": rel(IN_P803),
            "f803_relation_bundle_packet": rel(IN_F803),
            "f802_provider_action_target": rel(IN_F802),
            "p792_family_winner_relation": rel(IN_P792),
            "p797_extremum_relation": rel(IN_P797),
            "p794_bounded_boundary_relation": rel(IN_P794),
            "p789_candidate_family_probe": rel(IN_P789),
            "p710_nonstrict_proxy_to_gev_probe": rel(IN_P710),
            "proxy_to_gev_policy": rel(IN_POLICY),
        },
        "same_domain_compositional_alignment_bundle_candidate": {
            "relation_bundle_ref": relation_bundle.get("object_id"),
            "provider_action_target_ref": relation_bundle.get("provider_action_target_ref"),
            "shared_family_id": winner_family_id,
            "shared_max_m_proxy": extremum_max,
            "shared_top_boundary_point": 1.0,
            "alignment_refs": [
                rel(IN_P792),
                rel(IN_P797),
                rel(IN_P794),
                rel(IN_P789),
            ],
        },
        "clause_split": clause_split,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The current same-domain alpha_s lane already exports an aligned compositional chain beneath the missing provider action rule.",
            "That chain aligns one shared family winner, one shared maximum, and one shared top-boundary point 1 on the same lane.",
            "What remains missing is the active rule that would make this aligned chain act as a provider.",
        ],
        "recommended_next_packet": {
            "id": "F804_CURRENT_STRICT_ALPHA_S_SAME_DOMAIN_COMPOSITIONAL_ALIGNMENT_BUNDLE_PACKET",
            "goal": "Export the admitted same-domain compositional alignment bundle explicitly while keeping the provider action rule unresolved.",
            "export_object_id": "alpha_s_reference_scale_same_domain_compositional_alignment_bundle_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P804",
        "status": status,
        "as_of": AS_OF,
        "alignment_bundle_clause_status": clause_split["alignment_bundle_clause_status"],
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
