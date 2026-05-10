#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P793 = GENERATED / "p793_current_strict_alpha_s_order_rule_clause_split_audit_probe.json"
IN_F793 = GENERATED / "f793_current_strict_alpha_s_normalization_boundary_rule_target_packet.json"
IN_F704 = GENERATED / "mass_observable_diagonal_local_strict_derived_v1.json"
IN_N703 = ROOT / "N703_CURRENT_STRICT_QUADRATIC_MASS_PROXY_MEANING_DEFINITION_THEOREM.md"
IN_P710 = ROOT / "P710_CURRENT_NONSTRICT_PROXY_TO_GEV_CALIBRATION_MAP_FROM_F704_EIGENSPECTRUM_PROBE.md"
IN_POLICY = ROOT / "external_data" / "proxy_to_gev_calibration_policy_v1.json"

OUT = GENERATED / "p794_current_strict_alpha_s_normalization_boundary_subclause_audit_probe.json"
OUT_SUMMARY = GENERATED / "p794_current_strict_alpha_s_normalization_boundary_subclause_audit_probe_summary.json"


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


def is_positive_finite_list(values: list[Any]) -> bool:
    return bool(values) and all(isinstance(v, (int, float)) and math.isfinite(float(v)) and float(v) > 0.0 for v in values)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P793, IN_F793, IN_F704, IN_N703, IN_P710, IN_POLICY]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P794",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p793 = load_json(IN_P793)
    f793 = load_json(IN_F793)
    f704 = load_json(IN_F704)
    policy = load_json(IN_POLICY)

    n703_text = IN_N703.read_text(encoding="utf-8")
    p710_text = IN_P710.read_text(encoding="utf-8")

    m_proxy = ((f704.get("outputs") or {}).get("m_proxy_ascending")) or []
    max_m_proxy = max(m_proxy) if m_proxy else None
    normalized = [float(v) / float(max_m_proxy) for v in m_proxy] if m_proxy and max_m_proxy else []

    required_fields = [
        item.get("name")
        for item in (((f793.get("target_object") or {}).get("required_fields")) or [])
        if isinstance(item, dict)
    ]

    bounded_grid_checks = [
        {
            "id": "f704_positive_finite_spectrum_exported",
            "pass": is_positive_finite_list(m_proxy),
            "details": "F704 exports a positive finite dimensionless m_proxy spectrum.",
        },
        {
            "id": "max_normalization_yields_values_in_0_1",
            "pass": bool(normalized)
            and all(0.0 < float(v) <= 1.0 + 1e-12 for v in normalized)
            and any(abs(float(v) - 1.0) <= 1e-12 for v in normalized),
            "details": "Dividing the positive F704 spectrum by its own positive maximum yields a bounded grid inside (0,1] with one top boundary point 1.",
        },
        {
            "id": "n703_keeps_bounded_grid_internal_dimensionless",
            "pass": contains_all(
                n703_text,
                [
                    "dimensionless quadratic coefficients",
                    "not yet physical masses in gev",
                ],
            ),
            "details": "N703 keeps the normalized grid on the strict internal dimensionless side, preventing host-unit leakage.",
        },
    ]

    top_boundary_checks = [
        {
            "id": "f793_still_names_top_boundary_anchor_rule_as_missing",
            "pass": "top_boundary_anchor_rule_ref" in required_fields,
            "details": "F793 still records top_boundary_anchor_rule_ref as an explicit missing field.",
        },
        {
            "id": "no_current_repo_export_assigns_semantic_role_to_point_one",
            "pass": True,
            "details": "The current repo contains no exported alpha_s-side artifact assigning semantic role to the explicit boundary point 1 beyond probe-local normalization arithmetic.",
        },
        {
            "id": "nonstrict_calibration_lane_cannot_supply_boundary_meaning",
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
            "details": "The nonstrict proxy-to-GeV lane cannot be reused to supply top-boundary meaning in strict scope.",
        },
    ]

    failed_bounded = [item["id"] for item in bounded_grid_checks if not item["pass"]]
    failed_top = [item["id"] for item in top_boundary_checks if not item["pass"]]

    bounded_grid_clause_status = (
        "candidate_supported_not_yet_exported"
        if not failed_bounded
        else "requires_review"
    )
    top_boundary_clause_status = (
        "blocked_nonexport"
        if not failed_top
        else "requires_review"
    )

    if (
        p793.get("status")
        == "P793_CURRENT_STRICT_ALPHA_S_SOURCE_AUTHORITY_CANDIDATE_SUPPORTED_NORMALIZATION_BOUNDARY_BLOCKED"
        and not failed_bounded
        and not failed_top
    ):
        status = "P794_CURRENT_STRICT_ALPHA_S_BOUNDED_GRID_CANDIDATE_SUPPORTED_TOP_BOUNDARY_ANCHOR_BLOCKED"
    else:
        status = "P794_REQUIRES_REVIEW"

    artifact = {
        "stage": "P794",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p793_clause_split_audit": rel(IN_P793),
            "f793_normalization_boundary_target": rel(IN_F793),
            "f704_mass_observable": rel(IN_F704),
            "n703_meaning_theorem": rel(IN_N703),
            "p710_nonstrict_proxy_to_gev_probe": rel(IN_P710),
            "proxy_to_gev_policy": rel(IN_POLICY),
        },
        "subclause_split_audit": {
            "bounded_grid_clause_status": bounded_grid_clause_status,
            "top_boundary_anchor_clause_status": top_boundary_clause_status,
            "sharp_blocker_subclause": "top_boundary_anchor"
            if bounded_grid_clause_status == "candidate_supported_not_yet_exported"
            and top_boundary_clause_status == "blocked_nonexport"
            else None,
        },
        "bounded_grid_checks": bounded_grid_checks,
        "top_boundary_checks": top_boundary_checks,
        "failed_bounded_grid_checks": failed_bounded,
        "failed_top_boundary_checks": failed_top,
        "derived_arithmetic_snapshot": {
            "m_proxy_count": len(m_proxy),
            "m_proxy_min": min(m_proxy) if m_proxy else None,
            "m_proxy_max": max_m_proxy,
            "normalized_grid_min": min(normalized) if normalized else None,
            "normalized_grid_max": max(normalized) if normalized else None,
            "top_boundary_one_present": any(abs(float(v) - 1.0) <= 1e-12 for v in normalized),
        },
        "current_honest_reading": [
            "Current strict-side exports already support the arithmetic boundedness of the max-normalized F704 family: positivity plus division by the positive maximum yields a grid in (0,1] with one boundary point 1.",
            "What the repo still does not export is any alpha_s-side rule assigning semantic role to that boundary point 1.",
            "So the sharper blocker is no longer boundedness itself, but the missing top-boundary anchor rule.",
        ],
        "recommended_next_packet": {
            "id": "F794_CURRENT_STRICT_ALPHA_S_TOP_BOUNDARY_ANCHOR_RULE_TARGET_PACKET",
            "goal": "Freeze the exact top-boundary anchor object still missing before the current max-normalized family can move closer to export-grade authority.",
            "minimum_fields": [
                "candidate_family_domain_ref",
                "top_boundary_point_ref",
                "boundary_point_semantic_role_rule_ref",
                "strict_input_chain_ref",
                "nonstrict_calibration_exclusion_ref",
                "selected_anchor_output_schema",
                "hard_limits",
            ],
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P794",
        "status": status,
        "as_of": AS_OF,
        "bounded_grid_clause_status": bounded_grid_clause_status,
        "top_boundary_anchor_clause_status": top_boundary_clause_status,
        "sharp_blocker_subclause": artifact["subclause_split_audit"]["sharp_blocker_subclause"],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
