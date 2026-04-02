#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P789 = GENERATED / "p789_current_strict_alpha_s_normalized_boundary_interface_candidate_probe.json"
IN_P791 = GENERATED / "p791_current_strict_alpha_s_selection_principle_reuse_audit_probe.json"
IN_F704 = ROOT / "F704_CURRENT_STRICT_INVARIANT_MASS_OBSERVABLE_FROM_DIAGONAL_LOCAL_PSI_HESSIAN_EIGENSYSTEM_EXPORT_PACKET.md"
IN_P694 = ROOT / "P694_CURRENT_STRICT_PHYSICAL_COMPUTABILITY_MASS_SPECTRUM_PROXY_FROM_PROJECTIVE_SELECTOR_CLOSURE_PROBE.md"
IN_POLICY = ROOT / "external_data" / "proxy_to_gev_calibration_policy_v1.json"

OUT = GENERATED / "p792_current_strict_alpha_s_family_selection_order_rule_probe.json"
OUT_SUMMARY = GENERATED / "p792_current_strict_alpha_s_family_selection_order_rule_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def approx_equal(a: float, b: float, tol: float = 1e-12) -> bool:
    return abs(float(a) - float(b)) <= tol


def is_unit_bounded_grid(values: list[float]) -> bool:
    if not values:
        return False
    return all(float(v) > 0.0 and float(v) <= 1.0 + 1e-12 for v in values) and approx_equal(max(values), 1.0)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P789, IN_P791, IN_F704, IN_P694, IN_POLICY]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P792",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p789 = load_json(IN_P789)
    p791 = load_json(IN_P791)
    policy = load_json(IN_POLICY)

    f704_text = IN_F704.read_text(encoding="utf-8")
    p694_text = IN_P694.read_text(encoding="utf-8")

    families = p789.get("candidate_families") or {}
    domain_ids = list(families.keys())

    scored: list[dict[str, Any]] = []
    for family_id, family in families.items():
        values = [float(v) for v in family.get("normalized_validation_points_candidate") or []]
        source_object = str(family.get("source_object"))
        source_tier = 0
        basis_invariant_full_spectrum = False

        if source_object == "MassObservable_diagonal_local_strict_derived_v1":
            source_tier = 2
            basis_invariant_full_spectrum = True
        elif "pair-rayleigh" in source_object.lower():
            source_tier = 1

        score_tuple = (
            source_tier,
            1 if is_unit_bounded_grid(values) else 0,
            len(values),
        )

        scored.append(
            {
                "family_id": family_id,
                "source_object": source_object,
                "source_tier": source_tier,
                "basis_invariant_full_spectrum": basis_invariant_full_spectrum,
                "unit_bounded_grid": is_unit_bounded_grid(values),
                "validation_point_count": len(values),
                "max_validation_point": max(values) if values else None,
                "min_validation_point": min(values) if values else None,
                "normalization_rule": family.get("normalization_rule"),
                "score_tuple_lexicographic": score_tuple,
            }
        )

    scored_sorted = sorted(scored, key=lambda item: item["score_tuple_lexicographic"], reverse=True)
    top = scored_sorted[0]
    second = scored_sorted[1] if len(scored_sorted) > 1 else None
    unique_winner = second is None or top["score_tuple_lexicographic"] != second["score_tuple_lexicographic"]

    geometric_mean_policy_is_nonstrict = (
        policy.get("scope") == "nonstrict_physical_units_calibration_only"
        and "geometric_mean" in str(policy.get("policy_id"))
    )

    checks = [
        {
            "id": "candidate_family_domain_present",
            "pass": p789.get("candidate_family_constructible") is True and len(domain_ids) == 3,
            "details": "The current alpha_s family domain is the finite three-family domain exported by P789.",
        },
        {
            "id": "source_authority_rule_prefers_f704_over_p694",
            "pass": (
                "basis-invariant" in f704_text
                and "projective" in p694_text
                and any(item["family_id"] == "p694_pairmax_anchor_family" and item["source_tier"] == 1 for item in scored)
                and any(item["family_id"].startswith("f704_") and item["source_tier"] == 2 for item in scored)
            ),
            "details": "F704 exports a basis-invariant full-spectrum strict object, while P694 is a projective pair-summary probe; this supports a probe-local source-authority preference.",
        },
        {
            "id": "normalization_boundary_rule_prefers_unit_bounded_grid",
            "pass": (
                any(item["family_id"] == "f704_max_mode_anchor_family" and item["unit_bounded_grid"] for item in scored)
                and any(item["family_id"] == "f704_geometric_mean_anchor_family" and not item["unit_bounded_grid"] for item in scored)
            ),
            "details": "On the current family domain, the max-mode family stays inside (0,1] with explicit top boundary 1, while the geometric-mean family exceeds 1.",
        },
        {
            "id": "geometric_mean_already_has_explicit_nonstrict_role_only",
            "pass": geometric_mean_policy_is_nonstrict,
            "details": "The only explicit exported geometric-mean calibration role currently present in the repo is the nonstrict proxy-to-GeV policy layer.",
        },
        {
            "id": "unique_winner_exists_but_selection_export_remains_missing",
            "pass": (
                unique_winner
                and top["family_id"] == "f704_max_mode_anchor_family"
                and p791.get("status") == "P791_CURRENT_STRICT_ALPHA_S_SELECTION_PRINCIPLE_REUSE_NEGATIVE_ON_CURRENT_REPO_STATE"
            ),
            "details": "A unique winner exists under the explicit probe-local order tuple, but P791 still blocks any silent upgrade into an exported strict selection principle.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]

    if not failed_checks and unique_winner:
        status = "P792_CURRENT_STRICT_ALPHA_S_PROBE_LOCAL_ORDER_RULE_UNIQUE_WINNER_NONEXPORT"
    else:
        status = "P792_REQUIRES_REVIEW"

    artifact = {
        "stage": "P792",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p789_candidate_family_probe": rel(IN_P789),
            "p791_selection_reuse_audit": rel(IN_P791),
            "f704_invariant_mass_observable_packet": rel(IN_F704),
            "p694_projective_pair_proxy_probe": rel(IN_P694),
            "proxy_to_gev_policy": rel(IN_POLICY),
        },
        "candidate_family_domain_ids": domain_ids,
        "probe_local_order_rule_tuple": {
            "rule_id": "probe_local_alpha_s_family_order_rule_v1",
            "clauses_in_priority_order": [
                {
                    "id": "source_authority_preference",
                    "meaning": "Prefer strict basis-invariant full-spectrum source objects over weaker projective pair-summary proxies.",
                },
                {
                    "id": "unit_bounded_boundary_preference",
                    "meaning": "Among surviving families, prefer normalized validation grids confined to (0,1] with explicit top boundary point 1.",
                },
                {
                    "id": "coverage_count_tie_break",
                    "meaning": "Use validation-point coverage count only as a residual tie-break on the current finite domain.",
                },
            ],
            "export_status": "probe_local_only",
        },
        "family_evaluation": scored_sorted,
        "unique_winner_exists": unique_winner,
        "unique_winner_family_id": top["family_id"] if unique_winner else None,
        "supporting_boundary_notes": [
            "F704 is explicitly exported as the basis-invariant full-spectrum mass observable object.",
            "P694 remains a projective pair-summary probe-level proxy object.",
            "The exported geometric-mean calibration policy is explicitly nonstrict and physical-units-facing; this blocks silent reuse of geometric-mean language as a strict selection authority.",
        ],
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The current alpha_s family-space is no longer selection-chaotic: one explicit probe-local order tuple already isolates a unique winner.",
            "That winner is f704_max_mode_anchor_family, not because the repo already exported a strict selection principle, but because the current finite family domain plus explicit ranking clauses leave no tie.",
            "The blocker is now narrower still: what remains missing is not a winner on the current domain, but export-grade authority for the order rule itself.",
        ],
        "recommended_next_packet": {
            "id": "F792_CURRENT_STRICT_ALPHA_S_FAMILY_SELECTION_ORDER_RULE_TARGET_PACKET",
            "goal": "Freeze the exact export-grade order-rule object still missing before the probe-local winner can count as a real alpha_s family selection result.",
            "minimum_fields": [
                "candidate_family_domain_ref",
                "source_authority_rule_ref",
                "normalization_boundary_rule_ref",
                "residual_tie_break_rule_ref",
                "nonstrict_calibration_separation_ref",
                "selected_family_output_schema",
                "hard_limits",
            ],
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P792",
        "status": status,
        "as_of": AS_OF,
        "unique_winner_exists": unique_winner,
        "unique_winner_family_id": artifact["unique_winner_family_id"],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
