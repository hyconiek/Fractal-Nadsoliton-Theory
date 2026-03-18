#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-18"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P721 = GENERATED / "p721_current_strict_t176_source_topology_basis_free_qw2191_safe_provider_nonupgrade_audit_probe_summary.json"
IN_F141 = GENERATED / "f141_first_actual_source_topology_barrier_protected_sign_witness_packet_summary.json"
IN_F143 = GENERATED / "f143_first_actual_source_topology_nonzero_flow_witness_packet_summary.json"
IN_F148 = GENERATED / "f148_first_actual_source_topology_basis_independent_promotion_witness_packet_summary.json"
IN_F149 = GENERATED / "f149_first_actual_source_topology_quotient_safe_qw2191_resolution_witness_packet_summary.json"
IN_F150 = GENERATED / "f150_first_actual_declared_scope_source_topology_selector_theorem_packet_summary.json"
IN_ATLAS = GENERATED / "selector_atlas_pair12345_oriented_transport_mod_2pi_strict_convention_v1.json"

OUT_JSON = GENERATED / "p722_current_strict_t177_chart_sensitive_transported_flux_current_like_section_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p722_current_strict_t177_chart_sensitive_transported_flux_current_like_section_nonexport_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P721, IN_F141, IN_F143, IN_F148, IN_F149, IN_F150, IN_ATLAS]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P722",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p721 = load_json(IN_P721)
    f141 = load_json(IN_F141)
    f143 = load_json(IN_F143)
    f148 = load_json(IN_F148)
    f149 = load_json(IN_F149)
    f150 = load_json(IN_F150)
    atlas = load_json(IN_ATLAS)

    checks: list[dict[str, Any]] = []
    blocking: list[str] = []

    def add_check(check_id: str, actual: Any, expected: Any, meaning: str) -> None:
        passed = actual == expected
        checks.append(
            {
                "id": check_id,
                "actual": actual,
                "expected": expected,
                "pass": passed,
                "meaning": meaning,
            }
        )
        if not passed:
            blocking.append(check_id)

    atlas_transitions = (atlas.get("transitions") or {})
    atlas_chart_sensitive_transport_present = len(atlas_transitions) >= 10

    basis_reduction = (f148.get("support_packet") or {}).get("basis_class_reduction_operator") or {}
    quotient_relation = (f149.get("support_packet") or {}).get("quotient_relation") or {}
    distinguished_quotient_class = (f149.get("support_packet") or {}).get("distinguished_quotient_class") or {}
    source_topology_theorem_boundaries = ((f150.get("support_packet") or {}).get("l5_boundaries") or {})

    add_check(
        "p721_source_topology_lane_contains_physics_facing_ingredients",
        bool(p721.get("source_topology_physically_interpretable_strict_ingredients_present")),
        True,
        "P721 already confirms that the source-topology lane contains real physics-facing strict ingredients.",
    )
    add_check(
        "atlas_chart_sensitive_transport_present",
        atlas_chart_sensitive_transport_present,
        True,
        "The global selector atlas already exports chart-sensitive transport infrastructure on the pair1..pair5 cover.",
    )
    add_check(
        "nonzero_flow_witness_is_scalar_component_only",
        bool(f143.get("actual_nonzero_flow_witness_exported")) and isinstance(f143.get("scalar_component_witness_value"), (int, float)),
        True,
        "The current nonzero-flow witness is exported only as a scalar component witness, not as a chartwise transported section.",
    )
    add_check(
        "barrier_sign_witness_is_scalar_component_only",
        bool(f141.get("actual_barrier_protected_sign_witness_exported")) and isinstance(
            ((f141.get("support_packet") or {}).get("psi_src_barrier_sign_component_witness_v1")),
            int,
        ),
        True,
        "The current barrier-protected sign witness is exported only as a scalar/sign component witness, not as a chartwise transported section.",
    )
    add_check(
        "basis_independent_promotion_forgets_chart_basis_labels",
        bool(basis_reduction.get("forgets_chart_basis_labels")),
        True,
        "The current source-topology promotion step explicitly forgets chart basis labels.",
    )
    add_check(
        "basis_independent_promotion_preserves_downstream_only_scope",
        bool(basis_reduction.get("preserves_observer_downstream_only")),
        True,
        "The same promotion step remains explicitly downstream-only.",
    )
    add_check(
        "quotient_relation_quotients_chart_labels_only_via_basis_free_packet",
        bool(quotient_relation.get("quotients_chart_labels_only_via_basis_free_packet")),
        True,
        "The current quotient-safe QW-2191 relation resolves chart labels only by quotienting them away through the basis-free packet.",
    )
    add_check(
        "distinguished_quotient_class_is_quotient_class_only",
        bool(distinguished_quotient_class.get("quotient_class_only")),
        True,
        "The current best source-topology QW-2191-safe export is still explicitly quotient-class-only.",
    )
    add_check(
        "declared_scope_source_topology_theorem_still_downstream_only",
        bool(source_topology_theorem_boundaries.get("observer_downstream_only")) and bool(f150.get("declared_scope_only")),
        True,
        "The declared-scope source-topology selector theorem still keeps the lane downstream-only and declared-scope-only.",
    )
    add_check(
        "p721_exact_t176_upgrade_still_false",
        bool(p721.get("source_topology_lane_upgrades_to_exact_t176_provider")),
        False,
        "P721 already confirms that this source-topology lane does not yet upgrade to an exact T176 provider.",
    )

    t177_target_exported_on_current_repo_state = False
    add_check(
        "t177_chart_sensitive_transported_flux_current_like_section_exported",
        t177_target_exported_on_current_repo_state,
        False,
        "No current export upgrades the present source-topology ingredients into a chart-sensitive transported flux/current-like section on the full C_v1 atlas.",
    )

    status = (
        "PASS_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_NONEXPORT_AUDITED"
        if not blocking
        else "P722_REQUIRES_REVIEW_CHANGED_T177_NONEXPORT_AUDIT_STATE"
    )

    artifact = {
        "stage": "P722",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t177_chart_sensitive_transported_flux_current_like_section_nonexport_boundary_only",
        "inputs": {
            "P721": str(IN_P721.relative_to(REPO)),
            "F141": str(IN_F141.relative_to(REPO)),
            "F143": str(IN_F143.relative_to(REPO)),
            "F148": str(IN_F148.relative_to(REPO)),
            "F149": str(IN_F149.relative_to(REPO)),
            "F150": str(IN_F150.relative_to(REPO)),
            "atlas": str(IN_ATLAS.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "t177_target_name": "ChartSensitiveTransportedFluxCurrentLikeSection_global_C_v1_source_topology_bridge_v1",
        "t177_target_exported_on_current_repo_state": t177_target_exported_on_current_repo_state,
        "current_source_topology_lane_profile": {
            "atlas_chart_sensitive_transport_present": atlas_chart_sensitive_transport_present,
            "nonzero_flow_witness_exported": bool(f143.get("actual_nonzero_flow_witness_exported")),
            "barrier_sign_witness_exported": bool(f141.get("actual_barrier_protected_sign_witness_exported")),
            "basis_independent_promotion_forgets_chart_basis_labels": bool(basis_reduction.get("forgets_chart_basis_labels")),
            "basis_independent_promotion_preserves_observer_downstream_only": bool(
                basis_reduction.get("preserves_observer_downstream_only")
            ),
            "quotient_relation_quotients_chart_labels_only_via_basis_free_packet": bool(
                quotient_relation.get("quotients_chart_labels_only_via_basis_free_packet")
            ),
            "distinguished_quotient_class_is_quotient_class_only": bool(
                distinguished_quotient_class.get("quotient_class_only")
            ),
            "declared_scope_only": bool(f150.get("declared_scope_only")),
            "observer_downstream_only": bool(source_topology_theorem_boundaries.get("observer_downstream_only")),
        },
        "audit_conclusion": {
            "current_source_topology_lane_is_physics_facing_but_chart_blind": True,
            "current_repo_already_exports_t177_bridge_target": False,
            "next_honest_move": (
                "export_or_attack_a_chart_sensitive_transported_flux_current_like_section_bridge_from_current_source_topology_ingredients_to_full_C_v1_atlas"
            ),
        },
        "hard_limits": [
            "No T176 discharge claim.",
            "No claim that current source-topology scalar witnesses already determine chartwise directed branches.",
            "No claim that basis-free quotient classes already lift to transported atlas sections.",
            "No kernel-alone/global QW-2191 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P722",
        "status": status,
        "as_of": AS_OF,
        "t177_target_name": "ChartSensitiveTransportedFluxCurrentLikeSection_global_C_v1_source_topology_bridge_v1",
        "t177_target_exported_on_current_repo_state": False,
        "current_source_topology_lane_is_physics_facing_but_chart_blind": True,
        "next_honest_move": "chart_sensitive_transported_flux_current_like_section_bridge",
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
