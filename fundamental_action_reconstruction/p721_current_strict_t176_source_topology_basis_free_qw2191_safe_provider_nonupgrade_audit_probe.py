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

IN_P229 = GENERATED / "p229_current_actual_source_topology_barrier_protected_sign_witness_probe_summary.json"
IN_P231 = GENERATED / "p231_current_actual_source_topology_nonzero_flow_witness_probe_summary.json"
IN_P235 = GENERATED / "p235_current_actual_source_topology_selector_witness_probe_summary.json"
IN_P236 = GENERATED / "p236_current_actual_source_topology_basis_independent_promotion_witness_probe_summary.json"
IN_P237 = GENERATED / "p237_current_actual_source_topology_quotient_safe_qw2191_resolution_witness_probe_summary.json"
IN_F141 = GENERATED / "f141_first_actual_source_topology_barrier_protected_sign_witness_packet_summary.json"
IN_F143 = GENERATED / "f143_first_actual_source_topology_nonzero_flow_witness_packet_summary.json"
IN_F147 = GENERATED / "f147_first_actual_source_topology_selector_witness_packet_summary.json"
IN_F148 = GENERATED / "f148_first_actual_source_topology_basis_independent_promotion_witness_packet_summary.json"
IN_F149 = GENERATED / "f149_first_actual_source_topology_quotient_safe_qw2191_resolution_witness_packet_summary.json"
IN_P720 = GENERATED / "p720_current_strict_t176_observer_facing_signed_output_channel_projection_provider_class_audit_probe_summary.json"

OUT_JSON = GENERATED / "p721_current_strict_t176_source_topology_basis_free_qw2191_safe_provider_nonupgrade_audit_probe.json"
OUT_SUMMARY = GENERATED / "p721_current_strict_t176_source_topology_basis_free_qw2191_safe_provider_nonupgrade_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [
        IN_P229,
        IN_P231,
        IN_P235,
        IN_P236,
        IN_P237,
        IN_F141,
        IN_F143,
        IN_F147,
        IN_F148,
        IN_F149,
        IN_P720,
    ]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P721",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p229 = load_json(IN_P229)
    p231 = load_json(IN_P231)
    p235 = load_json(IN_P235)
    p236 = load_json(IN_P236)
    p237 = load_json(IN_P237)
    f141 = load_json(IN_F141)
    f143 = load_json(IN_F143)
    f147 = load_json(IN_F147)
    f148 = load_json(IN_F148)
    f149 = load_json(IN_F149)
    p720 = load_json(IN_P720)

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

    add_check(
        "actual_nonzero_flow_witness_exported",
        bool(p231.get("actual_nonzero_flow_witness_exported")) and bool(f143.get("actual_nonzero_flow_witness_exported")),
        True,
        "The physics-facing source-topology lane does export an actual nonzero-flow witness.",
    )
    add_check(
        "actual_barrier_protected_sign_witness_exported",
        bool(p229.get("actual_barrier_protected_sign_witness_exported")) and bool(f141.get("actual_barrier_protected_sign_witness_exported")),
        True,
        "The same lane exports an actual barrier-protected sign witness.",
    )
    add_check(
        "actual_selector_witness_exported",
        bool(p235.get("actual_selector_witness_exported")) and bool(f147.get("actual_selector_witness_exported")),
        True,
        "The same lane exports an actual selector witness.",
    )
    add_check(
        "basis_independent_promotion_discharged",
        bool(p236.get("basis_independence_discharged")) and bool(f148.get("basis_independence_discharged")),
        True,
        "The source-topology selector lane reaches a basis-independent promotion witness.",
    )
    add_check(
        "qw2191_quotient_safe_resolution_discharged",
        bool(p237.get("qw2191_quotient_safe_discharged")) and bool(f149.get("qw2191_quotient_safe_discharged")),
        True,
        "That same source-topology lane reaches a quotient-safe QW-2191 resolution witness.",
    )
    add_check(
        "source_topology_lane_observer_role_still_downstream_only",
        {
            "P229": p229.get("observer_role"),
            "P231": p231.get("observer_role"),
            "P235": p235.get("observer_role"),
            "P236": p236.get("observer_role"),
            "P237": p237.get("observer_role"),
        },
        {
            "P229": "downstream_only",
            "P231": "downstream_only",
            "P235": "downstream_only",
            "P236": "downstream_only",
            "P237": "downstream_only",
        },
        "All current source-topology witness exports remain explicitly downstream-only rather than promoted into the current strict selector closure carrier.",
    )
    add_check(
        "tau_src_not_identified_with_current_selector_carrier",
        {
            "P235": bool(p235.get("tau_src_identified_with_s_prelm")),
            "P236": bool(p236.get("tau_src_identified_with_s_prelm")),
            "P237": bool(p237.get("tau_src_identified_with_s_prelm")),
        },
        {
            "P235": False,
            "P236": False,
            "P237": False,
        },
        "The source-topology input tau_src is still not identified with the current strict selector carrier.",
    )
    add_check(
        "basis_free_qw2191_safe_output_is_quotient_class_only",
        {
            "distinguished_quotient_class_exported": bool(p237.get("distinguished_quotient_class_exported")),
            "raw_theta_uniqueness_claimed": bool(p237.get("raw_theta_uniqueness_claimed")),
            "quotient_class_only": bool(((f149.get("support_packet") or {}).get("distinguished_quotient_class") or {}).get("quotient_class_only")),
        },
        {
            "distinguished_quotient_class_exported": True,
            "raw_theta_uniqueness_claimed": False,
            "quotient_class_only": True,
        },
        "The current source-topology QW-2191-safe result remains a quotient-class export only, not one exact directed branch.",
    )
    add_check(
        "source_topology_lane_current_selector_closure_false",
        {
            "P236": bool(p236.get("current_selector_closure")),
            "P237": bool(p237.get("current_selector_closure")),
        },
        {
            "P236": False,
            "P237": False,
        },
        "Even after basis-free promotion and quotient-safe resolution, the source-topology lane still does not export current strict selector closure.",
    )
    add_check(
        "observer_facing_static_output_class_also_not_exact",
        bool(p720.get("exact_candidates_found")),
        False,
        "The nearest observer-facing static output-channel provider class is also already audited as non-exact, so there is no hidden exact upgrade on the already exported output carrier.",
    )

    status = (
        "PASS_SOURCE_TOPOLOGY_BASIS_FREE_QW2191_SAFE_PROVIDER_NONUPGRADE_AUDITED"
        if not blocking
        else "P721_REQUIRES_REVIEW_CHANGED_SOURCE_TOPOLOGY_PROVIDER_STATE"
    )

    artifact = {
        "stage": "P721",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": (
            "audit whether the current physics-facing source-topology lane (nonzero flow, barrier sign, selector witness, "
            "basis-free promotion, quotient-safe QW2191 resolution) already upgrades to an exact T176 provider on C_v1"
        ),
        "inputs": {
            "P229": str(IN_P229.relative_to(REPO)),
            "P231": str(IN_P231.relative_to(REPO)),
            "P235": str(IN_P235.relative_to(REPO)),
            "P236": str(IN_P236.relative_to(REPO)),
            "P237": str(IN_P237.relative_to(REPO)),
            "F141": str(IN_F141.relative_to(REPO)),
            "F143": str(IN_F143.relative_to(REPO)),
            "F147": str(IN_F147.relative_to(REPO)),
            "F148": str(IN_F148.relative_to(REPO)),
            "F149": str(IN_F149.relative_to(REPO)),
            "P720": str(IN_P720.relative_to(REPO)),
        },
        "current_source_topology_provider_lane_state": {
            "actual_nonzero_flow_witness_exported": bool(p231.get("actual_nonzero_flow_witness_exported")),
            "actual_barrier_protected_sign_witness_exported": bool(p229.get("actual_barrier_protected_sign_witness_exported")),
            "actual_selector_witness_exported": bool(p235.get("actual_selector_witness_exported")),
            "basis_independence_discharged": bool(p236.get("basis_independence_discharged")),
            "qw2191_quotient_safe_discharged": bool(p237.get("qw2191_quotient_safe_discharged")),
            "observer_role": {
                "P229": p229.get("observer_role"),
                "P231": p231.get("observer_role"),
                "P235": p235.get("observer_role"),
                "P236": p236.get("observer_role"),
                "P237": p237.get("observer_role"),
            },
            "tau_src_identified_with_s_prelm": {
                "P235": bool(p235.get("tau_src_identified_with_s_prelm")),
                "P236": bool(p236.get("tau_src_identified_with_s_prelm")),
                "P237": bool(p237.get("tau_src_identified_with_s_prelm")),
            },
            "distinguished_quotient_class_exported": bool(p237.get("distinguished_quotient_class_exported")),
            "raw_theta_uniqueness_claimed": bool(p237.get("raw_theta_uniqueness_claimed")),
            "quotient_class_only": bool(((f149.get("support_packet") or {}).get("distinguished_quotient_class") or {}).get("quotient_class_only")),
            "current_selector_closure": {
                "P236": bool(p236.get("current_selector_closure")),
                "P237": bool(p237.get("current_selector_closure")),
            },
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "result": {
            "source_topology_basis_free_qw2191_safe_lane_contains_physically_interpretable_strict_ingredients": True,
            "source_topology_basis_free_qw2191_safe_lane_upgrades_to_exact_t176_provider": False,
            "current_best_source_topology_output_is_basis_free_quotient_class_only": True,
            "missing_upgrade_step": (
                "bridge_from_downstream_basis_free_quotient_class_to_chartwise_transported_exact_directed_sign_section_on_full_C_v1_atlas"
            ),
            "preferred_next_direction": (
                "chart_sensitive_transported_flux_or_current_like_provider_class_over_further_static_basis_free_or_output_axis_classes"
            ),
        },
        "hard_limits": [
            "No T176 discharge claim.",
            "No claim that source-topology witnesses already populate the current strict selector carrier.",
            "No claim that quotient-safe source-topology resolution determines one exact directed branch.",
            "No kernel-alone/global QW-2191 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P721",
        "status": status,
        "as_of": AS_OF,
        "source_topology_physically_interpretable_strict_ingredients_present": True,
        "basis_independence_discharged": bool(p236.get("basis_independence_discharged")),
        "qw2191_quotient_safe_discharged": bool(p237.get("qw2191_quotient_safe_discharged")),
        "source_topology_lane_upgrades_to_exact_t176_provider": False,
        "current_best_source_topology_output_is_basis_free_quotient_class_only": True,
        "preferred_next_direction": (
            "chart_sensitive_transported_flux_or_current_like_provider_class_over_further_static_basis_free_or_output_axis_classes"
        ),
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
