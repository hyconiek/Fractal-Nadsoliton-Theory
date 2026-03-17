#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-17"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_N8 = GENERATED / "n8_current_strict_core_sigma_int_residual_datum_obstruction_after_target_slot_export_theorem_summary.json"
IN_P5 = GENERATED / "p5_strict_core_sigma_int_to_residual_datum_rerun_after_target_slot_export_summary.json"
IN_P451 = GENERATED / "p451_current_strict_sigma_int_slot_free_r1_target_slot_population_probe_summary.json"
IN_N680 = GENERATED / "n680_current_strict_t173_projective_strict_core_selector_closure_discharge_theorem_summary.json"
IN_N681 = GENERATED / "n681_current_strict_t173_directed_output_sign_lift_obstruction_boundary_theorem_summary.json"
IN_N686 = GENERATED / "n686_current_strict_t173_global_axis_only_transition_edge_sign_flip_boundary_theorem_summary.json"
IN_N687 = GENERATED / "n687_current_strict_t173_global_edge_sign_coherence_obstruction_boundary_theorem_summary.json"
IN_P709 = GENERATED / "p709_current_strict_release_7_os_residual_sign_gauge_irrelevance_audit_probe_summary.json"

OUT_JSON = GENERATED / "p711_current_strict_t173_previous_methodology_survival_and_global_gap_audit_probe.json"
OUT_SUMMARY = GENERATED / "p711_current_strict_t173_previous_methodology_survival_and_global_gap_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def theorem_result(obj: dict[str, Any]) -> dict[str, Any]:
    val = obj.get("theorem_result")
    return val if isinstance(val, dict) else {}


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_N8, IN_P5, IN_P451, IN_N680, IN_N681, IN_N686, IN_N687, IN_P709]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P711",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    n8 = load_json(IN_N8)
    p5 = load_json(IN_P5)
    p451 = load_json(IN_P451)
    n680 = load_json(IN_N680)
    n681 = load_json(IN_N681)
    n686 = load_json(IN_N686)
    n687 = load_json(IN_N687)
    p709 = load_json(IN_P709)

    n8_tr = theorem_result(n8)
    n680_tr = theorem_result(n680)
    n681_tr = theorem_result(n681)
    n686_tr = theorem_result(n686)
    n687_tr = theorem_result(n687)

    checks: list[dict[str, Any]] = []
    blocking: list[str] = []

    def add_check(check_id: str, actual: Any, expected: Any, meaning: str) -> None:
        ok = actual == expected
        checks.append(
            {
                "id": check_id,
                "actual": actual,
                "expected": expected,
                "pass": ok,
                "meaning": meaning,
            }
        )
        if not ok:
            blocking.append(check_id)

    add_check(
        "n8_local_sigma_int_bridge_route_discharged",
        bool(n8_tr.get("discharged")),
        True,
        "The previous topological sigma_int methodology survives locally: theorem-level bridge discharge exists on the strict R1 lane (N8).",
    )
    add_check(
        "n8_target_slot_population_derived",
        bool(n8_tr.get("strict_core_sigma_int_route_derives_target_slot_population")),
        True,
        "The strict sigma_int lane derives target-slot population in declared R1 scope (N8).",
    )
    add_check(
        "p5_route_computable",
        p5.get("status"),
        "PASS_COMPUTABLE_FROM_CURRENT_STRICT_CORE_RESIDUAL_DATUM_ROUTE_AFTER_TARGET_SLOT_EXPORT",
        "The rerun probe confirms the strict sigma_int residual-datum route is computable up to post-map object support (P5).",
    )
    add_check(
        "p451_slot_free_theta_population_audited",
        bool(p451.get("audits_pass")) and bool(p451.get("strict_core_promotion")),
        True,
        "The slot-free sigma_int theta-pair source yields an audited strict-core R1 inhabitant instance (P451).",
    )
    add_check(
        "n680_projective_global_closure_present",
        bool(n680_tr.get("strict_core_selector_closure")),
        True,
        "Projective global strict-core selector closure is exported (N680).",
    )
    add_check(
        "n680_kernel_alone_global_qw2191_still_open",
        bool(n680_tr.get("QW2191_kernel_alone_discharge")),
        False,
        "Even after those local/topological survivals, kernel-alone/global QW-2191 discharge remains false (N680).",
    )
    add_check(
        "n681_no_directed_output_sign_lift_in_strict_core",
        bool(n681_tr.get("directed_output_sign_lift_determined_in_strict_core")),
        False,
        "The previous methodology still does not determine a directed output sign-lift in strict core (N681).",
    )
    add_check(
        "n686_global_axis_only_edge_sign_flips_present",
        bool(n686_tr.get("global_axis_only_transition_edge_sign_flips_present")),
        True,
        "At global atlas level, axis-only transition representatives still force sign flips on some edges (N686).",
    )
    add_check(
        "n687_no_chart_sign_relift_solution",
        bool(n687_tr.get("global_edge_sign_coherence_solvable_by_chart_sign_relift")),
        False,
        "No per-chart Z2 relift solves the full-edge sign-coherence problem under fixed axis-only transitions (N687).",
    )
    add_check(
        "p709_release7_os_sign_gauge_irrelevance",
        p709.get("status"),
        "PASS_RELEASE_7_OS_RESIDUAL_SIGN_GAUGE_IRRELEVANCE_AUDITED",
        "Operationally, the Release-7 OS outputs used downstream are already audited as residual-sign gauge-irrelevant (P709).",
    )

    status = (
        "PASS_PREVIOUS_METHODOLOGY_SURVIVAL_AND_GLOBAL_GAP_AUDITED"
        if not blocking
        else "P711_REQUIRES_REVIEW_CHANGED_OR_INCOMPLETE_PREVIOUS_METHODOLOGY_AUDIT_STATE"
    )

    artifact = {
        "stage": "P711",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t173_previous_methodology_survival_and_global_gap_only",
        "inputs": {
            "N8": str(IN_N8.relative_to(REPO)),
            "P5": str(IN_P5.relative_to(REPO)),
            "P451": str(IN_P451.relative_to(REPO)),
            "N680": str(IN_N680.relative_to(REPO)),
            "N681": str(IN_N681.relative_to(REPO)),
            "N686": str(IN_N686.relative_to(REPO)),
            "N687": str(IN_N687.relative_to(REPO)),
            "P709": str(IN_P709.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "surviving_previous_methodology_components": {
            "strict_sigma_int_local_bridge_route": bool(
                n8_tr.get("discharged")
                and n8_tr.get("strict_core_sigma_int_route_derives_target_slot_population")
                and n8_tr.get("strict_core_sigma_int_route_derives_actual_object_support_above_export_map_object")
            ),
            "strict_r1_target_slot_population_from_slot_free_theta_source": bool(
                p451.get("audits_pass") and p451.get("strict_core_promotion")
            ),
            "projective_global_selector_closure": bool(n680_tr.get("strict_core_selector_closure")),
            "release_7_os_residual_sign_gauge_irrelevance": (
                p709.get("status") == "PASS_RELEASE_7_OS_RESIDUAL_SIGN_GAUGE_IRRELEVANCE_AUDITED"
            ),
        },
        "global_gap_after_previous_methodology_review": {
            "kernel_alone_global_qw2191_discharge": bool(n680_tr.get("QW2191_kernel_alone_discharge")),
            "directed_output_sign_lift_determined_in_strict_core": bool(
                n681_tr.get("directed_output_sign_lift_determined_in_strict_core")
            ),
            "global_axis_only_transition_edge_sign_flips_present": bool(
                n686_tr.get("global_axis_only_transition_edge_sign_flips_present")
            ),
            "global_edge_sign_coherence_solvable_by_chart_sign_relift": bool(
                n687_tr.get("global_edge_sign_coherence_solvable_by_chart_sign_relift")
            ),
        },
        "audit_conclusion": {
            "previous_methodology_contains_reusable_strict_ingredients": True,
            "previous_methodology_suffices_for_global_t173_discharge": False,
            "small_gap_identified": (
                "missing_global_provider_from_local_sigma_int_or_pair_scoped_selector_ingredients_to_chartwise_"
                "directed_sign_coherence_on_the_full_C_v1_atlas"
            ),
            "operational_release_7_continuation_can_remain_honest": True,
            "recommended_next_strict_move": (
                "seek_genuinely_new_global_provider_if_directed_physical_sign_datum_is_still_required;"
                "_otherwise_keep_residual_sign_frozen_as_gauge_for_release_7_os_scope"
            ),
        },
        "hard_limits": [
            "No kernel-alone/global QW-2191 discharge.",
            "No directed/sign-sensitive physical orientation datum promotion into strict core.",
            "No Standard Model identification claim.",
            "No strict physical-unit calibration claim.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P711",
        "status": status,
        "as_of": AS_OF,
        "previous_methodology_contains_reusable_strict_ingredients": True,
        "local_strict_sigma_int_bridge_route_available": bool(
            n8_tr.get("discharged") and n8_tr.get("strict_core_sigma_int_route_derives_target_slot_population")
        ),
        "previous_methodology_suffices_for_global_t173_discharge": False,
        "small_gap": "missing_global_provider_for_chartwise_directed_sign_coherence_on_full_C_v1_atlas",
        "release_7_os_scope_can_keep_residual_sign_as_gauge": (
            p709.get("status") == "PASS_RELEASE_7_OS_RESIDUAL_SIGN_GAUGE_IRRELEVANCE_AUDITED"
        ),
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
