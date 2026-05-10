#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

AS_OF = "2026-03-29"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_V1 = ROOT / "V1_INFORMATIONAL_VISCOSITY_HYPOTHESIS_AUDIT.md"
IN_N457 = ROOT / "N457_CURRENT_FIRST_STRICT_PHASE_12_AUT_Z12_QUOTIENT_ORBIT_SPACE_TOPOLOGICAL_HOLONOMY_TRIVIALITY_BOUNDARY_THEOREM.md"
IN_NEURAL = REPO / "raport_qw540_544_neural.md"
IN_REPORT99 = REPO / "report_99_quick_win.md"
IN_TEMP = REPO / "TEMP_FULL_REPORT_CONTEXT.md"
IN_P1091 = GENERATED / "p1091_current_strict_t173_t176_post_f965_failure_map_to_exported_noncyclic_provider_split_non_same_lane_upgrade_route_decision_audit_probe_summary.json"
IN_F966 = GENERATED / "f966_current_strict_t173_t176_post_f965_failure_map_to_exported_noncyclic_provider_split_non_same_lane_upgrade_route_packet_summary.json"

OUT_JSON = GENERATED / "p1092_current_strict_t173_t176_time_arrow_light_resonance_oriented_dephasing_competing_extension_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1092_current_strict_t173_t176_time_arrow_light_resonance_oriented_dephasing_competing_extension_audit_probe_summary.json"


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)
    prereq = [IN_V1, IN_N457, IN_NEURAL, IN_REPORT99, IN_TEMP, IN_P1091, IN_F966]
    missing = [rel(p) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P1092",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    v1 = IN_V1.read_text(encoding="utf-8")
    n457 = IN_N457.read_text(encoding="utf-8")
    neural = IN_NEURAL.read_text(encoding="utf-8")
    report99 = IN_REPORT99.read_text(encoding="utf-8")
    temp = IN_TEMP.read_text(encoding="utf-8")
    p1091 = load_json(IN_P1091)
    f966 = load_json(IN_F966)

    checks: list[dict[str, Any]] = []
    blocking: list[str] = []

    def add_check(check_id: str, actual: Any, expected: Any, meaning: str) -> None:
        passed = actual == expected
        checks.append({
            "id": check_id,
            "actual": actual,
            "expected": expected,
            "pass": passed,
            "meaning": meaning,
        })
        if not passed:
            blocking.append(check_id)

    add_check(
        "v1_contains_viscosity_retardation_but_not_selector_operator",
        (
            "informational viscosity" in v1
            and "retard_phase" in v1
            and "no explicit viscosity operator" in v1
            and "selector-sector reduction" in v1
            and "a proof that viscosity-like damping breaks the residual `O(2)` degeneracy" in v1
        ),
        True,
        "V1 already keeps damping/retardation as a plausible but unsolved extension lane.",
    )
    add_check(
        "report99_contains_light_resonance_phase_material",
        (
            "Rezonans światło–materia" in report99
            and "Rabi frequency" in report99
            and "oscylacji fazy" in report99
        ),
        True,
        "Report 99 already contains non-strict light-resonance and phase-based exploratory material.",
    )
    add_check(
        "neural_report_contains_old_arrow_of_time_motif",
        (
            "Strzałka Czasu" in neural
            and "Entropia wzrosła" in neural
        ),
        True,
        "The old neural report already contains an explicit arrow-of-time motif tied to entropy growth.",
    )
    add_check(
        "temp_contains_decoherence_and_resonator_stability_motif",
        (
            "Szybka dekoherencja" in temp
            and "rezonatorów" in temp
        ),
        True,
        "The repository already contains decoherence and resonator-stability language relevant to time-sensitive measurements.",
    )
    add_check(
        "n457_blocks_pure_phase_holonomy_without_extra_structure",
        (
            "additional typed" in n457
            and "new selector slot" in n457
            and "cannot support a nontrivial Berry holonomy" in n457
        ),
        True,
        "N457 blocks pure phase/holonomy-by-itself as a strict route without extra selector structure.",
    )
    add_check(
        "f966_only_allows_new_provider_class_or_non_same_lane_upgrade_route",
        (
            p1091.get("allowed_next_move_contract")
            == "search_one_genuinely_new_non_same_lane_upgrade_route_within_exported_noncyclic_provider_split_family_or_one_genuinely_new_inversion_sensitive_source_side_provider_class"
            and f966.get("preferred_search_family") == "Xi_nad12_sigma_residual_shannon_noncyclic_provider_split_target_v1"
            and bool(f966.get("counts_as_lawful_supplier")) is False
        ),
        True,
        "The active strict route contract still allows only a genuinely new provider class or new non-same-lane upgrade route.",
    )

    discharged = len(blocking) == 0
    status = (
        "PASS_CURRENT_STRICT_T173_T176_TIME_ARROW_LIGHT_RESONANCE_ORIENTED_DEPHASING_COMPETING_EXTENSION_AUDITED"
        if discharged
        else "P1092_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_HYPOTHESIS_STATE"
    )

    artifact = {
        "stage": "P1092",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "V1": rel(IN_V1),
            "N457": rel(IN_N457),
            "raport_qw540_544_neural": rel(IN_NEURAL),
            "report_99_quick_win": rel(IN_REPORT99),
            "TEMP_FULL_REPORT_CONTEXT": rel(IN_TEMP),
            "P1091": rel(IN_P1091),
            "F966": rel(IN_F966),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "classification": {
            "repo_suggests_time_arrow_light_resonance_oriented_dephasing_lane": discharged,
            "hypothesis_grade": "competing_extension_hypothesis_only" if discharged else "unresolved",
            "closest_existing_repo_lane": "informational_viscosity_competing_extension_hypothesis",
            "preferred_candidate_form": "oriented_nonreciprocal_dephasing_or_retard_phase_operator",
            "plain_time_constant_dephasing_scalar_counts_as_symmetry_breaker": False,
            "requires_additional_oriented_transport_or_selector_slot": True,
            "counts_as_lawful_supplier": False,
            "counts_as_strict_physical_orientation_datum": False,
            "active_missing_bridge": p1091.get("active_missing_bridge"),
            "next_honest_move": "test_oriented_nonreciprocal_dephasing_or_retard_phase_operator_as_new_provider_class_candidate_not_as_current_solution",
        },
        "hard_limits": [
            "No actual oriented dephasing operator export.",
            "No lawful supplier export.",
            "No strict physical orientation datum export.",
            "No T183 discharge.",
            "No T176 discharge.",
            "No QW-2191 discharge.",
            "No ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "hypothesis_grade": artifact["classification"]["hypothesis_grade"],
        "closest_existing_repo_lane": artifact["classification"]["closest_existing_repo_lane"],
        "preferred_candidate_form": artifact["classification"]["preferred_candidate_form"],
        "plain_time_constant_dephasing_scalar_counts_as_symmetry_breaker": artifact["classification"]["plain_time_constant_dephasing_scalar_counts_as_symmetry_breaker"],
        "requires_additional_oriented_transport_or_selector_slot": artifact["classification"]["requires_additional_oriented_transport_or_selector_slot"],
        "counts_as_lawful_supplier": False,
        "counts_as_strict_physical_orientation_datum": False,
        "active_missing_bridge": artifact["classification"]["active_missing_bridge"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
