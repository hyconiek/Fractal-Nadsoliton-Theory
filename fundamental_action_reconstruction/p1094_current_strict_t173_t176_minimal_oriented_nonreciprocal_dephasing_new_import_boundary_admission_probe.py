#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

AS_OF = "2026-03-29"
EXPECTED_BOUNDARY = "provenance_safe_orientation_anchor_plus_nonreciprocal_rate_asymmetry_plus_oriented_memory_rule"
EXPECTED_BRIDGE = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1093 = GENERATED / "p1093_current_strict_t173_t176_candidate_oriented_nonreciprocal_dephasing_operator_provenance_audit_probe_summary.json"
IN_F968 = GENERATED / "f968_current_strict_t173_t176_candidate_oriented_nonreciprocal_dephasing_operator_provenance_packet_summary.json"

OUT_JSON = GENERATED / "p1094_current_strict_t173_t176_minimal_oriented_nonreciprocal_dephasing_new_import_boundary_admission_probe.json"
OUT_SUMMARY = GENERATED / "p1094_current_strict_t173_t176_minimal_oriented_nonreciprocal_dephasing_new_import_boundary_admission_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)
    prereq = [IN_P1093, IN_F968]
    missing = [rel(p) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P1094",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1093 = load_json(IN_P1093)
    f968 = load_json(IN_F968)

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
        "p1093_passed_as_provenance_audit",
        p1093.get("status") == "PASS_CURRENT_STRICT_T173_T176_CANDIDATE_ORIENTED_NONRECIPROCAL_DEPHASING_OPERATOR_PROVENANCE_AUDITED",
        True,
        "P1093 must pass before the minimal boundary can be frozen honestly.",
    )
    add_check(
        "f968_passed_as_provenance_packet",
        f968.get("status") == "F968_EXECUTED_CURRENT_STRICT_T173_T176_CANDIDATE_ORIENTED_NONRECIPROCAL_DEPHASING_OPERATOR_PROVENANCE_PACKET_NO_FALSE_PASS",
        True,
        "F968 must pass before the minimal boundary can be admitted as a separate object.",
    )
    add_check(
        "minimal_boundary_string_is_stable",
        p1093.get("minimal_new_import_boundary"),
        EXPECTED_BOUNDARY,
        "The irreducible new-import boundary must remain stable and explicit.",
    )
    add_check(
        "full_operator_is_still_not_exported",
        bool(p1093.get("complete_candidate_operator_exported")),
        False,
        "The full operator must remain unexported here.",
    )
    add_check(
        "boundary_is_below_lawful_supplier_export",
        bool(p1093.get("counts_as_lawful_supplier")) or bool(f968.get("counts_as_lawful_supplier")),
        False,
        "The boundary must remain below lawful supplier export.",
    )
    add_check(
        "boundary_is_below_solution_export",
        bool(f968.get("counts_as_solution")),
        False,
        "The boundary must remain below solution export.",
    )
    add_check(
        "boundary_is_below_strict_physical_orientation_datum_export",
        bool(p1093.get("counts_as_strict_physical_orientation_datum")) or bool(f968.get("counts_as_strict_physical_orientation_datum")),
        False,
        "The boundary must remain below strict physical orientation datum export.",
    )
    add_check(
        "frontier_context_bridge_is_stable",
        p1093.get("active_missing_bridge"),
        EXPECTED_BRIDGE,
        "The boundary remains framed against the same active missing bridge.",
    )

    discharged = len(blocking) == 0
    status = (
        "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ORIENTED_NONRECIPROCAL_DEPHASING_NEW_IMPORT_BOUNDARY_ADMITTED"
        if discharged
        else "P1094_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_BOUNDARY_STATE"
    )

    artifact = {
        "stage": "P1094",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "p1093_provenance_audit_summary": rel(IN_P1093),
            "f968_provenance_packet_summary": rel(IN_F968),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "exportable_boundary": {
            "object_id": "MinimalOrientedNonreciprocalDephasingNewImportBoundary_v1",
            "grade": "candidate_provider_class_seed_only",
            "orientation_anchor_component": "provenance_safe_orientation_anchor",
            "nonreciprocal_rate_asymmetry_component": "nonreciprocal_rate_asymmetry_rule",
            "oriented_memory_rule_component": "oriented_memory_rule",
            "frontier_context_bridge": EXPECTED_BRIDGE,
            "exact_reduction_to_frontier_context_bridge_exported": False,
            "admissible_as_candidate_provider_class_seed_only": discharged,
            "counts_as_lawful_supplier": False,
            "counts_as_solution": False,
            "counts_as_strict_physical_orientation_datum": False,
        },
        "classification": {
            "full_operator_exported": False,
            "minimal_new_import_boundary": EXPECTED_BOUNDARY,
            "next_honest_move": "search_one_genuinely_new_provider_class_route_from_this_boundary_without_promoting_it_to_supplier",
        },
        "hard_limits": [
            "No full operator export.",
            "No lawful supplier export.",
            "No solution export.",
            "No strict physical orientation datum export.",
            "No exact reduction to the frontier context bridge.",
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
        "boundary_object_id": artifact["exportable_boundary"]["object_id"],
        "boundary_grade": artifact["exportable_boundary"]["grade"],
        "minimal_new_import_boundary": artifact["classification"]["minimal_new_import_boundary"],
        "frontier_context_bridge": artifact["exportable_boundary"]["frontier_context_bridge"],
        "exact_reduction_to_frontier_context_bridge_exported": False,
        "admissible_as_candidate_provider_class_seed_only": artifact["exportable_boundary"]["admissible_as_candidate_provider_class_seed_only"],
        "counts_as_lawful_supplier": False,
        "counts_as_solution": False,
        "counts_as_strict_physical_orientation_datum": False,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
