#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-04-27"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1134 = GENERATED / "p1134_current_strict_t173_t176_existing_t337_actual_realization_nonexport_audit_probe_summary.json"
IN_T338 = ROOT / "T338_CURRENT_STRICT_T173_T176_EXISTING_T337_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p1135_current_strict_t173_t176_existing_t337_actual_realization_attempt_export_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1135_current_strict_t173_t176_existing_t337_actual_realization_attempt_export_audit_probe_summary.json"

ACTIVE_BRIDGE = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"
T178_TARGET = "SourceTopologyToAtlasChartSeedSelectionBridge_global_C_v1_strict_v1"
T180_TARGET = "PositiveCorridorOuterInteriorChartSelectionBridge_global_C_v1_strict_v1"
TARGET_NAME = (
    "SourceTopologyToAtlasChartSeedSelectionBridge_current_best_actual_realization_attempt_exact_further_lower_boundary_target_actual_realization_attempt_exact_yet_further_lower_boundary_target_actual_realization_attempt_exact_even_further_lower_boundary_target_v1"
)
ATTEMPT_NAME = (
    "SourceTopologyToAtlasChartSeedSelectionBridge_current_best_actual_realization_attempt_exact_further_lower_boundary_target_actual_realization_attempt_exact_yet_further_lower_boundary_target_actual_realization_attempt_exact_even_further_lower_boundary_target_actual_realization_attempt_v1"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1134, IN_T338]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1135",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1134 = load_json(IN_P1134)
    t338_text = load_text(IN_T338)

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

    p1134_nonexport_audit_passed = (
        p1134.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_EXISTING_T337_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and p1134.get("active_missing_bridge") == ACTIVE_BRIDGE
        and p1134.get("target_candidate") == T178_TARGET
        and p1134.get("current_best_lower_support_family_target") == T180_TARGET
        and p1134.get("target_name") == TARGET_NAME
        and p1134.get("current_repo_has_exported_actual_realization_of_t337_target") is False
        and p1134.get("next_honest_move_is_exact_actual_realization_attempt_of_same_t337_target") is True
    )

    t338_attempt_shape_frozen = all(
        needle in t338_text
        for needle in [
            "T338_CURRENT_STRICT_T173_T176_EXISTING_T337_ACTUAL_REALIZATION_ATTEMPT_SPEC_NO_FALSE_PASS",
            ATTEMPT_NAME,
            TARGET_NAME,
            T180_TARGET,
            ACTIVE_BRIDGE,
            T178_TARGET,
            "attempt_must_not_promote_to_lawful_verdict_for_t329_by_fiat := yes",
            "attempt_must_not_promote_to_actual_t180_export_by_fiat := yes",
            "attempt_must_not_promote_to_actual_t179_export_by_fiat := yes",
            "attempt_must_not_promote_to_actual_t178_export_by_fiat := yes",
            "attempt_must_not_promote_to_lawful_supplier_by_fiat := yes",
            "attempt_must_not_promote_to_solution_by_fiat := yes",
            "attempt_must_not_promote_to_strict_physical_orientation_datum_by_fiat := yes",
        ]
    )

    attempt_exported_on_current_repo_state = (
        p1134_nonexport_audit_passed and t338_attempt_shape_frozen
    )

    add_check(
        "p1134_nonexport_audit_passed",
        p1134_nonexport_audit_passed,
        True,
        "P1134 already freezes that the exact T337 target still lacks actual realization on the current repo state.",
    )
    add_check(
        "t338_attempt_shape_frozen",
        t338_attempt_shape_frozen,
        True,
        "T338 freezes one exact first actual-realization attempt instance over the exact T337 target.",
    )
    add_check(
        "attempt_exported_on_current_repo_state",
        attempt_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one exact first actual-realization attempt over the exact T337 even-further lower-boundary target.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_EXISTING_T337_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        if not blocking and attempt_exported_on_current_repo_state
        else "FAIL_CURRENT_STRICT_T173_T176_EXISTING_T337_ACTUAL_REALIZATION_ATTEMPT_PROBE"
    )

    artifact = {
        "stage": "P1135",
        "status": status,
        "as_of": AS_OF,
        "attempt_name": ATTEMPT_NAME,
        "target_name": TARGET_NAME,
        "active_missing_bridge": ACTIVE_BRIDGE,
        "target_candidate": T178_TARGET,
        "current_best_lower_support_family_target": T180_TARGET,
        "attempt_exported_on_current_repo_state": attempt_exported_on_current_repo_state,
        "counts_as_lawful_verdict_for_t329": False,
        "counts_as_actual_t180_export": False,
        "counts_as_actual_t179_export": False,
        "counts_as_actual_t178_export": False,
        "counts_as_lawful_supplier": False,
        "counts_as_solution": False,
        "counts_as_strict_physical_orientation_datum": False,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "attempt_name": artifact["attempt_name"],
        "target_name": artifact["target_name"],
        "active_missing_bridge": artifact["active_missing_bridge"],
        "target_candidate": artifact["target_candidate"],
        "current_best_lower_support_family_target": artifact["current_best_lower_support_family_target"],
        "attempt_exported_on_current_repo_state": artifact["attempt_exported_on_current_repo_state"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
