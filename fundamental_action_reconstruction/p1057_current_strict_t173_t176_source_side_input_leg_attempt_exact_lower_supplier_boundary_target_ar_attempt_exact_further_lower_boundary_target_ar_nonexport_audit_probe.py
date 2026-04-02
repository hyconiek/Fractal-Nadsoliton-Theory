#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1055 = GENERATED / "p1055_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_attempt_verdict_or_exact_further_lower_boundary_nonexport_audit_probe_summary.json"
IN_P1056 = GENERATED / "p1056_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_attempt_exact_further_lower_boundary_target_probe_summary.json"
IN_P770 = GENERATED / "p770_current_strict_t224_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_attempt_probe_summary.json"
IN_T305 = ROOT / "T305_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p1057_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_ar_attempt_exact_further_lower_boundary_target_ar_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1057_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_ar_attempt_exact_further_lower_boundary_target_ar_nonexport_audit_probe_summary.json"

P1055_STATUS = (
    "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_FURTHER_LOWER_BOUNDARY_NONEXPORT_AUDITED"
)
P1056_STATUS = (
    "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_EXPORTED"
)
P770_STATUS = (
    "PASS_STRICT_T224_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
)

TARGET_NAME = (
    "source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_attempt_exact_further_lower_boundary_target_v1"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def scan_actual_t305_hits() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py")
    excluded = {
        Path(__file__).name,
        "P1057_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_AR_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_AR_NONEXPORT_AUDIT_PROBE.md",
        "N890_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_AR_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_AR_NONEXPORT_AUDIT_THEOREM.md",
        "T305_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_SPEC.md",
        "N889_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_THEOREM.md",
        "p1056_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_attempt_exact_further_lower_boundary_target_probe.py",
        "T306_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_AR_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_AR_ATTEMPT_SPEC.md",
        "N891_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_AR_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_AR_ATTEMPT_THEOREM.md",
        "p1058_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_ar_attempt_exact_further_lower_boundary_target_ar_attempt_probe.py",
    }
    hits: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if path in seen or path.name in excluded:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            if TARGET_NAME in text and "actual_realization" in text:
                hits.append(rel(path))
    return hits


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1055, IN_P1056, IN_P770, IN_T305]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1057",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1055 = load_json(IN_P1055)
    p1056 = load_json(IN_P1056)
    p770 = load_json(IN_P770)
    t305_text = load_text(IN_T305)
    actual_hits = scan_actual_t305_hits()

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

    t305_target_is_frozen = (
        p1055.get("status") == P1055_STATUS
        and p1056.get("status") == P1056_STATUS
        and p1056.get("target_name") == TARGET_NAME
        and p1056.get("target_exported_on_current_repo_state") is True
        and TARGET_NAME in t305_text
    )

    older_lower_support_attempt_remains_explicit = (
        p770.get("status") == P770_STATUS
        and p770.get("t224_attempt_exported_on_current_repo_state") is True
    )

    no_current_actual_realization_of_t305_target_found = len(actual_hits) == 0

    target_actual_realization_still_nonexport = (
        t305_target_is_frozen
        and older_lower_support_attempt_remains_explicit
        and no_current_actual_realization_of_t305_target_found
    )

    add_check(
        "t305_target_is_frozen",
        t305_target_is_frozen,
        True,
        "P1055/P1056/T305 already freeze one exact further lower-boundary target beneath T304.",
    )
    add_check(
        "older_lower_support_attempt_remains_explicit",
        older_lower_support_attempt_remains_explicit,
        True,
        "P770 still keeps one exact older route-local lower support attempt explicit.",
    )
    add_check(
        "no_current_actual_realization_of_t305_target_found",
        no_current_actual_realization_of_t305_target_found,
        True,
        "No current export yet upgrades the exact T305 target into one actual realized further lower-boundary object.",
    )
    add_check(
        "target_actual_realization_still_nonexport",
        target_actual_realization_still_nonexport,
        True,
        "Therefore the exact T305 target still remains future-only on the current repo state.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_AR_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_AR_NONEXPORT_AUDITED"
        if not blocking and target_actual_realization_still_nonexport
        else "FAIL_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_AR_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_AR_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P1057",
        "status": status,
        "as_of": AS_OF,
        "target_name": TARGET_NAME,
        "target_actual_realization_exported_on_current_repo_state": not target_actual_realization_still_nonexport,
        "current_repo_already_exports_actual_t305_target_hits": actual_hits,
        "next_honest_move_is_freeze_one_exact_actual_realization_attempt_of_t305_target": target_actual_realization_still_nonexport,
        "no_false_pass": True,
        "checks": checks,
        "blocking_checks": blocking,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "target_name": artifact["target_name"],
        "target_actual_realization_exported_on_current_repo_state": artifact[
            "target_actual_realization_exported_on_current_repo_state"
        ],
        "next_honest_move_is_freeze_one_exact_actual_realization_attempt_of_t305_target": artifact[
            "next_honest_move_is_freeze_one_exact_actual_realization_attempt_of_t305_target"
        ],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
