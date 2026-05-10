#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1053 = GENERATED / "p1053_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_nonexport_audit_probe_summary.json"
IN_P1054 = GENERATED / "p1054_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_attempt_probe_summary.json"
IN_P770 = GENERATED / "p770_current_strict_t224_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_attempt_probe_summary.json"
IN_T304 = ROOT / "T304_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p1055_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_attempt_verdict_or_exact_further_lower_boundary_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1055_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_attempt_verdict_or_exact_further_lower_boundary_nonexport_audit_probe_summary.json"

P1053_STATUS = (
    "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
)
P1054_STATUS = (
    "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
)
P770_STATUS = (
    "PASS_STRICT_T224_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
)

T304_ATTEMPT_NAME = (
    "W_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_attempt_v1"
)
T224_ATTEMPT_NAME = (
    "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_attempt_v1"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def scan_t304_followup_hits() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py")
    excluded = {
        Path(__file__).name,
        "P1055_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_FURTHER_LOWER_BOUNDARY_NONEXPORT_AUDIT_PROBE.md",
        "N888_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_FURTHER_LOWER_BOUNDARY_NONEXPORT_AUDIT_THEOREM.md",
        "T304_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p1054_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_attempt_probe.py",
        "p1056_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_attempt_exact_further_lower_boundary_target_probe.py",
        "T305_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_SPEC.md",
    }
    hits: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if path in seen or path.name in excluded:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            if T304_ATTEMPT_NAME not in text:
                continue
            if any(
                marker in text
                for marker in (
                    "verdict",
                    "further_lower_boundary",
                    "further lower boundary",
                    "exact further lower boundary",
                )
            ):
                hits.append(rel(path))
    return hits


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1053, IN_P1054, IN_P770, IN_T304]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1055",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1053 = load_json(IN_P1053)
    p1054 = load_json(IN_P1054)
    p770 = load_json(IN_P770)
    t304_text = load_text(IN_T304)
    followup_hits = scan_t304_followup_hits()

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

    t304_attempt_is_frozen = (
        p1053.get("status") == P1053_STATUS
        and p1054.get("status") == P1054_STATUS
        and p1054.get("attempt_name") == T304_ATTEMPT_NAME
        and p1054.get("attempt_exported_on_current_repo_state") is True
        and T304_ATTEMPT_NAME in t304_text
    )

    older_lower_support_attempt_remains_explicit = (
        p770.get("status") == P770_STATUS
        and p770.get("t224_attempt_name") == T224_ATTEMPT_NAME
        and p770.get("t224_attempt_exported_on_current_repo_state") is True
    )

    no_current_t304_verdict_or_exact_further_lower_boundary_export_found = len(followup_hits) == 0

    t304_attempt_still_has_neither_verdict_nor_exact_further_lower_boundary = (
        t304_attempt_is_frozen
        and older_lower_support_attempt_remains_explicit
        and no_current_t304_verdict_or_exact_further_lower_boundary_export_found
    )

    add_check(
        "t304_attempt_is_frozen",
        t304_attempt_is_frozen,
        True,
        "P1053/P1054/T304 already freeze one exact actual-realization attempt over the exact T303 lower supplier-boundary target.",
    )
    add_check(
        "older_lower_support_attempt_remains_explicit",
        older_lower_support_attempt_remains_explicit,
        True,
        "P770 still keeps one exact older route-local lower support attempt explicit at the same chart-label-retaining pair12 typed seed subinterface.",
    )
    add_check(
        "no_current_t304_verdict_or_exact_further_lower_boundary_export_found",
        no_current_t304_verdict_or_exact_further_lower_boundary_export_found,
        True,
        "No current export yet upgrades the exact T304 attempt into either a lawful verdict or one exact further lower boundary beneath it.",
    )
    add_check(
        "t304_attempt_still_has_neither_verdict_nor_exact_further_lower_boundary",
        t304_attempt_still_has_neither_verdict_nor_exact_further_lower_boundary,
        True,
        "Therefore the exact T304 attempt still has neither lawful verdict nor exact further lower boundary on the current repo state.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_FURTHER_LOWER_BOUNDARY_NONEXPORT_AUDITED"
        if not blocking and t304_attempt_still_has_neither_verdict_nor_exact_further_lower_boundary
        else "FAIL_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_FURTHER_LOWER_BOUNDARY_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P1055",
        "status": status,
        "as_of": AS_OF,
        "t304_attempt_name": T304_ATTEMPT_NAME,
        "t304_verdict_exported_on_current_repo_state": False,
        "t304_exact_further_lower_boundary_exported_on_current_repo_state": not no_current_t304_verdict_or_exact_further_lower_boundary_export_found,
        "current_repo_already_exports_t304_verdict_or_exact_further_lower_boundary_hits": followup_hits,
        "next_honest_move_is_freeze_exact_further_lower_boundary_target_below_t304": t304_attempt_still_has_neither_verdict_nor_exact_further_lower_boundary,
        "no_false_pass": True,
        "checks": checks,
        "blocking_checks": blocking,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "t304_attempt_name": artifact["t304_attempt_name"],
        "t304_verdict_exported_on_current_repo_state": artifact["t304_verdict_exported_on_current_repo_state"],
        "t304_exact_further_lower_boundary_exported_on_current_repo_state": artifact[
            "t304_exact_further_lower_boundary_exported_on_current_repo_state"
        ],
        "next_honest_move_is_freeze_exact_further_lower_boundary_target_below_t304": artifact[
            "next_honest_move_is_freeze_exact_further_lower_boundary_target_below_t304"
        ],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
