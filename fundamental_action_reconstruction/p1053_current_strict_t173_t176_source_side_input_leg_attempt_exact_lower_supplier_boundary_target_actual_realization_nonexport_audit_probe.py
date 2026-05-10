#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1051 = GENERATED / "p1051_current_strict_t173_t176_source_side_input_leg_target_actual_realization_attempt_verdict_or_exact_lower_supplier_boundary_nonexport_audit_probe_summary.json"
IN_P1052 = GENERATED / "p1052_current_strict_t173_t176_source_side_input_leg_target_actual_realization_attempt_exact_lower_supplier_boundary_target_probe_summary.json"
IN_P767 = GENERATED / "p767_current_strict_t221_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_nonexport_audit_probe_summary.json"
IN_P768 = GENERATED / "p768_current_strict_t222_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_target_probe_summary.json"
IN_P769 = GENERATED / "p769_current_strict_t223_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_nonexport_audit_probe_summary.json"
IN_P770 = GENERATED / "p770_current_strict_t224_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_attempt_probe_summary.json"
IN_T303 = ROOT / "T303_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p1053_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1053_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_nonexport_audit_probe_summary.json"

P1051_STATUS = (
    "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_LOWER_SUPPLIER_BOUNDARY_NONEXPORT_AUDITED"
)
P1052_STATUS = (
    "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_EXPORTED"
)
P767_STATUS = (
    "PASS_STRICT_T221_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_NONEXPORT_AUDITED"
)
P768_STATUS = (
    "PASS_STRICT_T222_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_TARGET_EXPORTED"
)
P769_STATUS = (
    "PASS_STRICT_T223_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
)
P770_STATUS = (
    "PASS_STRICT_T224_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
)

TARGET_NAME = (
    "inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_source_side_input_leg_actual_realization_attempt_exact_lower_supplier_boundary_target_v1"
)
TARGET_SUBINTERFACE = (
    "chart_label_retaining_pair12_typed_seed_from_Sigma_sel_src_target_v1_toward_the_surviving_F301_pair12_carrier_prior_to_Q_basis_sel_v1_terminal_collapse_and_prior_to_projector_only_local_pair12_atlas_collapse"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def scan_actual_lower_supplier_boundary_hits() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py")
    excluded = {
        Path(__file__).name,
        "P1053_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
        "N886_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md",
        "T303_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_SPEC.md",
        "N885_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_THEOREM.md",
        "p1052_current_strict_t173_t176_source_side_input_leg_target_actual_realization_attempt_exact_lower_supplier_boundary_target_probe.py",
        "p1054_current_strict_t173_t176_source_side_input_leg_attempt_exact_lower_supplier_boundary_target_actual_realization_attempt_probe.py",
        "T304_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "N887_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
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

    prerequisites = [IN_P1051, IN_P1052, IN_P767, IN_P768, IN_P769, IN_P770, IN_T303]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1053",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1051 = load_json(IN_P1051)
    p1052 = load_json(IN_P1052)
    p767 = load_json(IN_P767)
    p768 = load_json(IN_P768)
    p769 = load_json(IN_P769)
    p770 = load_json(IN_P770)
    t303_text = load_text(IN_T303)
    actual_hits = scan_actual_lower_supplier_boundary_hits()

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

    t303_target_is_frozen = (
        p1051.get("status") == P1051_STATUS
        and p1052.get("status") == P1052_STATUS
        and p1052.get("target_name") == TARGET_NAME
        and p1052.get("target_exported_on_current_repo_state") is True
        and TARGET_NAME in t303_text
        and TARGET_SUBINTERFACE in t303_text
    )

    older_lower_support_family_still_explicit = (
        p767.get("status") == P767_STATUS
        and p767.get("exact_named_missing_subinterface") == TARGET_SUBINTERFACE
        and p768.get("status") == P768_STATUS
        and p768.get("current_t222_target_is_future_route_only") is True
        and p769.get("status") == P769_STATUS
        and p769.get("current_repo_still_does_not_export_actual_realization_of_t222_target") is True
        and p770.get("status") == P770_STATUS
        and p770.get("t224_attempt_exported_on_current_repo_state") is True
    )

    no_current_actual_realization_of_t303_target_found = len(actual_hits) == 0

    target_actual_realization_still_nonexport = (
        t303_target_is_frozen
        and older_lower_support_family_still_explicit
        and no_current_actual_realization_of_t303_target_found
    )

    add_check(
        "t303_target_is_frozen",
        t303_target_is_frozen,
        True,
        "P1051/P1052/T303 already freeze one exact lower supplier-boundary target beneath T302.",
    )
    add_check(
        "older_lower_support_family_still_explicit",
        older_lower_support_family_still_explicit,
        True,
        "P767/P768/P769/P770 still freeze the best current older route-local lower support family around the same exact chart-label-retaining pair12 typed seed subinterface.",
    )
    add_check(
        "no_current_actual_realization_of_t303_target_found",
        no_current_actual_realization_of_t303_target_found,
        True,
        "No current export lawfully upgrades the exact T303 target into one actual realized lower supplier-boundary.",
    )
    add_check(
        "target_actual_realization_still_nonexport",
        target_actual_realization_still_nonexport,
        True,
        "Therefore the exact T303 lower supplier-boundary target still remains future-only on the current repo state.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking and target_actual_realization_still_nonexport
        else "FAIL_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P1053",
        "status": status,
        "as_of": AS_OF,
        "target_name": TARGET_NAME,
        "target_actual_realization_exported_on_current_repo_state": not target_actual_realization_still_nonexport,
        "current_repo_already_exports_actual_lower_supplier_boundary_hits": actual_hits,
        "next_honest_move_is_freeze_one_exact_actual_realization_attempt_of_t303_target": target_actual_realization_still_nonexport,
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
        "next_honest_move_is_freeze_one_exact_actual_realization_attempt_of_t303_target": artifact[
            "next_honest_move_is_freeze_one_exact_actual_realization_attempt_of_t303_target"
        ],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
