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

IN_P742 = GENERATED / "p742_current_strict_t196_actual_source_topology_selector_witness_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_P747 = GENERATED / "p747_current_strict_t201_actual_source_topology_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge_nonexport_audit_probe_summary.json"
IN_P763 = GENERATED / "p763_current_strict_t217_pair12_source_side_branch_selection_provider_actual_realization_attempt_immediate_missing_interface_nonexport_audit_probe_summary.json"
IN_P764 = GENERATED / "p764_current_strict_t218_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_probe_summary.json"
IN_T216 = ROOT / "T216_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"
IN_T218 = ROOT / "T218_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p765_current_strict_t219_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p765_current_strict_t219_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_nonexport_audit_probe_summary.json"

T219_TARGET_NAME = (
    "Pair12SourceSideBranchSelectionProviderActualRealizationAttemptChartSensitivePair12TypedDescentInterface_strict_v1"
)
T218_TARGET_SYMBOL = (
    "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_v1"
)
EXACT_INTERFACE_NAME = (
    "chart_sensitive_pair12_typed_descent_from_Sigma_sel_src_target_v1_to_the_surviving_F301_pair12_carrier_without_Q_basis_sel_v1_terminal_collapse_and_without_projector_only_atlas_collapse"
)
CURRENT_THEOREM_FILE = (
    "N761_CURRENT_STRICT_T219_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_actual_interface_realization_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "T216_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p762_current_strict_t216_pair12_source_side_branch_selection_provider_actual_realization_attempt_probe.py",
        "N758_CURRENT_STRICT_T216_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p763_current_strict_t217_pair12_source_side_branch_selection_provider_actual_realization_attempt_immediate_missing_interface_nonexport_audit_probe.py",
        "N759_CURRENT_STRICT_T217_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_INTERFACE_NONEXPORT_AUDIT_THEOREM.md",
        "T218_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_TARGET_SPEC.md",
        "p764_current_strict_t218_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_probe.py",
        "N760_CURRENT_STRICT_T218_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_TARGET_THEOREM.md",
    }
    candidates: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if path.name in excluded_names or path in seen:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            if (
                EXACT_INTERFACE_NAME in text
                and "Sigma_sel_src_target_v1" in text
                and ("pair1/pair2" in text or "F301" in text)
            ):
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P742, IN_P747, IN_P763, IN_P764, IN_T216, IN_T218]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P765",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p742 = load_json(IN_P742)
    p747 = load_json(IN_P747)
    p763 = load_json(IN_P763)
    p764 = load_json(IN_P764)
    t216_text = load_text(IN_T216)
    t218_text = load_text(IN_T218)
    positive_actual_interface_realization_candidates = (
        scan_positive_actual_interface_realization_candidates()
    )

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

    p764_t218_target_already_exported_only_at_future_only_strength = (
        bool(p764.get("t218_target_exported_on_current_repo_state"))
        and bool(p764.get("current_t218_target_is_future_route_only"))
        and bool(p764.get("current_t218_target_freezes_exact_t216_immediate_missing_interface"))
        and bool(
            p764.get(
                "current_t218_target_remains_below_actual_interface_export_and_below_t176_discharge"
            )
        )
        and bool(
            p764.get(
                "next_honest_move_is_actual_export_of_frozen_exact_missing_interface_or_attempt_level_failure_boundary"
            )
        )
    )

    p763_exact_missing_interface_nonexport_boundary_already_exported = (
        bool(p763.get("t217_boundary_exported_on_current_repo_state"))
        and bool(p763.get("current_t216_attempt_immediate_missing_interface_is_still_unexported"))
        and bool(p763.get("current_t216_attempt_stalls_exactly_at_the_named_missing_interface"))
        and str(p763.get("exact_named_missing_interface") or "") == EXACT_INTERFACE_NAME
    )

    current_actual_selector_witness_codomain_still_lacks_actual_chart_sensitive_pair12_typed_descent_interface = (
        bool(
            p742.get(
                "current_actual_selector_witness_codomain_has_exported_basis_free_chart_label_forgetting_continuation"
            )
        )
        and bool(
            p742.get(
                "current_actual_selector_witness_codomain_continuation_remains_basis_free_not_pair12_typed"
            )
        )
        and not bool(
            p742.get(
                "current_actual_selector_witness_codomain_has_exported_pair12_typed_residual_datum_continuation"
            )
        )
    )

    current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_descent_interface = (
        bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_exported"))
        and bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe"))
        and bool(p747.get("current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas"))
        and not bool(p747.get("current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge"))
    )

    t216_t218_same_exact_interface_still_frozen = all(
        needle in t216_text
        for needle in [
            "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_v1",
            EXACT_INTERFACE_NAME,
            "attempt_immediate_missing_interface := " + EXACT_INTERFACE_NAME,
        ]
    ) and all(
        needle in t218_text
        for needle in [
            T218_TARGET_SYMBOL,
            EXACT_INTERFACE_NAME,
            "target_avoids_Q_basis_sel_v1_terminal_collapse := yes",
            "target_avoids_projector_only_atlas_collapse := yes",
        ]
    )

    current_repo_has_exported_actual_chart_sensitive_pair12_typed_descent_interface = (
        len(positive_actual_interface_realization_candidates) > 0
    )

    current_repo_still_does_not_export_actual_realization_of_t218_target = (
        p764_t218_target_already_exported_only_at_future_only_strength
        and p763_exact_missing_interface_nonexport_boundary_already_exported
        and current_actual_selector_witness_codomain_still_lacks_actual_chart_sensitive_pair12_typed_descent_interface
        and current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_descent_interface
        and t216_t218_same_exact_interface_still_frozen
        and not current_repo_has_exported_actual_chart_sensitive_pair12_typed_descent_interface
        and len(positive_actual_interface_realization_candidates) == 0
    )

    next_honest_move_is_actual_t218_interface_realization_attempt_or_attempt_level_failure_boundary = (
        current_repo_still_does_not_export_actual_realization_of_t218_target
    )

    add_check(
        "p764_t218_target_already_exported_only_at_future_only_strength",
        p764_t218_target_already_exported_only_at_future_only_strength,
        True,
        "P764 already freezes that T218 exists only at future-only target strength for the exact missing chart-sensitive pair1/pair2 typed descent interface.",
    )
    add_check(
        "p763_exact_missing_interface_nonexport_boundary_already_exported",
        p763_exact_missing_interface_nonexport_boundary_already_exported,
        True,
        "P763 already freezes that the same exact T216 missing interface remains unexported on the current repo state.",
    )
    add_check(
        "current_actual_selector_witness_codomain_still_lacks_actual_chart_sensitive_pair12_typed_descent_interface",
        current_actual_selector_witness_codomain_still_lacks_actual_chart_sensitive_pair12_typed_descent_interface,
        True,
        "The strongest current exported codomain continuation out of Sigma_sel_src_target_v1 still does not export one actual chart-sensitive pair1/pair2 typed descent interface (P742).",
    )
    add_check(
        "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_descent_interface",
        current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_descent_interface,
        True,
        "The strongest current local pair1/pair2 atlas-side lane still does not export one nonprojector chart-sensitive descent interface (P747).",
    )
    add_check(
        "t216_t218_same_exact_interface_still_frozen",
        t216_t218_same_exact_interface_still_frozen,
        True,
        "T216 and T218 still freeze the same exact missing interface on the same exact attempt route.",
    )
    add_check(
        "current_repo_has_exported_actual_chart_sensitive_pair12_typed_descent_interface",
        current_repo_has_exported_actual_chart_sensitive_pair12_typed_descent_interface,
        False,
        "No current repo artifact positively exports one actual realization of that exact chart-sensitive pair1/pair2 typed descent interface.",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t218_target",
        current_repo_still_does_not_export_actual_realization_of_t218_target,
        True,
        "Therefore the current repo still does not export one actual realization of the T218 interface target.",
    )
    add_check(
        "next_honest_move_is_actual_t218_interface_realization_attempt_or_attempt_level_failure_boundary",
        next_honest_move_is_actual_t218_interface_realization_attempt_or_attempt_level_failure_boundary,
        True,
        "Hence the next honest move is now one actual realization attempt of T218, or, only if that route later stalls, an attempt-level failure boundary for the same exact T216 attempt.",
    )

    status = (
        "PASS_STRICT_T219_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking and current_repo_still_does_not_export_actual_realization_of_t218_target
        else "FAIL_STRICT_T219_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P765",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "theorem_result": {
            "t219_target_name": T219_TARGET_NAME,
            "t219_target_exported_on_current_repo_state": False,
            "current_repo_still_does_not_export_actual_realization_of_t218_target": current_repo_still_does_not_export_actual_realization_of_t218_target,
            "current_actual_selector_witness_codomain_still_lacks_actual_chart_sensitive_pair12_typed_descent_interface": current_actual_selector_witness_codomain_still_lacks_actual_chart_sensitive_pair12_typed_descent_interface,
            "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_descent_interface": current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_descent_interface,
            "current_exact_t216_missing_interface_still_only_future_target_not_actual_export": current_repo_still_does_not_export_actual_realization_of_t218_target,
            "next_honest_move_is_actual_t218_interface_realization_attempt_or_attempt_level_failure_boundary": next_honest_move_is_actual_t218_interface_realization_attempt_or_attempt_level_failure_boundary,
            "positive_actual_interface_realization_candidates": positive_actual_interface_realization_candidates,
            "no_false_pass": True,
        },
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P765",
        "status": status,
        "as_of": AS_OF,
        "t219_target_name": artifact["theorem_result"]["t219_target_name"],
        "t219_target_exported_on_current_repo_state": artifact["theorem_result"][
            "t219_target_exported_on_current_repo_state"
        ],
        "current_repo_still_does_not_export_actual_realization_of_t218_target": artifact[
            "theorem_result"
        ]["current_repo_still_does_not_export_actual_realization_of_t218_target"],
        "current_actual_selector_witness_codomain_still_lacks_actual_chart_sensitive_pair12_typed_descent_interface": artifact[
            "theorem_result"
        ]["current_actual_selector_witness_codomain_still_lacks_actual_chart_sensitive_pair12_typed_descent_interface"],
        "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_descent_interface": artifact[
            "theorem_result"
        ]["current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_descent_interface"],
        "current_exact_t216_missing_interface_still_only_future_target_not_actual_export": artifact[
            "theorem_result"
        ]["current_exact_t216_missing_interface_still_only_future_target_not_actual_export"],
        "next_honest_move_is_actual_t218_interface_realization_attempt_or_attempt_level_failure_boundary": artifact[
            "theorem_result"
        ]["next_honest_move_is_actual_t218_interface_realization_attempt_or_attempt_level_failure_boundary"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
