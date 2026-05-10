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
IN_P767 = GENERATED / "p767_current_strict_t221_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_nonexport_audit_probe_summary.json"
IN_P768 = GENERATED / "p768_current_strict_t222_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_target_probe_summary.json"
IN_T220 = ROOT / "T220_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"
IN_T222 = ROOT / "T222_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p769_current_strict_t223_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p769_current_strict_t223_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_nonexport_audit_probe_summary.json"

T223_TARGET_NAME = (
    "Pair12SourceSideBranchSelectionProviderActualRealizationAttemptChartSensitivePair12TypedDescentInterfaceActualRealizationAttemptImmediateMissingSubinterface_strict_v1"
)
T222_TARGET_SYMBOL = (
    "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_target_v1"
)
EXACT_SUBINTERFACE_NAME = (
    "chart_label_retaining_pair12_typed_seed_from_Sigma_sel_src_target_v1_"
    "toward_the_surviving_F301_pair12_carrier_prior_to_Q_basis_sel_v1_terminal_"
    "collapse_and_prior_to_projector_only_local_pair12_atlas_collapse"
)
CURRENT_THEOREM_FILE = (
    "N765_CURRENT_STRICT_T223_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_actual_subinterface_realization_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "T220_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p766_current_strict_t220_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_probe.py",
        "N762_CURRENT_STRICT_T220_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p767_current_strict_t221_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_nonexport_audit_probe.py",
        "N763_CURRENT_STRICT_T221_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_NONEXPORT_AUDIT_THEOREM.md",
        "T222_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_TARGET_SPEC.md",
        "p768_current_strict_t222_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_target_probe.py",
        "N764_CURRENT_STRICT_T222_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_TARGET_THEOREM.md",
        "T224_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p770_current_strict_t224_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_attempt_probe.py",
        "N766_CURRENT_STRICT_T224_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p771_current_strict_t225_pair12_seed_slot_subsubinterface_nonexport_audit_probe.py",
        "N767_CURRENT_STRICT_T225_PAIR12_SEED_SLOT_SUBSUBINTERFACE_NONEXPORT_AUDIT_THEOREM.md",
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
                EXACT_SUBINTERFACE_NAME in text
                and "Sigma_sel_src_target_v1" in text
                and ("pair1/pair2" in text or "F301" in text)
            ):
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P742, IN_P747, IN_P767, IN_P768, IN_T220, IN_T222]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P769",
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
    p767 = load_json(IN_P767)
    p768 = load_json(IN_P768)
    t220_text = load_text(IN_T220)
    t222_text = load_text(IN_T222)
    positive_actual_subinterface_realization_candidates = (
        scan_positive_actual_subinterface_realization_candidates()
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

    p768_t222_target_already_exported_only_at_future_only_strength = (
        bool(p768.get("t222_target_exported_on_current_repo_state"))
        and bool(p768.get("current_t222_target_is_future_route_only"))
        and bool(p768.get("current_t222_target_freezes_exact_t220_immediate_missing_subinterface"))
        and bool(
            p768.get(
                "current_t222_target_remains_below_actual_subinterface_export_interface_export_and_t176_discharge"
            )
        )
        and bool(
            p768.get(
                "next_honest_move_is_actual_export_of_frozen_exact_missing_subinterface_or_exact_failure_localization_below_it"
            )
        )
    )

    p767_exact_missing_subinterface_nonexport_boundary_already_exported = (
        bool(p767.get("t221_boundary_exported_on_current_repo_state"))
        and bool(
            p767.get(
                "current_repo_still_does_not_export_actual_realization_of_t220_immediate_missing_subinterface"
            )
        )
        and bool(p767.get("current_t220_attempt_stalls_exactly_at_the_named_missing_subinterface"))
        and str(p767.get("exact_named_missing_subinterface") or "") == EXACT_SUBINTERFACE_NAME
    )

    current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_subinterface = (
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

    current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_subinterface = (
        bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_exported"))
        and bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe"))
        and bool(p747.get("current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas"))
        and not bool(p747.get("current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge"))
    )

    t220_t222_same_exact_subinterface_still_frozen = all(
        needle in t220_text
        for needle in [
            "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_v1",
            EXACT_SUBINTERFACE_NAME,
            "attempt_immediate_missing_subinterface := " + EXACT_SUBINTERFACE_NAME,
        ]
    ) and all(
        needle in t222_text
        for needle in [
            T222_TARGET_SYMBOL,
            EXACT_SUBINTERFACE_NAME,
            "target_precedes_Q_basis_sel_v1_terminal_collapse := yes",
            "target_precedes_projector_only_local_pair12_atlas_collapse := yes",
        ]
    )

    current_repo_has_exported_actual_chart_label_retaining_pair12_typed_seed_subinterface = (
        len(positive_actual_subinterface_realization_candidates) > 0
    )

    current_repo_still_does_not_export_actual_realization_of_t222_target = (
        p768_t222_target_already_exported_only_at_future_only_strength
        and p767_exact_missing_subinterface_nonexport_boundary_already_exported
        and current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_subinterface
        and current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_subinterface
        and t220_t222_same_exact_subinterface_still_frozen
        and not current_repo_has_exported_actual_chart_label_retaining_pair12_typed_seed_subinterface
        and len(positive_actual_subinterface_realization_candidates) == 0
    )

    next_honest_move_is_actual_t222_subinterface_realization_attempt_or_exact_failure_localization_below_it = (
        current_repo_still_does_not_export_actual_realization_of_t222_target
    )

    add_check(
        "p768_t222_target_already_exported_only_at_future_only_strength",
        p768_t222_target_already_exported_only_at_future_only_strength,
        True,
        "P768 already freezes that T222 exists only at future-only target strength for the exact missing chart-label-retaining pair1/pair2 typed seed subinterface.",
    )
    add_check(
        "p767_exact_missing_subinterface_nonexport_boundary_already_exported",
        p767_exact_missing_subinterface_nonexport_boundary_already_exported,
        True,
        "P767 already freezes that the same exact T220 missing subinterface remains unexported on the current repo state.",
    )
    add_check(
        "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_subinterface",
        current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_subinterface,
        True,
        "The strongest current exported codomain continuation out of Sigma_sel_src_target_v1 still does not export one actual chart-label-retaining pair1/pair2 typed seed subinterface (P742).",
    )
    add_check(
        "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_subinterface",
        current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_subinterface,
        True,
        "The strongest current local pair1/pair2 atlas-side lane still does not export one nonprojector chart-label-retaining seed subinterface (P747).",
    )
    add_check(
        "t220_t222_same_exact_subinterface_still_frozen",
        t220_t222_same_exact_subinterface_still_frozen,
        True,
        "T220 and T222 still freeze the same exact missing seed subinterface on the same exact attempt route.",
    )
    add_check(
        "current_repo_has_exported_actual_chart_label_retaining_pair12_typed_seed_subinterface",
        current_repo_has_exported_actual_chart_label_retaining_pair12_typed_seed_subinterface,
        False,
        "No current repo artifact positively exports one actual realization of that exact chart-label-retaining pair1/pair2 typed seed subinterface.",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t222_target",
        current_repo_still_does_not_export_actual_realization_of_t222_target,
        True,
        "Therefore the current repo still does not export one actual realization of the T222 subinterface target.",
    )
    add_check(
        "next_honest_move_is_actual_t222_subinterface_realization_attempt_or_exact_failure_localization_below_it",
        next_honest_move_is_actual_t222_subinterface_realization_attempt_or_exact_failure_localization_below_it,
        True,
        "Hence the next honest move is now one actual realization attempt of T222, or, only if that same route stalls later, exact failure localization below the same subinterface.",
    )

    status = (
        "PASS_STRICT_T223_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking and current_repo_still_does_not_export_actual_realization_of_t222_target
        else "FAIL_STRICT_T223_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P769",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "theorem_result": {
            "t223_target_name": T223_TARGET_NAME,
            "t223_target_exported_on_current_repo_state": False,
            "current_repo_still_does_not_export_actual_realization_of_t222_target": current_repo_still_does_not_export_actual_realization_of_t222_target,
            "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_subinterface": current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_subinterface,
            "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_subinterface": current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_subinterface,
            "current_exact_t220_missing_subinterface_still_only_future_target_not_actual_export": current_repo_still_does_not_export_actual_realization_of_t222_target,
            "next_honest_move_is_actual_t222_subinterface_realization_attempt_or_exact_failure_localization_below_it": next_honest_move_is_actual_t222_subinterface_realization_attempt_or_exact_failure_localization_below_it,
            "no_false_pass": True,
        },
        "positive_actual_subinterface_realization_candidates": positive_actual_subinterface_realization_candidates,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P769",
        "status": status,
        "as_of": AS_OF,
        "t223_target_name": artifact["theorem_result"]["t223_target_name"],
        "t223_target_exported_on_current_repo_state": artifact["theorem_result"][
            "t223_target_exported_on_current_repo_state"
        ],
        "current_repo_still_does_not_export_actual_realization_of_t222_target": artifact[
            "theorem_result"
        ]["current_repo_still_does_not_export_actual_realization_of_t222_target"],
        "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_subinterface": artifact[
            "theorem_result"
        ][
            "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_subinterface"
        ],
        "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_subinterface": artifact[
            "theorem_result"
        ]["current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_subinterface"],
        "current_exact_t220_missing_subinterface_still_only_future_target_not_actual_export": artifact[
            "theorem_result"
        ]["current_exact_t220_missing_subinterface_still_only_future_target_not_actual_export"],
        "next_honest_move_is_actual_t222_subinterface_realization_attempt_or_exact_failure_localization_below_it": artifact[
            "theorem_result"
        ]["next_honest_move_is_actual_t222_subinterface_realization_attempt_or_exact_failure_localization_below_it"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
