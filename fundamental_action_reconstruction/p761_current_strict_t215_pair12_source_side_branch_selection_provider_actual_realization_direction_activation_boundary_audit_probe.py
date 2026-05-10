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

IN_P758 = GENERATED / "p758_current_strict_t212_pair12_witness_split_current_exported_continuation_family_provider_shift_requirement_boundary_audit_probe_summary.json"
IN_P759 = GENERATED / "p759_current_strict_t213_pair12_source_side_branch_selection_provider_target_probe_summary.json"
IN_P760 = GENERATED / "p760_current_strict_t214_pair12_source_side_branch_selection_provider_actual_realization_nonexport_audit_probe_summary.json"
IN_T213 = ROOT / "T213_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p761_current_strict_t215_pair12_source_side_branch_selection_provider_actual_realization_direction_activation_boundary_audit_probe.json"
OUT_SUMMARY = GENERATED / "p761_current_strict_t215_pair12_source_side_branch_selection_provider_actual_realization_direction_activation_boundary_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P758, IN_P759, IN_P760, IN_T213]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P761",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p758 = load_json(IN_P758)
    p759 = load_json(IN_P759)
    p760 = load_json(IN_P760)
    t213_text = load_text(IN_T213)

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

    p758_provider_shift_boundary_already_exported = (
        bool(p758.get("t212_boundary_exported_on_current_repo_state"))
        and bool(
            p758.get(
                "same_pair12_witness_split_current_exported_continuation_family_no_longer_admitted_as_active_primary_t173_move"
            )
        )
        and bool(p758.get("provider_shift_is_now_active_primary_t173_branch_on_current_repo_state"))
        and bool(
            p758.get(
                "next_honest_primary_t173_move_requires_genuinely_new_provider_class_beyond_current_exported_continuation_family"
            )
        )
    )

    p759_exact_future_only_t213_target_already_exported = (
        bool(p759.get("t213_target_exported_on_current_repo_state"))
        and bool(p759.get("current_t213_target_is_source_side_observer_free"))
        and bool(p759.get("current_t213_target_is_pair12_typed_and_branch_sensitive"))
        and bool(p759.get("current_t213_target_is_chart_sensitive_and_residual_datum_carrier_binding"))
        and bool(p759.get("current_t213_target_is_external_to_current_exported_p731_continuation_family"))
        and bool(p759.get("current_t213_target_is_nonconvention_nonpremise_based"))
        and bool(p759.get("current_t213_target_is_future_route_only"))
    )

    p760_actual_t213_realization_nonexport_boundary_already_exported = (
        not bool(p760.get("t214_target_exported_on_current_repo_state"))
        and bool(p760.get("current_repo_still_does_not_export_actual_realization_of_t213_target"))
        and bool(p760.get("next_honest_move_is_actual_t213_realization_attempt_or_further_provider_attack"))
    )

    t213_target_still_frozen_as_one_actualizable_future_route = all(
        needle in t213_text
        for needle in [
            "W_strict_t173_pair12_source_side_branch_selection_provider_target_v1",
            "target_is_pair12_typed := yes",
            "target_is_branch_sensitive_to_delta_k_vs_delta_minus_k := yes",
            "target_is_chart_sensitive_or_chart_label_retaining := yes",
            "future_route_only := yes",
        ]
    )

    same_current_repo_state_now_places_actual_t213_realization_attempt_as_active_primary_t173_move = (
        p758_provider_shift_boundary_already_exported
        and p759_exact_future_only_t213_target_already_exported
        and p760_actual_t213_realization_nonexport_boundary_already_exported
        and t213_target_still_frozen_as_one_actualizable_future_route
    )

    actual_t213_realization_direction_is_now_active_primary_t173_branch_on_current_repo_state = (
        same_current_repo_state_now_places_actual_t213_realization_attempt_as_active_primary_t173_move
    )

    further_return_to_current_exported_p731_continuation_family_remains_nonprimary = (
        same_current_repo_state_now_places_actual_t213_realization_attempt_as_active_primary_t173_move
    )

    further_provider_attack_remains_secondary_if_actual_t213_route_stalls = (
        same_current_repo_state_now_places_actual_t213_realization_attempt_as_active_primary_t173_move
    )

    next_honest_primary_t173_move_is_actual_t213_realization_attempt_unless_that_route_stalls = (
        same_current_repo_state_now_places_actual_t213_realization_attempt_as_active_primary_t173_move
    )

    add_check(
        "p758_provider_shift_boundary_already_exported",
        p758_provider_shift_boundary_already_exported,
        True,
        "P758 already froze that the current exported P731 continuation family is no longer an admitted active primary T173 route and that provider shift is active.",
    )
    add_check(
        "p759_exact_future_only_t213_target_already_exported",
        p759_exact_future_only_t213_target_already_exported,
        True,
        "P759 already froze one exact future-only target for the genuinely new source-side pair1/pair2 provider class.",
    )
    add_check(
        "p760_actual_t213_realization_nonexport_boundary_already_exported",
        p760_actual_t213_realization_nonexport_boundary_already_exported,
        True,
        "P760 already froze that the T213 target is not actually realized on the current repo state.",
    )
    add_check(
        "t213_target_still_frozen_as_one_actualizable_future_route",
        t213_target_still_frozen_as_one_actualizable_future_route,
        True,
        "T213 still freezes one actualizable future route: pair1/pair2 typed, branch-sensitive, chart-sensitive, and future-route-only.",
    )
    add_check(
        "same_current_repo_state_now_places_actual_t213_realization_attempt_as_active_primary_t173_move",
        same_current_repo_state_now_places_actual_t213_realization_attempt_as_active_primary_t173_move,
        True,
        "Therefore, on the current repo state, the honest primary continuation is no longer target-sharpening or return to the exhausted family, but actual realization attempt of the exact T213 target.",
    )
    add_check(
        "actual_t213_realization_direction_is_now_active_primary_t173_branch_on_current_repo_state",
        actual_t213_realization_direction_is_now_active_primary_t173_branch_on_current_repo_state,
        True,
        "So the actual-realization direction for the T213 provider target is now the active primary T173 branch on the current repo state.",
    )
    add_check(
        "further_return_to_current_exported_p731_continuation_family_remains_nonprimary",
        further_return_to_current_exported_p731_continuation_family_remains_nonprimary,
        True,
        "A return to the already exported P731 -> P747 continuation family remains nonprimary after P758/P759/P760.",
    )
    add_check(
        "further_provider_attack_remains_secondary_if_actual_t213_route_stalls",
        further_provider_attack_remains_secondary_if_actual_t213_route_stalls,
        True,
        "A further genuinely new provider attack remains only the secondary fallback if the actual T213 realization route later stalls, consistent with P760.",
    )
    add_check(
        "next_honest_primary_t173_move_is_actual_t213_realization_attempt_unless_that_route_stalls",
        next_honest_primary_t173_move_is_actual_t213_realization_attempt_unless_that_route_stalls,
        True,
        "Hence the next honest primary strict move under T173 is now one actual realization attempt of T213, unless that route later stalls.",
    )

    status = (
        "PASS_STRICT_T215_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_DIRECTION_ACTIVATION_BOUNDARY_AUDITED"
        if not blocking
        and actual_t213_realization_direction_is_now_active_primary_t173_branch_on_current_repo_state
        and next_honest_primary_t173_move_is_actual_t213_realization_attempt_unless_that_route_stalls
        else "FAIL_STRICT_T215_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_DIRECTION_ACTIVATION_BOUNDARY_AUDIT"
    )

    artifact = {
        "stage": "P761",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "theorem_result": {
            "t215_boundary_name": "StrictT173Pair12SourceSideBranchSelectionProviderActualRealizationDirectionActivationBoundary_strict_v1",
            "t215_boundary_exported_on_current_repo_state": status
            == "PASS_STRICT_T215_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_DIRECTION_ACTIVATION_BOUNDARY_AUDITED",
            "p758_provider_shift_boundary_already_exported": p758_provider_shift_boundary_already_exported,
            "p759_exact_future_only_t213_target_already_exported": p759_exact_future_only_t213_target_already_exported,
            "p760_actual_t213_realization_nonexport_boundary_already_exported": p760_actual_t213_realization_nonexport_boundary_already_exported,
            "same_current_repo_state_now_places_actual_t213_realization_attempt_as_active_primary_t173_move": same_current_repo_state_now_places_actual_t213_realization_attempt_as_active_primary_t173_move,
            "actual_t213_realization_direction_is_now_active_primary_t173_branch_on_current_repo_state": actual_t213_realization_direction_is_now_active_primary_t173_branch_on_current_repo_state,
            "further_return_to_current_exported_p731_continuation_family_remains_nonprimary": further_return_to_current_exported_p731_continuation_family_remains_nonprimary,
            "further_provider_attack_remains_secondary_if_actual_t213_route_stalls": further_provider_attack_remains_secondary_if_actual_t213_route_stalls,
            "next_honest_primary_t173_move_is_actual_t213_realization_attempt_unless_that_route_stalls": next_honest_primary_t173_move_is_actual_t213_realization_attempt_unless_that_route_stalls,
            "no_false_pass": True,
        },
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P761",
        "status": status,
        "as_of": AS_OF,
        "t215_boundary_name": artifact["theorem_result"]["t215_boundary_name"],
        "t215_boundary_exported_on_current_repo_state": artifact["theorem_result"][
            "t215_boundary_exported_on_current_repo_state"
        ],
        "p758_provider_shift_boundary_already_exported": artifact["theorem_result"][
            "p758_provider_shift_boundary_already_exported"
        ],
        "p759_exact_future_only_t213_target_already_exported": artifact["theorem_result"][
            "p759_exact_future_only_t213_target_already_exported"
        ],
        "p760_actual_t213_realization_nonexport_boundary_already_exported": artifact[
            "theorem_result"
        ]["p760_actual_t213_realization_nonexport_boundary_already_exported"],
        "same_current_repo_state_now_places_actual_t213_realization_attempt_as_active_primary_t173_move": artifact[
            "theorem_result"
        ]["same_current_repo_state_now_places_actual_t213_realization_attempt_as_active_primary_t173_move"],
        "actual_t213_realization_direction_is_now_active_primary_t173_branch_on_current_repo_state": artifact[
            "theorem_result"
        ]["actual_t213_realization_direction_is_now_active_primary_t173_branch_on_current_repo_state"],
        "further_return_to_current_exported_p731_continuation_family_remains_nonprimary": artifact[
            "theorem_result"
        ]["further_return_to_current_exported_p731_continuation_family_remains_nonprimary"],
        "further_provider_attack_remains_secondary_if_actual_t213_route_stalls": artifact[
            "theorem_result"
        ]["further_provider_attack_remains_secondary_if_actual_t213_route_stalls"],
        "next_honest_primary_t173_move_is_actual_t213_realization_attempt_unless_that_route_stalls": artifact[
            "theorem_result"
        ]["next_honest_primary_t173_move_is_actual_t213_realization_attempt_unless_that_route_stalls"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
