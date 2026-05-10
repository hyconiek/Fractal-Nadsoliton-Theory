#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P708 = GENERATED / "p708_current_strict_t173_frontier_dashboard_probe_summary.json"
IN_P978 = GENERATED / "p978_current_strict_t269_t268_attempt_verdict_or_post_even_deeper_boundary_nonexport_probe_summary.json"
IN_P980 = GENERATED / "p980_current_strict_t271_t270_post_even_deeper_boundary_actual_nonexport_probe_summary.json"
IN_P981 = GENERATED / "p981_current_strict_t272_t270_post_even_deeper_boundary_actual_attempt_probe_summary.json"
IN_F18 = ROOT / "sandbox_strict_core_ingredient_attempt/F18_SANDBOX_STRICT_CORE_THETA_SUPPLY_POPULATION_LOOP_INCOMPATIBILITY_BOUNDARY_PACKET.md"
IN_T18 = ROOT / "sandbox_strict_core_ingredient_attempt/T18_SANDBOX_STRICT_CORE_THETA_SUPPLY_POPULATION_LOOP_INCOMPATIBILITY_BOUNDARY_SCOPE.md"
IN_P881 = ROOT / "P881_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_SAME_LANE_EXHAUSTION_BOUNDARY_AUDIT_PROBE.md"
IN_N355 = ROOT / "N355_CURRENT_FIRST_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_NONCYCLIC_PROVIDER_SPLIT_TARGET_THEOREM.md"
IN_N361 = ROOT / "N361_CURRENT_FIRST_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_PAIR_REALIZATION_SIDE_PROVIDER_SUPPORT_PACKET_THEOREM.md"
IN_N362 = ROOT / "N362_CURRENT_FIRST_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_PAIR_REALIZATION_SIDE_PROVIDER_SUPPORT_WITNESS_THEOREM.md"
IN_N463 = ROOT / "N463_CURRENT_FIRST_STRICT_QW2191_SHANNON_SITE_AMPLITUDE_ENTROPY_O2_CUT_NONUNIQUENESS_BOUNDARY_THEOREM.md"
IN_N464 = ROOT / "N464_CURRENT_FIRST_STRICT_QW2191_SITE_DISTRIBUTION_PERMUTATION_INVARIANT_OBJECTIVE_O2_CUT_NONUNIQUENESS_BOUNDARY_THEOREM.md"

OUT_JSON = GENERATED / "p982_current_strict_qw2191_pair12_entry_point_same_lane_exhaustion_and_noncyclic_pivot_audit_probe.json"
OUT_SUMMARY = GENERATED / "p982_current_strict_qw2191_pair12_entry_point_same_lane_exhaustion_and_noncyclic_pivot_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [
        IN_P708,
        IN_P978,
        IN_P980,
        IN_P981,
        IN_F18,
        IN_T18,
        IN_P881,
        IN_N355,
        IN_N361,
        IN_N362,
        IN_N463,
        IN_N464,
    ]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P982",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p708 = load_json(IN_P708)
    p978 = load_json(IN_P978)
    p980 = load_json(IN_P980)
    p981 = load_json(IN_P981)
    f18_text = load_text(IN_F18)
    t18_text = load_text(IN_T18)
    p881_text = load_text(IN_P881)
    n355_text = load_text(IN_N355)
    n361_text = load_text(IN_N361)
    n362_text = load_text(IN_N362)
    n463_text = load_text(IN_N463)
    n464_text = load_text(IN_N464)

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

    current_lane_repeats_same_target_nonexport_attempt_pattern_under_same_blocker_cut = (
        p978.get("status")
        == "PASS_STRICT_T269_PAIR12_ENTRY_POINT_EXACT_EVEN_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_POST_EVEN_DEEPER_BOUNDARY_REFINEMENT_NONEXPORT_AUDITED"
        and p980.get("status")
        == "PASS_STRICT_T271_PAIR12_ENTRY_POINT_EXACT_POST_EVEN_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and p981.get("status")
        == "PASS_STRICT_T272_PAIR12_ENTRY_POINT_EXACT_POST_EVEN_DEEPER_BOUNDARY_REFINEMENT_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        and bool(p981.get("first_actual_t270_exact_post_even_deeper_boundary_refinement_realization_attempt_keeps_verdict_and_next_object_open"))
    )

    current_lane_still_has_no_global_provider_upgrade_and_no_new_selector_source = (
        not bool(p708.get("t176_global_provider_exported"))
        and not bool(p708.get("directed_sign_sensitive_physical_orientation_in_strict_core"))
        and not bool(p708.get("QW2191_kernel_alone_discharge"))
    )

    repo_noncyclic_guardrail_already_rejects_same_blocker_cut_positive_repetition = (
        "same_blocker_cut_recurs_if_positive_lift_repeats := yes" in f18_text
        and "noncyclic_anchor_present_inside_current_sandbox := no" in f18_text
        and "further same-loop positive lifting is not the honest next move" in f18_text
        and "further positive lifting inside the same loop would be cyclic repetition" in t18_text
    )

    repo_already_admits_same_lane_exhaustion_boundary_auditing = (
        "same-lane passive groundwork is exhausted" in p881_text
        and "noncyclic continuation discipline forbids fabricated repeated splitting" in p881_text
    )

    strict_qw2191_naive_shannon_uniqueness_route_already_closed = (
        "cannot yield a unique O(2)-cut" in n463_text
        and "cannot yield a unique O(2)-cut" in n464_text
        and "discharge of `QW-2191`" in n463_text
        and "discharge of `QW-2191`" in n464_text
    )

    exported_noncyclic_provider_split_family_already_exists_for_pivot = (
        "Xi_nad12_sigma_residual_shannon_noncyclic_provider_split_target_v1" in n355_text
        and "Lambda_nad12_sigma_residual_shannon_pair_realization_side_provider_support_witness_v1" in n362_text
        and "Upsilon_nad12_sigma_residual_shannon_pair_realization_side_provider_support_packet_v1" in n361_text
        and "the next honest move on this route, if continued, must therefore be:" in n355_text
    )

    same_lane_exhaustion_boundary_reached = (
        current_lane_repeats_same_target_nonexport_attempt_pattern_under_same_blocker_cut
        and current_lane_still_has_no_global_provider_upgrade_and_no_new_selector_source
        and repo_noncyclic_guardrail_already_rejects_same_blocker_cut_positive_repetition
        and repo_already_admits_same_lane_exhaustion_boundary_auditing
        and strict_qw2191_naive_shannon_uniqueness_route_already_closed
        and exported_noncyclic_provider_split_family_already_exists_for_pivot
    )

    further_same_lane_t274_style_descent_is_not_honest_primary_move = (
        same_lane_exhaustion_boundary_reached
    )

    next_honest_move_is_noncyclic_pivot_to_exported_provider_split_family = (
        same_lane_exhaustion_boundary_reached
    )

    preferred_noncyclic_pivot_family = (
        "nad12_sigma_residual_shannon_noncyclic_provider_split_family"
        if next_honest_move_is_noncyclic_pivot_to_exported_provider_split_family
        else None
    )

    preferred_first_pivot_branch = (
        "pair_realization_side_provider_support_witness_family"
        if next_honest_move_is_noncyclic_pivot_to_exported_provider_split_family
        else None
    )

    add_check(
        "current_lane_repeats_same_target_nonexport_attempt_pattern_under_same_blocker_cut",
        current_lane_repeats_same_target_nonexport_attempt_pattern_under_same_blocker_cut,
        True,
        "The current QW-2191 local lane repeats the same target/nonexport/attempt grammar under one blocker-cut.",
    )
    add_check(
        "current_lane_still_has_no_global_provider_upgrade_and_no_new_selector_source",
        current_lane_still_has_no_global_provider_upgrade_and_no_new_selector_source,
        True,
        "The lane still has no T176 upgrade and no new strict internal selector source.",
    )
    add_check(
        "repo_noncyclic_guardrail_already_rejects_same_blocker_cut_positive_repetition",
        repo_noncyclic_guardrail_already_rejects_same_blocker_cut_positive_repetition,
        True,
        "Repo noncyclic guardrails already reject repeated positive lifting under the same blocker-cut as a primary strategy.",
    )
    add_check(
        "repo_already_admits_same_lane_exhaustion_boundary_auditing",
        repo_already_admits_same_lane_exhaustion_boundary_auditing,
        True,
        "Repo already admits same-lane exhaustion auditing as an honest boundary idiom.",
    )
    add_check(
        "strict_qw2191_naive_shannon_uniqueness_route_already_closed",
        strict_qw2191_naive_shannon_uniqueness_route_already_closed,
        True,
        "Naive permutation-invariant Shannon objective families are already closed as strict QW-2191 uniqueness routes.",
    )
    add_check(
        "exported_noncyclic_provider_split_family_already_exists_for_pivot",
        exported_noncyclic_provider_split_family_already_exists_for_pivot,
        True,
        "The repo already exports a noncyclic provider-split family available for pivot.",
    )
    add_check(
        "same_lane_exhaustion_boundary_reached",
        same_lane_exhaustion_boundary_reached,
        True,
        "Therefore the current same-lane QW-2191 pair12 entry-point descent has reached its honest same-lane exhaustion boundary.",
    )
    add_check(
        "further_same_lane_t274_style_descent_is_not_honest_primary_move",
        further_same_lane_t274_style_descent_is_not_honest_primary_move,
        True,
        "One more T274-style same-lane descent is no longer the honest primary move.",
    )
    add_check(
        "next_honest_move_is_noncyclic_pivot_to_exported_provider_split_family",
        next_honest_move_is_noncyclic_pivot_to_exported_provider_split_family,
        True,
        "The next honest move is now a noncyclic pivot to the exported provider-split family.",
    )

    status = (
        "PASS_STRICT_QW2191_PAIR12_ENTRY_POINT_SAME_LANE_EXHAUSTION_AND_NONCYCLIC_PIVOT_AUDITED"
        if not blocking and same_lane_exhaustion_boundary_reached
        else "FAIL_STRICT_QW2191_PAIR12_ENTRY_POINT_SAME_LANE_EXHAUSTION_AND_NONCYCLIC_PIVOT_AUDIT"
    )

    artifact = {
        "stage": "P982",
        "status": status,
        "as_of": AS_OF,
        "lane": "qw2191_pair12_entry_point_same_lane_exhaustion_only",
        "same_lane_exhaustion_boundary_reached": same_lane_exhaustion_boundary_reached,
        "further_same_lane_t274_style_descent_is_not_honest_primary_move": further_same_lane_t274_style_descent_is_not_honest_primary_move,
        "next_honest_move_is_noncyclic_pivot_to_exported_provider_split_family": next_honest_move_is_noncyclic_pivot_to_exported_provider_split_family,
        "preferred_noncyclic_pivot_family": preferred_noncyclic_pivot_family,
        "preferred_first_pivot_branch": preferred_first_pivot_branch,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "lane": artifact["lane"],
        "same_lane_exhaustion_boundary_reached": artifact["same_lane_exhaustion_boundary_reached"],
        "further_same_lane_t274_style_descent_is_not_honest_primary_move": artifact["further_same_lane_t274_style_descent_is_not_honest_primary_move"],
        "next_honest_move_is_noncyclic_pivot_to_exported_provider_split_family": artifact["next_honest_move_is_noncyclic_pivot_to_exported_provider_split_family"],
        "preferred_noncyclic_pivot_family": artifact["preferred_noncyclic_pivot_family"],
        "preferred_first_pivot_branch": artifact["preferred_first_pivot_branch"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
