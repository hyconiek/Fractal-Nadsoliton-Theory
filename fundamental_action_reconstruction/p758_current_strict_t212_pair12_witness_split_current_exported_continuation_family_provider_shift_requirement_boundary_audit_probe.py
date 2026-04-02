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

IN_P709 = GENERATED / "p709_current_strict_release_7_os_residual_sign_gauge_irrelevance_audit_probe_summary.json"
IN_P731 = GENERATED / "p731_current_strict_t185_w_break_witness_payload_residual_datum_pair12_orbit_direction_promotion_bridge_nonexport_audit_probe_summary.json"
IN_P732 = GENERATED / "p732_current_strict_t186_pair1_rooted_convention_state_pair12_witness_split_descent_bridge_nonexport_audit_probe_summary.json"
IN_P733 = GENERATED / "p733_current_strict_t187_convention_layer_pair12_witness_split_transport_descent_bridge_nonexport_audit_probe_summary.json"
IN_P734 = GENERATED / "p734_current_strict_t188_declared_scope_source_topology_selector_theorem_pair12_orbit_direction_descent_bridge_nonexport_audit_probe_summary.json"
IN_P735 = GENERATED / "p735_current_strict_t189_local_source_side_scalar_bind_pair12_orbit_direction_descent_bridge_nonexport_audit_probe_summary.json"
IN_P736 = GENERATED / "p736_current_strict_t190_local_provider_operator_shift_direction_pair12_orbit_direction_descent_bridge_nonexport_audit_probe_summary.json"
IN_P737 = GENERATED / "p737_current_strict_t191_local_pair12_projector_atlas_glue_orbit_direction_descent_bridge_nonexport_audit_probe_summary.json"
IN_P738 = GENERATED / "p738_current_strict_t192_global_projective_selector_state_pair12_orbit_direction_descent_bridge_nonexport_audit_probe_summary.json"
IN_P739 = GENERATED / "p739_current_strict_t193_global_premise_based_directed_selector_state_pair12_witness_split_strict_core_upgrade_bridge_nonexport_audit_probe_summary.json"
IN_P740 = GENERATED / "p740_current_strict_t194_global_sign_fixed_directed_closure_pair12_witness_split_strict_core_upgrade_bridge_nonexport_audit_probe_summary.json"
IN_P741 = GENERATED / "p741_current_strict_t195_actual_source_topology_selector_witness_pair12_witness_split_promotion_bridge_nonexport_audit_probe_summary.json"
IN_P742 = GENERATED / "p742_current_strict_t196_actual_source_topology_selector_witness_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_P743 = GENERATED / "p743_current_strict_t197_source_topology_quotient_safe_qw2191_resolution_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_P744 = GENERATED / "p744_current_strict_t198_declared_scope_source_topology_selector_theorem_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_P745 = GENERATED / "p745_current_strict_t199_declared_scope_source_topology_selector_theorem_target_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_P746 = GENERATED / "p746_current_strict_t200_axiom_augmented_declared_scope_selector_closure_to_residual_datum_pair12_typed_carrier_strict_core_upgrade_bridge_nonexport_audit_probe_summary.json"
IN_P747 = GENERATED / "p747_current_strict_t201_actual_source_topology_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge_nonexport_audit_probe_summary.json"

OUT_JSON = GENERATED / "p758_current_strict_t212_pair12_witness_split_current_exported_continuation_family_provider_shift_requirement_boundary_audit_probe.json"
OUT_SUMMARY = GENERATED / "p758_current_strict_t212_pair12_witness_split_current_exported_continuation_family_provider_shift_requirement_boundary_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [
        IN_P709,
        IN_P731,
        IN_P732,
        IN_P733,
        IN_P734,
        IN_P735,
        IN_P736,
        IN_P737,
        IN_P738,
        IN_P739,
        IN_P740,
        IN_P741,
        IN_P742,
        IN_P743,
        IN_P744,
        IN_P745,
        IN_P746,
        IN_P747,
    ]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P758",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p709 = load_json(IN_P709)
    p731 = load_json(IN_P731)
    p732 = load_json(IN_P732)
    p733 = load_json(IN_P733)
    p734 = load_json(IN_P734)
    p735 = load_json(IN_P735)
    p736 = load_json(IN_P736)
    p737 = load_json(IN_P737)
    p738 = load_json(IN_P738)
    p739 = load_json(IN_P739)
    p740 = load_json(IN_P740)
    p741 = load_json(IN_P741)
    p742 = load_json(IN_P742)
    p743 = load_json(IN_P743)
    p744 = load_json(IN_P744)
    p745 = load_json(IN_P745)
    p746 = load_json(IN_P746)
    p747 = load_json(IN_P747)

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

    p731_pair12_witness_split_remains_live = (
        not bool(p731.get("t185_target_exported_on_current_repo_state"))
        and bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches"))
        and p731.get("pair1_w_break_branch_score_sign") == "negative"
        and p731.get("pair2_w_break_branch_score_sign") == "positive"
        and bool(p731.get("w_break_pair12_branch_scores_are_antisymmetric"))
    )

    current_pair12_witness_split_current_exported_continuation_family_named_members_all_real = all(
        [
            bool(p732.get("current_pair1_rooted_convention_state_exists")),
            bool(p733.get("current_convention_layer_pair12_transport_is_positive_under_all_exported_lifts")),
            bool(p734.get("current_declared_scope_source_topology_selector_theorem_exported")),
            bool(p735.get("current_local_source_side_scalar_witness_family_factors_through_shared_cos_phi_data")),
            bool(p736.get("current_local_provider_operator_shift_direction_lane_realizes_both_pair12_branches_symmetrically")),
            bool(p737.get("current_local_pair12_projector_atlas_glue_lane_exported")),
            bool(p738.get("current_global_projective_selector_state_lane_exported")),
            bool(p739.get("current_global_premise_based_directed_selector_state_lane_exported")),
            bool(p740.get("current_global_sign_fixed_directed_closure_lane_exported")),
            bool(p741.get("current_actual_source_topology_selector_witness_exported")),
            bool(p742.get("current_actual_selector_witness_codomain_has_exported_basis_free_chart_label_forgetting_continuation")),
            bool(p743.get("current_actual_source_topology_quotient_safe_qw2191_resolution_exported")),
            bool(p744.get("current_declared_scope_source_topology_selector_theorem_exported")),
            bool(p745.get("current_declared_scope_source_topology_selector_theorem_target_exported")),
            bool(p746.get("current_actual_nonstrict_declared_scope_selector_closure_exported")),
            bool(p747.get("current_actual_source_topology_selector_witness_target_exported")),
            bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_exported")),
        ]
    )

    current_pair12_witness_split_current_exported_continuation_family_named_members_all_negative = all(
        [
            not bool(p732.get("t186_target_exported_on_current_repo_state"))
            and not bool(p732.get("p731_pair12_witness_split_descends_to_current_pair1_rooted_convention_state")),
            not bool(p733.get("t187_target_exported_on_current_repo_state"))
            and not bool(p733.get("p731_pair12_witness_split_descends_to_current_convention_layer_transport")),
            not bool(p734.get("t188_target_exported_on_current_repo_state"))
            and not bool(p734.get("p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem")),
            not bool(p735.get("t189_target_exported_on_current_repo_state"))
            and not bool(p735.get("p731_pair12_witness_split_descends_to_current_local_source_side_scalar_bind")),
            not bool(p736.get("t190_target_exported_on_current_repo_state"))
            and not bool(p736.get("p731_pair12_witness_split_descends_to_current_local_provider_operator_shift_direction_lane")),
            not bool(p737.get("t191_target_exported_on_current_repo_state"))
            and not bool(p737.get("p731_pair12_witness_split_descends_to_current_local_pair12_projector_atlas_glue_lane")),
            not bool(p738.get("t192_target_exported_on_current_repo_state"))
            and not bool(p738.get("p731_pair12_witness_split_descends_to_current_global_projective_selector_state_lane")),
            not bool(p739.get("t193_target_exported_on_current_repo_state"))
            and not bool(
                p739.get(
                    "p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_premise_based_directed_selector_state_lane"
                )
            ),
            not bool(p740.get("t194_target_exported_on_current_repo_state"))
            and not bool(
                p740.get(
                    "p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_sign_fixed_directed_closure_lane"
                )
            ),
            not bool(p741.get("t195_target_exported_on_current_repo_state"))
            and not bool(p741.get("p731_pair12_witness_split_descends_to_current_actual_source_topology_selector_witness")),
            not bool(p742.get("t196_target_exported_on_current_repo_state"))
            and not bool(
                p742.get(
                    "p731_pair12_witness_split_descends_to_current_actual_selector_witness_to_residual_datum_typed_carrier_bridge"
                )
            ),
            not bool(p743.get("t197_target_exported_on_current_repo_state"))
            and not bool(
                p743.get(
                    "p731_pair12_witness_split_descends_to_current_source_topology_quotient_safe_qw2191_resolution_to_residual_datum_typed_bridge"
                )
            ),
            not bool(p744.get("t198_target_exported_on_current_repo_state"))
            and not bool(
                p744.get(
                    "p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem_to_residual_datum_typed_bridge"
                )
            ),
            not bool(p745.get("t199_target_exported_on_current_repo_state"))
            and not bool(
                p745.get(
                    "p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem_target_to_residual_datum_typed_carrier_bridge"
                )
            ),
            not bool(p746.get("t200_target_exported_on_current_repo_state"))
            and not bool(
                p746.get(
                    "p731_pair12_witness_split_upgrades_to_strict_core_via_current_axiom_augmented_declared_scope_selector_closure"
                )
            ),
            not bool(p747.get("t201_target_exported_on_current_repo_state"))
            and not bool(
                p747.get(
                    "p731_pair12_witness_split_descends_to_current_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge"
                )
            ),
        ]
    )

    release7_os_residual_sign_gauge_irrelevance_already_audited = (
        bool(p709.get("baseline_ok"))
        and bool(p709.get("sign_ok"))
        and float(p709.get("max_abs_diff_p694_m2_under_sign_patterns", 1.0)) == 0.0
        and float(p709.get("max_abs_diff_p696_channel_m2_under_sign_patterns", 1.0)) == 0.0
        and float(p709.get("max_abs_diff_p696_offdiag_to_diag_ratio_under_sign_patterns", 1.0)) == 0.0
        and float(p709.get("max_abs_diff_p696_offblock_max_fro_under_sign_patterns", 1.0)) == 0.0
    )

    same_pair12_witness_split_current_exported_continuation_family_no_longer_admitted_as_active_primary_t173_move = (
        p731_pair12_witness_split_remains_live
        and current_pair12_witness_split_current_exported_continuation_family_named_members_all_real
        and current_pair12_witness_split_current_exported_continuation_family_named_members_all_negative
    )

    provider_shift_is_now_active_primary_t173_branch_on_current_repo_state = (
        same_pair12_witness_split_current_exported_continuation_family_no_longer_admitted_as_active_primary_t173_move
    )

    next_honest_primary_t173_move_requires_genuinely_new_provider_class_beyond_current_exported_continuation_family = (
        same_pair12_witness_split_current_exported_continuation_family_no_longer_admitted_as_active_primary_t173_move
    )

    explicit_release7_os_residual_sign_gauge_freeze_remains_admissible_fallback_if_provider_shift_stalls = (
        same_pair12_witness_split_current_exported_continuation_family_no_longer_admitted_as_active_primary_t173_move
        and release7_os_residual_sign_gauge_irrelevance_already_audited
    )

    add_check(
        "p731_pair12_witness_split_remains_live",
        {
            "t185_target_exported_on_current_repo_state": bool(
                p731.get("t185_target_exported_on_current_repo_state")
            ),
            "current_w_break_witness_payload_separates_pair12_orbit_direction_branches": bool(
                p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches")
            ),
            "pair1_w_break_branch_score_sign": p731.get("pair1_w_break_branch_score_sign"),
            "pair2_w_break_branch_score_sign": p731.get("pair2_w_break_branch_score_sign"),
            "w_break_pair12_branch_scores_are_antisymmetric": bool(
                p731.get("w_break_pair12_branch_scores_are_antisymmetric")
            ),
        },
        {
            "t185_target_exported_on_current_repo_state": False,
            "current_w_break_witness_payload_separates_pair12_orbit_direction_branches": True,
            "pair1_w_break_branch_score_sign": "negative",
            "pair2_w_break_branch_score_sign": "positive",
            "w_break_pair12_branch_scores_are_antisymmetric": True,
        },
        "P731 already freezes one real surviving pair1/pair2 witness split: the current w_break payload separates the opposite residual-datum branches by antisymmetric nonzero signs, but still without one exported promotion bridge.",
    )
    add_check(
        "current_pair12_witness_split_current_exported_continuation_family_named_members_all_real",
        current_pair12_witness_split_current_exported_continuation_family_named_members_all_real,
        True,
        "The currently named exported continuation family after P731 is real enough to count as tested: convention, theorem, local, global, source-topology, theorem-target, nonstrict-upgrade, and local-atlas-side members are already exported on current repo state.",
    )
    add_check(
        "current_pair12_witness_split_current_exported_continuation_family_named_members_all_negative",
        current_pair12_witness_split_current_exported_continuation_family_named_members_all_negative,
        True,
        "Each currently named exported continuation-family member audited in P732-P747 remains negative on the exact question of descending or upgrading the opposite P731 pair1/pair2 witness split into one typed strict branch distinction.",
    )
    add_check(
        "release7_os_residual_sign_gauge_irrelevance_already_audited",
        release7_os_residual_sign_gauge_irrelevance_already_audited,
        True,
        "P709/N706 already audit that residual sign remains gauge-irrelevant for the concrete Release-7 OS observables, so explicit gauge freeze remains an honest fallback for those observables if provider shift later stalls.",
    )
    add_check(
        "same_pair12_witness_split_current_exported_continuation_family_no_longer_admitted_as_active_primary_t173_move",
        same_pair12_witness_split_current_exported_continuation_family_no_longer_admitted_as_active_primary_t173_move,
        True,
        "Therefore the same currently named exported continuation family after P731 may no longer be treated as the active primary T173 move on the current repo state.",
    )
    add_check(
        "provider_shift_is_now_active_primary_t173_branch_on_current_repo_state",
        provider_shift_is_now_active_primary_t173_branch_on_current_repo_state,
        True,
        "So a genuinely new provider-class shift beyond the current exported P731 continuation family is now the active primary T173 branch on current repo state.",
    )
    add_check(
        "next_honest_primary_t173_move_requires_genuinely_new_provider_class_beyond_current_exported_continuation_family",
        next_honest_primary_t173_move_requires_genuinely_new_provider_class_beyond_current_exported_continuation_family,
        True,
        "Hence the next honest primary strict move under T173 now requires a genuinely new provider class beyond the current exported continuation family rather than another continuation inside it.",
    )
    add_check(
        "explicit_release7_os_residual_sign_gauge_freeze_remains_admissible_fallback_if_provider_shift_stalls",
        explicit_release7_os_residual_sign_gauge_freeze_remains_admissible_fallback_if_provider_shift_stalls,
        True,
        "If that provider shift later stalls too, explicit residual-sign gauge freeze remains an admissible fallback for the already audited Release-7 OS observables because their outputs are already sign-gauge-invariant.",
    )

    status = (
        "PASS_STRICT_T212_PAIR12_WITNESS_SPLIT_CURRENT_EXPORTED_CONTINUATION_FAMILY_PROVIDER_SHIFT_REQUIREMENT_BOUNDARY_AUDITED"
        if not blocking
        and same_pair12_witness_split_current_exported_continuation_family_no_longer_admitted_as_active_primary_t173_move
        and provider_shift_is_now_active_primary_t173_branch_on_current_repo_state
        and next_honest_primary_t173_move_requires_genuinely_new_provider_class_beyond_current_exported_continuation_family
        else "FAIL_STRICT_T212_PAIR12_WITNESS_SPLIT_CURRENT_EXPORTED_CONTINUATION_FAMILY_PROVIDER_SHIFT_REQUIREMENT_BOUNDARY_AUDIT"
    )

    artifact = {
        "stage": "P758",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "theorem_result": {
            "t212_boundary_name": "Pair12WitnessSplitCurrentExportedContinuationFamilyProviderShiftRequirementBoundary_strict_v1",
            "t212_boundary_exported_on_current_repo_state": status
            == "PASS_STRICT_T212_PAIR12_WITNESS_SPLIT_CURRENT_EXPORTED_CONTINUATION_FAMILY_PROVIDER_SHIFT_REQUIREMENT_BOUNDARY_AUDITED",
            "p731_pair12_witness_split_remains_live": p731_pair12_witness_split_remains_live,
            "current_pair12_witness_split_current_exported_continuation_family_named_members_all_real": current_pair12_witness_split_current_exported_continuation_family_named_members_all_real,
            "current_pair12_witness_split_current_exported_continuation_family_named_members_all_negative": current_pair12_witness_split_current_exported_continuation_family_named_members_all_negative,
            "release7_os_residual_sign_gauge_irrelevance_already_audited": release7_os_residual_sign_gauge_irrelevance_already_audited,
            "same_pair12_witness_split_current_exported_continuation_family_no_longer_admitted_as_active_primary_t173_move": same_pair12_witness_split_current_exported_continuation_family_no_longer_admitted_as_active_primary_t173_move,
            "provider_shift_is_now_active_primary_t173_branch_on_current_repo_state": provider_shift_is_now_active_primary_t173_branch_on_current_repo_state,
            "next_honest_primary_t173_move_requires_genuinely_new_provider_class_beyond_current_exported_continuation_family": next_honest_primary_t173_move_requires_genuinely_new_provider_class_beyond_current_exported_continuation_family,
            "explicit_release7_os_residual_sign_gauge_freeze_remains_admissible_fallback_if_provider_shift_stalls": explicit_release7_os_residual_sign_gauge_freeze_remains_admissible_fallback_if_provider_shift_stalls,
            "no_false_pass": True,
        },
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P758",
        "status": status,
        "as_of": AS_OF,
        "t212_boundary_name": artifact["theorem_result"]["t212_boundary_name"],
        "t212_boundary_exported_on_current_repo_state": artifact["theorem_result"][
            "t212_boundary_exported_on_current_repo_state"
        ],
        "p731_pair12_witness_split_remains_live": artifact["theorem_result"][
            "p731_pair12_witness_split_remains_live"
        ],
        "current_pair12_witness_split_current_exported_continuation_family_named_members_all_real": artifact[
            "theorem_result"
        ]["current_pair12_witness_split_current_exported_continuation_family_named_members_all_real"],
        "current_pair12_witness_split_current_exported_continuation_family_named_members_all_negative": artifact[
            "theorem_result"
        ]["current_pair12_witness_split_current_exported_continuation_family_named_members_all_negative"],
        "release7_os_residual_sign_gauge_irrelevance_already_audited": artifact["theorem_result"][
            "release7_os_residual_sign_gauge_irrelevance_already_audited"
        ],
        "same_pair12_witness_split_current_exported_continuation_family_no_longer_admitted_as_active_primary_t173_move": artifact[
            "theorem_result"
        ]["same_pair12_witness_split_current_exported_continuation_family_no_longer_admitted_as_active_primary_t173_move"],
        "provider_shift_is_now_active_primary_t173_branch_on_current_repo_state": artifact[
            "theorem_result"
        ]["provider_shift_is_now_active_primary_t173_branch_on_current_repo_state"],
        "next_honest_primary_t173_move_requires_genuinely_new_provider_class_beyond_current_exported_continuation_family": artifact[
            "theorem_result"
        ]["next_honest_primary_t173_move_requires_genuinely_new_provider_class_beyond_current_exported_continuation_family"],
        "explicit_release7_os_residual_sign_gauge_freeze_remains_admissible_fallback_if_provider_shift_stalls": artifact[
            "theorem_result"
        ]["explicit_release7_os_residual_sign_gauge_freeze_remains_admissible_fallback_if_provider_shift_stalls"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
