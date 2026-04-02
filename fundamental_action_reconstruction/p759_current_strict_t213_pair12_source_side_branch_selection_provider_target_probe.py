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

IN_P729 = GENERATED / "p729_current_strict_t183_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe_summary.json"
IN_P731 = GENERATED / "p731_current_strict_t185_w_break_witness_payload_residual_datum_pair12_orbit_direction_promotion_bridge_nonexport_audit_probe_summary.json"
IN_P758 = GENERATED / "p758_current_strict_t212_pair12_witness_split_current_exported_continuation_family_provider_shift_requirement_boundary_audit_probe_summary.json"
IN_T173 = ROOT / "T173_CURRENT_STRICT_CORE_SELECTOR_CLOSURE_AND_KERNEL_ALONE_QW2191_DISCHARGE_TARGET_SPEC.md"
IN_T213 = ROOT / "T213_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p759_current_strict_t213_pair12_source_side_branch_selection_provider_target_probe.json"
OUT_SUMMARY = GENERATED / "p759_current_strict_t213_pair12_source_side_branch_selection_provider_target_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P729, IN_P731, IN_P758, IN_T173, IN_T213]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P759",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p729 = load_json(IN_P729)
    p731 = load_json(IN_P731)
    p758 = load_json(IN_P758)
    t173_text = load_text(IN_T173)
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

    p729_pair12_split_localized_as_opposite_orbit_directions = (
        not bool(p729.get("t183_target_exported_on_current_repo_state"))
        and bool(p729.get("remaining_pair12_split_localized_as_opposite_orbit_directions"))
        and p729.get("pair1_orbit_branch_kind") == "delta_k_positive_index_branch"
        and p729.get("pair2_orbit_branch_kind") == "delta_minus_k_negative_index_branch"
    )

    p731_w_break_witness_split_already_separates_pair12_branches = (
        not bool(p731.get("t185_target_exported_on_current_repo_state"))
        and bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches"))
        and p731.get("pair1_w_break_branch_score_sign") == "negative"
        and p731.get("pair2_w_break_branch_score_sign") == "positive"
        and bool(p731.get("w_break_pair12_branch_scores_are_antisymmetric"))
    )

    p758_provider_shift_boundary_already_exports_need_for_genuinely_new_provider_class = (
        bool(p758.get("t212_boundary_exported_on_current_repo_state"))
        and bool(p758.get("same_pair12_witness_split_current_exported_continuation_family_no_longer_admitted_as_active_primary_t173_move"))
        and bool(p758.get("provider_shift_is_now_active_primary_t173_branch_on_current_repo_state"))
        and bool(
            p758.get(
                "next_honest_primary_t173_move_requires_genuinely_new_provider_class_beyond_current_exported_continuation_family"
            )
        )
    )

    t173_discipline_already_frozen = all(
        needle in t173_text
        for needle in [
            "QW2191_kernel_alone_discharge = false",
            "physical orientation datum in strict core remain explicitly unclaimed",
            "no silent import of `QW-2192` into strict core",
        ]
    )

    t213_needles = {
        "target_symbol": "W_strict_t173_pair12_source_side_branch_selection_provider_target_v1"
        in t213_text,
        "source_side": "target_is_source_side := yes" in t213_text,
        "observer_free": "target_is_observer_free := yes" in t213_text,
        "pair12_typed": "target_is_pair12_typed := yes" in t213_text,
        "residual_datum_carrier_binding": "target_is_residual_datum_carrier_binding := yes"
        in t213_text,
        "branch_sensitive": "target_is_branch_sensitive_to_delta_k_vs_delta_minus_k := yes"
        in t213_text,
        "chart_sensitive": "target_is_chart_sensitive_or_chart_label_retaining := yes"
        in t213_text,
        "external_to_current_exported_family": "target_is_external_to_current_exported_p731_continuation_family := yes"
        in t213_text,
        "nonconvention_nonpremise_based": "target_is_nonconvention_nonpremise_based := yes"
        in t213_text,
        "below_t176": "target_remains_below_global_t176_discharge := yes" in t213_text,
        "below_physical_sign_claim": "target_remains_below_directed_physical_orientation_datum_claim := yes"
        in t213_text,
        "future_route_only": "future_route_only := yes" in t213_text,
    }

    current_t213_target_is_source_side_observer_free = (
        t213_needles["source_side"] and t213_needles["observer_free"]
    )
    current_t213_target_is_pair12_typed_and_branch_sensitive = (
        t213_needles["pair12_typed"] and t213_needles["branch_sensitive"]
    )
    current_t213_target_is_chart_sensitive_and_residual_datum_carrier_binding = (
        t213_needles["chart_sensitive"] and t213_needles["residual_datum_carrier_binding"]
    )
    current_t213_target_is_external_to_current_exported_p731_continuation_family = t213_needles[
        "external_to_current_exported_family"
    ]
    current_t213_target_is_nonconvention_nonpremise_based = t213_needles[
        "nonconvention_nonpremise_based"
    ]
    current_t213_target_remains_below_t176_discharge_and_physical_sign_datum_claim = (
        t213_needles["below_t176"] and t213_needles["below_physical_sign_claim"]
    )
    current_t213_target_is_future_route_only = t213_needles["future_route_only"]

    add_check(
        "p729_pair12_split_localized_as_opposite_orbit_directions",
        {
            "t183_target_exported_on_current_repo_state": bool(
                p729.get("t183_target_exported_on_current_repo_state")
            ),
            "remaining_pair12_split_localized_as_opposite_orbit_directions": bool(
                p729.get("remaining_pair12_split_localized_as_opposite_orbit_directions")
            ),
            "pair1_orbit_branch_kind": p729.get("pair1_orbit_branch_kind"),
            "pair2_orbit_branch_kind": p729.get("pair2_orbit_branch_kind"),
        },
        {
            "t183_target_exported_on_current_repo_state": False,
            "remaining_pair12_split_localized_as_opposite_orbit_directions": True,
            "pair1_orbit_branch_kind": "delta_k_positive_index_branch",
            "pair2_orbit_branch_kind": "delta_minus_k_negative_index_branch",
        },
        "P729 already freezes that the surviving pair1/pair2 ambiguity is exactly the opposite orbit-direction split delta_k versus delta_-k.",
    )
    add_check(
        "p731_w_break_witness_split_already_separates_pair12_branches",
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
        "P731 already freezes that the surviving pair1/pair2 split is not invisible anymore: w_break separates the two branches by opposite nonzero signs.",
    )
    add_check(
        "p758_provider_shift_boundary_already_exports_need_for_genuinely_new_provider_class",
        {
            "t212_boundary_exported_on_current_repo_state": bool(
                p758.get("t212_boundary_exported_on_current_repo_state")
            ),
            "same_pair12_witness_split_current_exported_continuation_family_no_longer_admitted_as_active_primary_t173_move": bool(
                p758.get(
                    "same_pair12_witness_split_current_exported_continuation_family_no_longer_admitted_as_active_primary_t173_move"
                )
            ),
            "provider_shift_is_now_active_primary_t173_branch_on_current_repo_state": bool(
                p758.get("provider_shift_is_now_active_primary_t173_branch_on_current_repo_state")
            ),
            "next_honest_primary_t173_move_requires_genuinely_new_provider_class_beyond_current_exported_continuation_family": bool(
                p758.get(
                    "next_honest_primary_t173_move_requires_genuinely_new_provider_class_beyond_current_exported_continuation_family"
                )
            ),
        },
        {
            "t212_boundary_exported_on_current_repo_state": True,
            "same_pair12_witness_split_current_exported_continuation_family_no_longer_admitted_as_active_primary_t173_move": True,
            "provider_shift_is_now_active_primary_t173_branch_on_current_repo_state": True,
            "next_honest_primary_t173_move_requires_genuinely_new_provider_class_beyond_current_exported_continuation_family": True,
        },
        "P758 already freezes that the next honest T173 move must now leave the current exported continuation family and attack a genuinely new provider class.",
    )
    add_check(
        "t173_discipline_already_frozen",
        t173_discipline_already_frozen,
        True,
        "T173 already freezes no-false-pass discipline: no silent kernel-alone QW-2191 discharge and no silent physical sign datum claim in strict core.",
    )
    add_check(
        "t213_target_spec_exported_with_required_properties",
        t213_needles,
        {
            "target_symbol": True,
            "source_side": True,
            "observer_free": True,
            "pair12_typed": True,
            "residual_datum_carrier_binding": True,
            "branch_sensitive": True,
            "chart_sensitive": True,
            "external_to_current_exported_family": True,
            "nonconvention_nonpremise_based": True,
            "below_t176": True,
            "below_physical_sign_claim": True,
            "future_route_only": True,
        },
        "T213 exports one exact future-only target for the genuinely new provider class now required beyond the current exported P731 continuation family.",
    )
    add_check(
        "current_t213_target_is_source_side_observer_free",
        current_t213_target_is_source_side_observer_free,
        True,
        "The T213 target remains source-side and observer-free.",
    )
    add_check(
        "current_t213_target_is_pair12_typed_and_branch_sensitive",
        current_t213_target_is_pair12_typed_and_branch_sensitive,
        True,
        "The T213 target is genuinely typed on pair1/pair2 and branch-sensitive to delta_k versus delta_-k.",
    )
    add_check(
        "current_t213_target_is_chart_sensitive_and_residual_datum_carrier_binding",
        current_t213_target_is_chart_sensitive_and_residual_datum_carrier_binding,
        True,
        "The T213 target remains chart-sensitive / chart-label-retaining and binds the surviving residual-datum carrier frontier directly.",
    )
    add_check(
        "current_t213_target_is_external_to_current_exported_p731_continuation_family",
        current_t213_target_is_external_to_current_exported_p731_continuation_family,
        True,
        "The T213 target is explicitly external to the currently exported P731 -> P747 continuation family.",
    )
    add_check(
        "current_t213_target_is_nonconvention_nonpremise_based",
        current_t213_target_is_nonconvention_nonpremise_based,
        True,
        "The T213 target is explicitly not just another convention-layer or premise-based reformulation.",
    )
    add_check(
        "current_t213_target_remains_below_t176_discharge_and_physical_sign_datum_claim",
        current_t213_target_remains_below_t176_discharge_and_physical_sign_datum_claim,
        True,
        "The T213 target remains below global T176 discharge and below any directed physical sign-datum claim in strict core.",
    )
    add_check(
        "current_t213_target_is_future_route_only",
        current_t213_target_is_future_route_only,
        True,
        "The T213 target remains explicitly future-route-only.",
    )

    status = (
        "PASS_STRICT_T213_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_TARGET_EXPORTED"
        if not blocking
        and p729_pair12_split_localized_as_opposite_orbit_directions
        and p731_w_break_witness_split_already_separates_pair12_branches
        and p758_provider_shift_boundary_already_exports_need_for_genuinely_new_provider_class
        and t173_discipline_already_frozen
        and current_t213_target_is_source_side_observer_free
        and current_t213_target_is_pair12_typed_and_branch_sensitive
        and current_t213_target_is_chart_sensitive_and_residual_datum_carrier_binding
        and current_t213_target_is_external_to_current_exported_p731_continuation_family
        and current_t213_target_is_nonconvention_nonpremise_based
        and current_t213_target_remains_below_t176_discharge_and_physical_sign_datum_claim
        and current_t213_target_is_future_route_only
        else "FAIL_STRICT_T213_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_TARGET_EXPORT"
    )

    artifact = {
        "stage": "P759",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "theorem_result": {
            "t213_target_name": "W_strict_t173_pair12_source_side_branch_selection_provider_target_v1",
            "t213_target_exported_on_current_repo_state": status
            == "PASS_STRICT_T213_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_TARGET_EXPORTED",
            "p729_pair12_split_localized_as_opposite_orbit_directions": p729_pair12_split_localized_as_opposite_orbit_directions,
            "p731_w_break_witness_split_already_separates_pair12_branches": p731_w_break_witness_split_already_separates_pair12_branches,
            "p758_provider_shift_boundary_already_exports_need_for_genuinely_new_provider_class": p758_provider_shift_boundary_already_exports_need_for_genuinely_new_provider_class,
            "t173_discipline_already_frozen": t173_discipline_already_frozen,
            "current_t213_target_is_source_side_observer_free": current_t213_target_is_source_side_observer_free,
            "current_t213_target_is_pair12_typed_and_branch_sensitive": current_t213_target_is_pair12_typed_and_branch_sensitive,
            "current_t213_target_is_chart_sensitive_and_residual_datum_carrier_binding": current_t213_target_is_chart_sensitive_and_residual_datum_carrier_binding,
            "current_t213_target_is_external_to_current_exported_p731_continuation_family": current_t213_target_is_external_to_current_exported_p731_continuation_family,
            "current_t213_target_is_nonconvention_nonpremise_based": current_t213_target_is_nonconvention_nonpremise_based,
            "current_t213_target_remains_below_t176_discharge_and_physical_sign_datum_claim": current_t213_target_remains_below_t176_discharge_and_physical_sign_datum_claim,
            "current_t213_target_is_future_route_only": current_t213_target_is_future_route_only,
            "no_false_pass": True,
        },
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P759",
        "status": status,
        "as_of": AS_OF,
        "t213_target_name": artifact["theorem_result"]["t213_target_name"],
        "t213_target_exported_on_current_repo_state": artifact["theorem_result"][
            "t213_target_exported_on_current_repo_state"
        ],
        "p729_pair12_split_localized_as_opposite_orbit_directions": artifact["theorem_result"][
            "p729_pair12_split_localized_as_opposite_orbit_directions"
        ],
        "p731_w_break_witness_split_already_separates_pair12_branches": artifact["theorem_result"][
            "p731_w_break_witness_split_already_separates_pair12_branches"
        ],
        "p758_provider_shift_boundary_already_exports_need_for_genuinely_new_provider_class": artifact[
            "theorem_result"
        ]["p758_provider_shift_boundary_already_exports_need_for_genuinely_new_provider_class"],
        "t173_discipline_already_frozen": artifact["theorem_result"][
            "t173_discipline_already_frozen"
        ],
        "current_t213_target_is_source_side_observer_free": artifact["theorem_result"][
            "current_t213_target_is_source_side_observer_free"
        ],
        "current_t213_target_is_pair12_typed_and_branch_sensitive": artifact["theorem_result"][
            "current_t213_target_is_pair12_typed_and_branch_sensitive"
        ],
        "current_t213_target_is_chart_sensitive_and_residual_datum_carrier_binding": artifact[
            "theorem_result"
        ]["current_t213_target_is_chart_sensitive_and_residual_datum_carrier_binding"],
        "current_t213_target_is_external_to_current_exported_p731_continuation_family": artifact[
            "theorem_result"
        ]["current_t213_target_is_external_to_current_exported_p731_continuation_family"],
        "current_t213_target_is_nonconvention_nonpremise_based": artifact["theorem_result"][
            "current_t213_target_is_nonconvention_nonpremise_based"
        ],
        "current_t213_target_remains_below_t176_discharge_and_physical_sign_datum_claim": artifact[
            "theorem_result"
        ]["current_t213_target_remains_below_t176_discharge_and_physical_sign_datum_claim"],
        "current_t213_target_is_future_route_only": artifact["theorem_result"][
            "current_t213_target_is_future_route_only"
        ],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
