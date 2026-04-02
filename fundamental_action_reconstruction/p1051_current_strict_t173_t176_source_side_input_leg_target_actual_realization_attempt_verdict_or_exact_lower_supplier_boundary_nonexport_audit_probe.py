#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1050 = GENERATED / "p1050_current_strict_t173_t176_source_side_input_leg_target_actual_realization_attempt_probe_summary.json"
IN_F949 = GENERATED / "f949_first_current_strict_qw2191_pair12_entry_point_same_lane_exhaustion_and_noncyclic_pivot_packet_summary.json"
IN_P766 = GENERATED / "p766_current_strict_t220_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_probe_summary.json"
IN_P767 = GENERATED / "p767_current_strict_t221_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_nonexport_audit_probe_summary.json"
IN_P768 = GENERATED / "p768_current_strict_t222_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_target_probe_summary.json"
IN_P769 = GENERATED / "p769_current_strict_t223_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_nonexport_audit_probe_summary.json"
IN_P770 = GENERATED / "p770_current_strict_t224_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_attempt_probe_summary.json"
IN_T302 = ROOT / "T302_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p1051_current_strict_t173_t176_source_side_input_leg_target_actual_realization_attempt_verdict_or_exact_lower_supplier_boundary_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1051_current_strict_t173_t176_source_side_input_leg_target_actual_realization_attempt_verdict_or_exact_lower_supplier_boundary_nonexport_audit_probe_summary.json"

P1050_STATUS = "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
F949_STATUS = (
    "PASS_CURRENT_STRICT_QW2191_PAIR12_ENTRY_POINT_SAME_LANE_EXHAUSTION_AND_NONCYCLIC_PIVOT_PACKET_EXPORTED"
)
P766_STATUS = (
    "PASS_STRICT_T220_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
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

T302_ATTEMPT_NAME = (
    "W_strict_t173_t176_inversion_sensitive_pair12_branch_separation_bridge_source_side_input_leg_actual_realization_attempt_v1"
)
T220_ATTEMPT_NAME = (
    "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_v1"
)
T224_ATTEMPT_NAME = (
    "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_attempt_v1"
)
EXACT_NAMED_SUBINTERFACE = (
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


def scan_t302_followup_hits() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py")
    excluded = {
        Path(__file__).name,
        "P1051_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_LOWER_SUPPLIER_BOUNDARY_NONEXPORT_AUDIT_PROBE.md",
        "N884_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_LOWER_SUPPLIER_BOUNDARY_NONEXPORT_AUDIT_THEOREM.md",
        "T302_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p1050_current_strict_t173_t176_source_side_input_leg_target_actual_realization_attempt_probe.py",
    }
    hits: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if path in seen or path.name in excluded:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            if T302_ATTEMPT_NAME not in text:
                continue
            if any(
                marker in text
                for marker in (
                    "verdict",
                    "lower_supplier_boundary",
                    "lower supplier-boundary",
                    "lower supplier boundary",
                )
            ):
                hits.append(rel(path))
    return hits


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1050, IN_F949, IN_P766, IN_P767, IN_P768, IN_P769, IN_P770, IN_T302]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1051",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1050 = load_json(IN_P1050)
    f949 = load_json(IN_F949)
    p766 = load_json(IN_P766)
    p767 = load_json(IN_P767)
    p768 = load_json(IN_P768)
    p769 = load_json(IN_P769)
    p770 = load_json(IN_P770)
    t302_text = load_text(IN_T302)
    followup_hits = scan_t302_followup_hits()

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

    t302_exact_attempt_is_frozen = (
        p1050.get("status") == P1050_STATUS
        and p1050.get("attempt_name") == T302_ATTEMPT_NAME
        and p1050.get("attempt_exported_on_current_repo_state") is True
        and p1050.get("uses_noncyclic_guardrail") is True
        and T302_ATTEMPT_NAME in t302_text
        and "attempt_must_not_reenter_exhausted_pair12_entry_point_same_lane_as_primary_descent := yes" in t302_text
    )

    noncyclic_guardrail_remains_frozen = (
        f949.get("status") == F949_STATUS
        and f949.get("same_lane_exhaustion_boundary_reached") is True
    )

    exact_lower_route_local_support_family_is_frozen = (
        p766.get("status") == P766_STATUS
        and p766.get("t220_attempt_exported_on_current_repo_state") is True
        and p766.get("first_actual_t218_interface_realization_attempt", {}).get("immediate_missing_subinterface")
        == EXACT_NAMED_SUBINTERFACE
        and p767.get("status") == P767_STATUS
        and p767.get("current_t220_attempt_stalls_exactly_at_the_named_missing_subinterface") is True
        and p767.get("exact_named_missing_subinterface") == EXACT_NAMED_SUBINTERFACE
        and p768.get("status") == P768_STATUS
        and p768.get("t222_target_exported_on_current_repo_state") is True
        and p768.get("current_t222_target_is_future_route_only") is True
        and p769.get("status") == P769_STATUS
        and p769.get("current_repo_still_does_not_export_actual_realization_of_t222_target") is True
        and p770.get("status") == P770_STATUS
        and p770.get("t224_attempt_exported_on_current_repo_state") is True
        and p770.get("first_actual_t222_subinterface_realization_attempt", {}).get("targeted_subinterface")
        == EXACT_NAMED_SUBINTERFACE
    )

    no_current_t302_verdict_or_exact_lower_supplier_boundary_export_found = len(followup_hits) == 0

    t302_attempt_still_has_neither_verdict_nor_exact_lower_supplier_boundary = (
        t302_exact_attempt_is_frozen
        and noncyclic_guardrail_remains_frozen
        and exact_lower_route_local_support_family_is_frozen
        and no_current_t302_verdict_or_exact_lower_supplier_boundary_export_found
    )

    add_check(
        "t302_exact_attempt_is_frozen",
        t302_exact_attempt_is_frozen,
        True,
        "P1050/T302 already freeze one exact first noncyclic actual-realization attempt over the exact F947 source-side input-leg target.",
    )
    add_check(
        "noncyclic_guardrail_remains_frozen",
        noncyclic_guardrail_remains_frozen,
        True,
        "F949 still freezes that the old pair12 same-lane descent is exhausted as a primary route.",
    )
    add_check(
        "exact_lower_route_local_support_family_is_frozen",
        exact_lower_route_local_support_family_is_frozen,
        True,
        "P766/P767/P768/P769/P770 already freeze one exact older route-local lower support family at the named chart-label-retaining pair12 typed seed subinterface.",
    )
    add_check(
        "no_current_t302_verdict_or_exact_lower_supplier_boundary_export_found",
        no_current_t302_verdict_or_exact_lower_supplier_boundary_export_found,
        True,
        "No current export yet upgrades the exact T302 attempt into either a lawful verdict or one exact lower supplier-boundary frozen explicitly under T302.",
    )
    add_check(
        "t302_attempt_still_has_neither_verdict_nor_exact_lower_supplier_boundary",
        t302_attempt_still_has_neither_verdict_nor_exact_lower_supplier_boundary,
        True,
        "Therefore the exact T302 attempt still has neither lawful verdict nor exact lower supplier-boundary on the current repo state.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_LOWER_SUPPLIER_BOUNDARY_NONEXPORT_AUDITED"
        if not blocking and t302_attempt_still_has_neither_verdict_nor_exact_lower_supplier_boundary
        else "FAIL_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_LOWER_SUPPLIER_BOUNDARY_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P1051",
        "status": status,
        "as_of": AS_OF,
        "t302_attempt_name": T302_ATTEMPT_NAME,
        "t220_attempt_name": T220_ATTEMPT_NAME,
        "t224_attempt_name": T224_ATTEMPT_NAME,
        "exact_named_lower_route_local_subinterface": EXACT_NAMED_SUBINTERFACE,
        "t302_verdict_exported_on_current_repo_state": False,
        "t302_exact_lower_supplier_boundary_exported_on_current_repo_state": not no_current_t302_verdict_or_exact_lower_supplier_boundary_export_found,
        "current_repo_already_exports_t302_verdict_or_exact_lower_supplier_boundary_hits": followup_hits,
        "next_honest_move_is_freeze_exact_lower_supplier_boundary_target_below_t302": t302_attempt_still_has_neither_verdict_nor_exact_lower_supplier_boundary,
        "no_false_pass": True,
        "checks": checks,
        "blocking_checks": blocking,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "t302_attempt_name": artifact["t302_attempt_name"],
        "t302_verdict_exported_on_current_repo_state": artifact["t302_verdict_exported_on_current_repo_state"],
        "t302_exact_lower_supplier_boundary_exported_on_current_repo_state": artifact[
            "t302_exact_lower_supplier_boundary_exported_on_current_repo_state"
        ],
        "exact_named_lower_route_local_subinterface": artifact[
            "exact_named_lower_route_local_subinterface"
        ],
        "next_honest_move_is_freeze_exact_lower_supplier_boundary_target_below_t302": artifact[
            "next_honest_move_is_freeze_exact_lower_supplier_boundary_target_below_t302"
        ],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
