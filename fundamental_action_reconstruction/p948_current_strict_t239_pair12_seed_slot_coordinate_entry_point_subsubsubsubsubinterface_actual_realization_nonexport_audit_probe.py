#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-21"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P742 = GENERATED / "p742_current_strict_t196_actual_source_topology_selector_witness_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_P747 = GENERATED / "p747_current_strict_t201_actual_source_topology_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge_nonexport_audit_probe_summary.json"
IN_P783 = GENERATED / "p783_current_strict_t237_pair12_seed_slot_coordinate_entry_subsubsubsubinterface_actual_realization_attempt_immediate_missing_subsubsubsubsubinterface_nonexport_audit_probe_summary.json"
IN_P784 = GENERATED / "p784_current_strict_t238_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_target_probe_summary.json"
IN_T236 = ROOT / "T236_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_SUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"
IN_T238 = ROOT / "T238_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p948_current_strict_t239_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p948_current_strict_t239_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_nonexport_audit_probe_summary.json"

T239_TARGET_NAME = "Pair12SeedSlotCoordinateEntryPointSubsubsubsubsubinterface_strict_v1"
T238_TARGET_SYMBOL = "W_strict_t173_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_target_v1"
T236_ATTEMPT_SYMBOL = "W_strict_t173_pair12_seed_slot_coordinate_entry_subsubsubsubinterface_actual_realization_attempt_v1"
EXACT_SUBSUBSUBSUBINTERFACE_NAME = (
    "chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_on_Sigma_sel_src_target_v1_"
    "prior_to_surviving_F301_pair12_carrier_binding_and_prior_to_Q_basis_sel_v1_"
    "terminal_collapse_and_prior_to_projector_only_local_pair12_atlas_collapse"
)
EXACT_SUBSUBSUBSUBSUBINTERFACE_NAME = (
    "chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_point_on_Sigma_sel_src_target_v1_"
    "prior_to_surviving_F301_pair12_carrier_binding_and_prior_to_Q_basis_sel_v1_"
    "terminal_collapse_and_prior_to_projector_only_local_pair12_atlas_collapse"
)
CURRENT_THEOREM_FILE = (
    "N781_CURRENT_STRICT_T239_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_actual_subsubsubsubsubinterface_realization_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "T236_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_SUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p782_current_strict_t236_pair12_seed_slot_coordinate_entry_subsubsubsubinterface_actual_realization_attempt_probe.py",
        "N778_CURRENT_STRICT_T236_PAIR12_SEED_SLOT_COORDINATE_ENTRY_SUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p783_current_strict_t237_pair12_seed_slot_coordinate_entry_subsubsubsubinterface_actual_realization_attempt_immediate_missing_subsubsubsubsubinterface_nonexport_audit_probe.py",
        "N779_CURRENT_STRICT_T237_PAIR12_SEED_SLOT_COORDINATE_ENTRY_SUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBSUBSUBSUBSUBINTERFACE_NONEXPORT_AUDIT_THEOREM.md",
        "T238_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_TARGET_SPEC.md",
        "p784_current_strict_t238_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_target_probe.py",
        "N780_CURRENT_STRICT_T238_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_TARGET_THEOREM.md",
        "T240_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p949_current_strict_t240_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_probe.py",
        "N782_CURRENT_STRICT_T240_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
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
                EXACT_SUBSUBSUBSUBSUBINTERFACE_NAME in text
                and "Sigma_sel_src_target_v1" in text
                and ("pair1/pair2" in text or "F301" in text)
            ):
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P742, IN_P747, IN_P783, IN_P784, IN_T236, IN_T238]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P948",
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
    p783 = load_json(IN_P783)
    p784 = load_json(IN_P784)
    t236_text = load_text(IN_T236)
    t238_text = load_text(IN_T238)
    positive_actual_subsubsubsubsubinterface_realization_candidates = (
        scan_positive_actual_subsubsubsubsubinterface_realization_candidates()
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

    p784_t238_target_already_exported_only_at_future_only_strength = (
        bool(p784.get("t238_target_exported_on_current_repo_state"))
        and bool(p784.get("current_t238_target_is_future_route_only"))
        and bool(p784.get("current_t238_target_freezes_exact_t236_immediate_missing_subsubsubsubsubinterface"))
        and bool(
            p784.get(
                "current_t238_target_remains_below_actual_subsubsubsubsubinterface_export_subsubsubsubinterface_export_subsubsubinterface_export_interface_export_and_t176_discharge"
            )
        )
        and bool(
            p784.get(
                "next_honest_move_is_actual_export_of_frozen_exact_missing_subsubsubsubsubinterface_or_exact_failure_localization_below_it"
            )
        )
    )

    p783_exact_missing_subsubsubsubsubinterface_nonexport_boundary_already_exported = (
        bool(p783.get("t237_boundary_exported_on_current_repo_state"))
        and bool(
            p783.get(
                "current_repo_still_does_not_export_actual_realization_of_t236_immediate_missing_subsubsubsubsubinterface"
            )
        )
        and bool(
            p783.get("current_t236_attempt_stalls_exactly_at_the_named_missing_subsubsubsubsubinterface")
        )
        and str(p783.get("exact_named_missing_subsubsubsubsubinterface") or "")
        == EXACT_SUBSUBSUBSUBSUBINTERFACE_NAME
    )

    current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_point_subsubsubsubsubinterface = (
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

    current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_entry_point_subsubsubsubsubinterface = (
        bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_exported"))
        and bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe"))
        and bool(p747.get("current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas"))
        and not bool(
            p747.get("current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge")
        )
    )

    t236_t238_same_exact_subsubsubsubsubinterface_route_still_frozen = all(
        needle in t236_text
        for needle in [
            T236_ATTEMPT_SYMBOL,
            EXACT_SUBSUBSUBSUBINTERFACE_NAME,
            "attempt_precedes_surviving_F301_pair12_carrier_binding := yes",
            "attempt_precedes_Q_basis_sel_v1_terminal_collapse := yes",
            "attempt_precedes_projector_only_local_pair12_atlas_collapse := yes",
            "attempt_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "attempt_must_remain_below_success_verdict := yes",
        ]
    ) and all(
        needle in t238_text
        for needle in [
            T238_TARGET_SYMBOL,
            EXACT_SUBSUBSUBSUBSUBINTERFACE_NAME,
            "target_precedes_surviving_F301_pair12_carrier_binding := yes",
            "target_precedes_Q_basis_sel_v1_terminal_collapse := yes",
            "target_precedes_projector_only_local_pair12_atlas_collapse := yes",
            "target_remains_below_actual_subsubsubsubsubinterface_export := yes",
        ]
    )

    current_repo_has_exported_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_point_subsubsubsubsubinterface = (
        len(positive_actual_subsubsubsubsubinterface_realization_candidates) > 0
    )

    current_repo_still_does_not_export_actual_realization_of_t238_target = (
        p784_t238_target_already_exported_only_at_future_only_strength
        and p783_exact_missing_subsubsubsubsubinterface_nonexport_boundary_already_exported
        and current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_point_subsubsubsubsubinterface
        and current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_entry_point_subsubsubsubsubinterface
        and t236_t238_same_exact_subsubsubsubsubinterface_route_still_frozen
        and not current_repo_has_exported_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_point_subsubsubsubsubinterface
        and len(positive_actual_subsubsubsubsubinterface_realization_candidates) == 0
    )

    next_honest_move_is_actual_t238_subsubsubsubsubinterface_realization_attempt_or_exact_failure_localization_below_it = (
        current_repo_still_does_not_export_actual_realization_of_t238_target
    )

    add_check(
        "p784_t238_target_already_exported_only_at_future_only_strength",
        p784_t238_target_already_exported_only_at_future_only_strength,
        True,
        "P784 already freezes that T238 exists only at future-only target strength for the exact missing chart-label-retaining pair1/pair2 typed seed-slot coordinate entry point subsubsubsubsubinterface.",
    )
    add_check(
        "p783_exact_missing_subsubsubsubsubinterface_nonexport_boundary_already_exported",
        p783_exact_missing_subsubsubsubsubinterface_nonexport_boundary_already_exported,
        True,
        "P783 already freezes that the same exact T236 missing subsubsubsubsubinterface remains unexported on the current repo state.",
    )
    add_check(
        "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_point_subsubsubsubsubinterface",
        current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_point_subsubsubsubsubinterface,
        True,
        "The strongest current exported codomain continuation out of Sigma_sel_src_target_v1 still does not export one actual chart-label-retaining pair1/pair2 typed seed-slot coordinate entry point subsubsubsubsubinterface (P742).",
    )
    add_check(
        "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_entry_point_subsubsubsubsubinterface",
        current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_entry_point_subsubsubsubsubinterface,
        True,
        "The strongest current local pair1/pair2 atlas lane still does not export one nonprojector chart-label-retaining seed-slot coordinate entry point subsubsubsubsubinterface (P747).",
    )
    add_check(
        "t236_t238_same_exact_subsubsubsubsubinterface_route_still_frozen",
        t236_t238_same_exact_subsubsubsubsubinterface_route_still_frozen,
        True,
        "T236 and T238 still freeze the same exact lower seed-slot coordinate entry point subsubsubsubsubinterface on the same exact attempt route.",
    )
    add_check(
        "current_repo_has_exported_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_point_subsubsubsubsubinterface",
        current_repo_has_exported_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_point_subsubsubsubsubinterface,
        False,
        "No current repo artifact positively exports one actual realization of the exact chart-label-retaining pair1/pair2 typed seed-slot coordinate entry point subsubsubsubsubinterface named below T238.",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t238_target",
        current_repo_still_does_not_export_actual_realization_of_t238_target,
        True,
        "Therefore the current repo still does not export one actual realization of the T238 target.",
    )
    add_check(
        "next_honest_move_is_actual_t238_subsubsubsubsubinterface_realization_attempt_or_exact_failure_localization_below_it",
        next_honest_move_is_actual_t238_subsubsubsubsubinterface_realization_attempt_or_exact_failure_localization_below_it,
        True,
        "Hence the next honest move is now one actual realization attempt of T238, or, only if that same route later stalls, exact failure localization below it.",
    )

    status = (
        "PASS_STRICT_T239_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking and current_repo_still_does_not_export_actual_realization_of_t238_target
        else "FAIL_STRICT_T239_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P948",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "theorem_result": {
            "t239_target_name": T239_TARGET_NAME,
            "t239_target_exported_on_current_repo_state": False,
            "current_repo_still_does_not_export_actual_realization_of_t238_target": current_repo_still_does_not_export_actual_realization_of_t238_target,
            "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_point_subsubsubsubsubinterface": current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_point_subsubsubsubsubinterface,
            "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_entry_point_subsubsubsubsubinterface": current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_entry_point_subsubsubsubsubinterface,
            "current_exact_t236_missing_subsubsubsubsubinterface_still_only_future_target_not_actual_export": current_repo_still_does_not_export_actual_realization_of_t238_target,
            "next_honest_move_is_actual_t238_subsubsubsubsubinterface_realization_attempt_or_exact_failure_localization_below_it": next_honest_move_is_actual_t238_subsubsubsubsubinterface_realization_attempt_or_exact_failure_localization_below_it,
            "no_false_pass": True,
        },
        "positive_actual_subsubsubsubsubinterface_realization_candidates": positive_actual_subsubsubsubsubinterface_realization_candidates,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P948",
        "status": status,
        "as_of": AS_OF,
        "t239_target_name": T239_TARGET_NAME,
        "t239_target_exported_on_current_repo_state": False,
        "current_repo_still_does_not_export_actual_realization_of_t238_target": current_repo_still_does_not_export_actual_realization_of_t238_target,
        "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_point_subsubsubsubsubinterface": current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_entry_point_subsubsubsubsubinterface,
        "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_entry_point_subsubsubsubsubinterface": current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_entry_point_subsubsubsubsubinterface,
        "current_exact_t236_missing_subsubsubsubsubinterface_still_only_future_target_not_actual_export": current_repo_still_does_not_export_actual_realization_of_t238_target,
        "next_honest_move_is_actual_t238_subsubsubsubsubinterface_realization_attempt_or_exact_failure_localization_below_it": next_honest_move_is_actual_t238_subsubsubsubsubinterface_realization_attempt_or_exact_failure_localization_below_it,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
