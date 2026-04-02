#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P742 = GENERATED / "p742_current_strict_t196_actual_source_topology_selector_witness_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_P747 = GENERATED / "p747_current_strict_t201_actual_source_topology_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge_nonexport_audit_probe_summary.json"
IN_P775 = GENERATED / "p775_current_strict_t229_pair12_seed_slot_subsubinterface_actual_realization_attempt_immediate_missing_subsubsubinterface_nonexport_audit_probe_summary.json"
IN_P776 = GENERATED / "p776_current_strict_t230_pair12_seed_slot_coordinate_subsubsubinterface_target_probe_summary.json"
IN_T228 = ROOT / "T228_CURRENT_STRICT_PAIR12_SEED_SLOT_SUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"
IN_T230 = ROOT / "T230_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_SUBSUBSUBINTERFACE_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p777_current_strict_t231_pair12_seed_slot_coordinate_subsubsubinterface_actual_realization_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p777_current_strict_t231_pair12_seed_slot_coordinate_subsubsubinterface_actual_realization_nonexport_audit_probe_summary.json"

T231_TARGET_NAME = "Pair12SeedSlotCoordinateSubsubsubinterface_strict_v1"
T230_TARGET_SYMBOL = "W_strict_t173_pair12_seed_slot_coordinate_subsubsubinterface_target_v1"
T228_ATTEMPT_SYMBOL = "W_strict_t173_pair12_seed_slot_subsubinterface_actual_realization_attempt_v1"
EXACT_SUBSUBINTERFACE_NAME = (
    "chart_label_retaining_pair12_typed_seed_slot_on_Sigma_sel_src_target_v1_"
    "prior_to_surviving_F301_pair12_carrier_binding_and_prior_to_Q_basis_sel_v1_"
    "terminal_collapse_and_prior_to_projector_only_local_pair12_atlas_collapse"
)
EXACT_SUBSUBSUBINTERFACE_NAME = (
    "chart_label_retaining_pair12_typed_seed_slot_coordinate_on_Sigma_sel_src_target_v1_"
    "prior_to_surviving_F301_pair12_carrier_binding_and_prior_to_Q_basis_sel_v1_"
    "terminal_collapse_and_prior_to_projector_only_local_pair12_atlas_collapse"
)
CURRENT_THEOREM_FILE = (
    "N773_CURRENT_STRICT_T231_PAIR12_SEED_SLOT_COORDINATE_SUBSUBSUBINTERFACE_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_actual_subsubsubinterface_realization_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "T228_CURRENT_STRICT_PAIR12_SEED_SLOT_SUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p774_current_strict_t228_pair12_seed_slot_subsubinterface_actual_realization_attempt_probe.py",
        "N770_CURRENT_STRICT_T228_PAIR12_SEED_SLOT_SUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p775_current_strict_t229_pair12_seed_slot_subsubinterface_actual_realization_attempt_immediate_missing_subsubsubinterface_nonexport_audit_probe.py",
        "N771_CURRENT_STRICT_T229_PAIR12_SEED_SLOT_SUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBSUBSUBINTERFACE_NONEXPORT_AUDIT_THEOREM.md",
        "T230_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_SUBSUBSUBINTERFACE_TARGET_SPEC.md",
        "p776_current_strict_t230_pair12_seed_slot_coordinate_subsubsubinterface_target_probe.py",
        "N772_CURRENT_STRICT_T230_PAIR12_SEED_SLOT_COORDINATE_SUBSUBSUBINTERFACE_TARGET_THEOREM.md",
        "T230_CURRENT_STRICT_PAIR12_SEED_SLOT_SUBSUBSUBINTERFACE_TARGET_SPEC.md",
        "p776_current_strict_t230_pair12_seed_slot_subsubsubinterface_target_probe.py",
        "N772_CURRENT_STRICT_T230_PAIR12_SEED_SLOT_SUBSUBSUBINTERFACE_TARGET_THEOREM.md",
        "T232_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_SUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p778_current_strict_t232_pair12_seed_slot_coordinate_subsubsubinterface_actual_realization_attempt_probe.py",
        "N774_CURRENT_STRICT_T232_PAIR12_SEED_SLOT_COORDINATE_SUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "T232_CURRENT_STRICT_PAIR12_SEED_SLOT_SUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p778_current_strict_t232_pair12_seed_slot_subsubsubinterface_actual_realization_attempt_probe.py",
        "N774_CURRENT_STRICT_T232_PAIR12_SEED_SLOT_SUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p779_current_strict_t233_pair12_seed_slot_coordinate_subsubsubinterface_actual_realization_attempt_immediate_missing_subsubsubsubinterface_nonexport_audit_probe.py",
        "N775_CURRENT_STRICT_T233_PAIR12_SEED_SLOT_COORDINATE_SUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBSUBSUBSUBINTERFACE_NONEXPORT_AUDIT_THEOREM.md",
        "T234_CURRENT_STRICT_PAIR12_SEED_SLOT_COORDINATE_ENTRY_SUBSUBSUBSUBINTERFACE_TARGET_SPEC.md",
        "p780_current_strict_t234_pair12_seed_slot_coordinate_entry_subsubsubsubinterface_target_probe.py",
        "N776_CURRENT_STRICT_T234_PAIR12_SEED_SLOT_COORDINATE_ENTRY_SUBSUBSUBSUBINTERFACE_TARGET_THEOREM.md",
        "T234_CURRENT_STRICT_PAIR12_SEED_SLOT_ENTRY_SUBSUBSUBSUBINTERFACE_TARGET_SPEC.md",
        "p780_current_strict_t234_pair12_seed_slot_entry_subsubsubsubinterface_target_probe.py",
        "N776_CURRENT_STRICT_T234_PAIR12_SEED_SLOT_ENTRY_SUBSUBSUBSUBINTERFACE_TARGET_THEOREM.md",
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
                EXACT_SUBSUBSUBINTERFACE_NAME in text
                and "Sigma_sel_src_target_v1" in text
                and ("pair1/pair2" in text or "F301" in text)
            ):
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P742, IN_P747, IN_P775, IN_P776, IN_T228, IN_T230]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P777",
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
    p775 = load_json(IN_P775)
    p776 = load_json(IN_P776)
    t228_text = load_text(IN_T228)
    t230_text = load_text(IN_T230)
    positive_actual_subsubsubinterface_realization_candidates = (
        scan_positive_actual_subsubsubinterface_realization_candidates()
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

    p776_t230_target_already_exported_only_at_future_only_strength = (
        bool(p776.get("t230_target_exported_on_current_repo_state"))
        and bool(p776.get("current_t230_target_is_future_route_only"))
        and bool(p776.get("current_t230_target_freezes_exact_t228_immediate_missing_subsubsubinterface"))
        and bool(
            p776.get(
                "current_t230_target_remains_below_actual_subsubsubinterface_export_subsubinterface_export_interface_export_and_t176_discharge"
            )
        )
        and bool(
            p776.get(
                "next_honest_move_is_actual_export_of_frozen_exact_missing_subsubsubinterface_or_exact_failure_localization_below_it"
            )
        )
    )

    p775_exact_missing_subsubsubinterface_nonexport_boundary_already_exported = (
        bool(p775.get("t229_boundary_exported_on_current_repo_state"))
        and bool(
            p775.get(
                "current_repo_still_does_not_export_actual_realization_of_t228_immediate_missing_subsubsubinterface"
            )
        )
        and bool(p775.get("current_t228_attempt_stalls_exactly_at_the_named_missing_subsubsubinterface"))
        and str(p775.get("exact_named_missing_subsubsubinterface") or "") == EXACT_SUBSUBSUBINTERFACE_NAME
    )

    current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_subsubsubinterface = (
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

    current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_subsubsubinterface = (
        bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_exported"))
        and bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe"))
        and bool(p747.get("current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas"))
        and not bool(
            p747.get("current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge")
        )
    )

    t228_t230_same_exact_subsubsubinterface_route_still_frozen = all(
        needle in t228_text
        for needle in [
            T228_ATTEMPT_SYMBOL,
            EXACT_SUBSUBINTERFACE_NAME,
            "attempt_precedes_surviving_F301_pair12_carrier_binding := yes",
            "attempt_precedes_Q_basis_sel_v1_terminal_collapse := yes",
            "attempt_precedes_projector_only_local_pair12_atlas_collapse := yes",
            "attempt_retains_branch_relevance_to_delta_k_vs_delta_minus_k := yes",
            "attempt_must_remain_below_success_verdict := yes",
        ]
    ) and all(
        needle in t230_text
        for needle in [
            T230_TARGET_SYMBOL,
            EXACT_SUBSUBSUBINTERFACE_NAME,
            "target_precedes_surviving_F301_pair12_carrier_binding := yes",
            "target_precedes_Q_basis_sel_v1_terminal_collapse := yes",
            "target_precedes_projector_only_local_pair12_atlas_collapse := yes",
            "target_remains_below_actual_subsubsubinterface_export := yes",
        ]
    )

    current_repo_has_exported_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_subsubsubinterface = (
        len(positive_actual_subsubsubinterface_realization_candidates) > 0
    )

    current_repo_still_does_not_export_actual_realization_of_t230_target = (
        p776_t230_target_already_exported_only_at_future_only_strength
        and p775_exact_missing_subsubsubinterface_nonexport_boundary_already_exported
        and current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_subsubsubinterface
        and current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_subsubsubinterface
        and t228_t230_same_exact_subsubsubinterface_route_still_frozen
        and not current_repo_has_exported_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_subsubsubinterface
        and len(positive_actual_subsubsubinterface_realization_candidates) == 0
    )

    next_honest_move_is_actual_t230_subsubsubinterface_realization_attempt_or_exact_failure_localization_below_it = (
        current_repo_still_does_not_export_actual_realization_of_t230_target
    )

    add_check(
        "p776_t230_target_already_exported_only_at_future_only_strength",
        p776_t230_target_already_exported_only_at_future_only_strength,
        True,
        "P776 already freezes that T230 exists only at future-only target strength for the exact missing chart-label-retaining pair1/pair2 typed seed-slot coordinate subsubsubinterface.",
    )
    add_check(
        "p775_exact_missing_subsubsubinterface_nonexport_boundary_already_exported",
        p775_exact_missing_subsubsubinterface_nonexport_boundary_already_exported,
        True,
        "P775 already freezes that the same exact T228 missing subsubsubinterface remains unexported on the current repo state.",
    )
    add_check(
        "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_subsubsubinterface",
        current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_subsubsubinterface,
        True,
        "The strongest current exported codomain continuation out of Sigma_sel_src_target_v1 still does not export one actual chart-label-retaining pair1/pair2 typed seed-slot coordinate subsubsubinterface (P742).",
    )
    add_check(
        "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_subsubsubinterface",
        current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_subsubsubinterface,
        True,
        "The strongest current local pair1/pair2 atlas-side lane still does not export one nonprojector chart-label-retaining seed-slot coordinate subsubsubinterface (P747).",
    )
    add_check(
        "t228_t230_same_exact_subsubsubinterface_route_still_frozen",
        t228_t230_same_exact_subsubsubinterface_route_still_frozen,
        True,
        "T228 and T230 still freeze the same exact lower seed-slot coordinate subsubsubinterface on the same exact attempt route.",
    )
    add_check(
        "current_repo_has_exported_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_subsubsubinterface",
        current_repo_has_exported_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_subsubsubinterface,
        False,
        "No current repo artifact positively exports one actual realization of the exact chart-label-retaining pair1/pair2 typed seed-slot coordinate subsubsubinterface named below T230.",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t230_target",
        current_repo_still_does_not_export_actual_realization_of_t230_target,
        True,
        "Therefore the current repo still does not export one actual realization of the T230 target.",
    )
    add_check(
        "next_honest_move_is_actual_t230_subsubsubinterface_realization_attempt_or_exact_failure_localization_below_it",
        next_honest_move_is_actual_t230_subsubsubinterface_realization_attempt_or_exact_failure_localization_below_it,
        True,
        "Hence the next honest move is now one actual realization attempt of T230, or, only if that same route later stalls, exact failure localization below it.",
    )

    status = (
        "PASS_STRICT_T231_PAIR12_SEED_SLOT_COORDINATE_SUBSUBSUBINTERFACE_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking and current_repo_still_does_not_export_actual_realization_of_t230_target
        else "FAIL_STRICT_T231_PAIR12_SEED_SLOT_COORDINATE_SUBSUBSUBINTERFACE_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P777",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "theorem_result": {
            "t231_target_name": T231_TARGET_NAME,
            "t231_target_exported_on_current_repo_state": False,
            "current_repo_still_does_not_export_actual_realization_of_t230_target": current_repo_still_does_not_export_actual_realization_of_t230_target,
            "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_subsubsubinterface": current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_subsubsubinterface,
            "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_subsubsubinterface": current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_subsubsubinterface,
            "current_exact_t228_missing_subsubsubinterface_still_only_future_target_not_actual_export": current_repo_still_does_not_export_actual_realization_of_t230_target,
            "next_honest_move_is_actual_t230_subsubsubinterface_realization_attempt_or_exact_failure_localization_below_it": next_honest_move_is_actual_t230_subsubsubinterface_realization_attempt_or_exact_failure_localization_below_it,
            "no_false_pass": True,
        },
        "positive_actual_subsubsubinterface_realization_candidates": positive_actual_subsubsubinterface_realization_candidates,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P777",
        "status": status,
        "as_of": AS_OF,
        "t231_target_name": T231_TARGET_NAME,
        "t231_target_exported_on_current_repo_state": False,
        "current_repo_still_does_not_export_actual_realization_of_t230_target": current_repo_still_does_not_export_actual_realization_of_t230_target,
        "current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_subsubsubinterface": current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_slot_coordinate_subsubsubinterface,
        "current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_subsubsubinterface": current_local_pair12_chart_sensitive_atlas_lane_still_lacks_nonprojector_seed_slot_coordinate_subsubsubinterface,
        "current_exact_t228_missing_subsubsubinterface_still_only_future_target_not_actual_export": current_repo_still_does_not_export_actual_realization_of_t230_target,
        "next_honest_move_is_actual_t230_subsubsubinterface_realization_attempt_or_exact_failure_localization_below_it": next_honest_move_is_actual_t230_subsubsubinterface_realization_attempt_or_exact_failure_localization_below_it,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
