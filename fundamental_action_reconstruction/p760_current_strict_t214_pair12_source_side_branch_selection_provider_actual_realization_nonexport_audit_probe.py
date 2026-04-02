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

IN_P759 = GENERATED / "p759_current_strict_t213_pair12_source_side_branch_selection_provider_target_probe_summary.json"
IN_P741 = GENERATED / "p741_current_strict_t195_actual_source_topology_selector_witness_pair12_witness_split_promotion_bridge_nonexport_audit_probe_summary.json"
IN_P742 = GENERATED / "p742_current_strict_t196_actual_source_topology_selector_witness_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_P743 = GENERATED / "p743_current_strict_t197_source_topology_quotient_safe_qw2191_resolution_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_P747 = GENERATED / "p747_current_strict_t201_actual_source_topology_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge_nonexport_audit_probe_summary.json"
IN_T213 = ROOT / "T213_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p760_current_strict_t214_pair12_source_side_branch_selection_provider_actual_realization_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p760_current_strict_t214_pair12_source_side_branch_selection_provider_actual_realization_nonexport_audit_probe_summary.json"

ACTUAL_BRANCH_SELECTION_PROVIDER = (
    "W_strict_t173_pair12_source_side_branch_selection_provider_v1"
)
CURRENT_THEOREM_FILE = (
    "N756_CURRENT_STRICT_T214_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_actual_realization_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "T213_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_TARGET_SPEC.md",
        "p759_current_strict_t213_pair12_source_side_branch_selection_provider_target_probe.py",
        "N755_CURRENT_STRICT_T213_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_TARGET_THEOREM.md",
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
                ACTUAL_BRANCH_SELECTION_PROVIDER in text
                and "tau_src_candidate_v1" in text
                and "pair1/pair2" in text
            ):
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P759, IN_P741, IN_P742, IN_P743, IN_P747, IN_T213]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P760",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p759 = load_json(IN_P759)
    p741 = load_json(IN_P741)
    p742 = load_json(IN_P742)
    p743 = load_json(IN_P743)
    p747 = load_json(IN_P747)
    t213_text = load_text(IN_T213)
    positive_actual_realization_candidates = scan_positive_actual_realization_candidates()

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

    p759_t213_target_already_exported_only_at_future_only_strength = (
        bool(p759.get("t213_target_exported_on_current_repo_state"))
        and bool(p759.get("current_t213_target_is_source_side_observer_free"))
        and bool(p759.get("current_t213_target_is_pair12_typed_and_branch_sensitive"))
        and bool(p759.get("current_t213_target_is_chart_sensitive_and_residual_datum_carrier_binding"))
        and bool(p759.get("current_t213_target_is_external_to_current_exported_p731_continuation_family"))
        and bool(p759.get("current_t213_target_is_nonconvention_nonpremise_based"))
        and bool(p759.get("current_t213_target_remains_below_t176_discharge_and_physical_sign_datum_claim"))
        and bool(p759.get("current_t213_target_is_future_route_only"))
    )

    current_actual_source_topology_selector_witness_still_not_pair12_typed = (
        not bool(p741.get("t195_target_exported_on_current_repo_state"))
        and bool(p741.get("current_actual_source_topology_selector_witness_exported"))
        and bool(p741.get("current_actual_source_topology_selector_witness_binds_same_tau_src_packet_as_pair12_carrier"))
        and bool(p741.get("current_actual_source_topology_selector_witness_is_chart_bound_preobserver_only"))
        and bool(p741.get("current_actual_source_topology_selector_witness_remains_prelm_not_pair12_typed"))
        and not bool(p741.get("p731_pair12_witness_split_descends_to_current_actual_source_topology_selector_witness"))
    )

    current_actual_selector_witness_codomain_still_lacks_pair12_typed_carrier_bridge = (
        not bool(p742.get("t196_target_exported_on_current_repo_state"))
        and bool(p742.get("current_actual_selector_witness_codomain_has_exported_basis_free_chart_label_forgetting_continuation"))
        and bool(p742.get("surviving_pair12_residual_datum_carrier_remains_selector_neutral"))
        and bool(p742.get("current_actual_selector_witness_codomain_continuation_remains_basis_free_not_pair12_typed"))
        and not bool(p742.get("current_actual_selector_witness_codomain_has_exported_pair12_typed_residual_datum_continuation"))
        and not bool(
            p742.get("p731_pair12_witness_split_descends_to_current_actual_selector_witness_to_residual_datum_typed_carrier_bridge")
        )
    )

    current_qw2191_safe_resolution_still_lacks_pair12_typed_branch_provider = (
        not bool(p743.get("t197_target_exported_on_current_repo_state"))
        and bool(p743.get("current_actual_source_topology_quotient_safe_qw2191_resolution_exported"))
        and bool(p743.get("current_actual_source_topology_quotient_safe_qw2191_resolution_binds_same_tau_src_packet_as_pair12_carrier"))
        and bool(p743.get("current_actual_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only"))
        and bool(p743.get("current_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only_not_pair12_typed"))
        and not bool(p743.get("current_source_topology_quotient_safe_qw2191_resolution_has_exported_pair12_typed_residual_datum_continuation"))
        and not bool(
            p743.get("p731_pair12_witness_split_descends_to_current_source_topology_quotient_safe_qw2191_resolution_to_residual_datum_typed_bridge")
        )
    )

    current_selector_witness_target_still_lacks_local_chart_sensitive_pair12_bridge = (
        not bool(p747.get("t201_target_exported_on_current_repo_state"))
        and bool(p747.get("current_actual_source_topology_selector_witness_target_exported"))
        and bool(p747.get("current_actual_source_topology_selector_witness_target_remains_chart_bound_prelm"))
        and bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_exported"))
        and bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe"))
        and bool(p747.get("current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas"))
        and not bool(p747.get("current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge"))
        and not bool(
            p747.get("p731_pair12_witness_split_descends_to_current_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge")
        )
    )

    t213_future_only_target_symbol_frozen = (
        "W_strict_t173_pair12_source_side_branch_selection_provider_target_v1" in t213_text
        and "target_is_pair12_typed := yes" in t213_text
        and "target_is_branch_sensitive_to_delta_k_vs_delta_minus_k := yes" in t213_text
        and "future_route_only := yes" in t213_text
    )

    current_repo_has_exported_actual_pair12_source_side_branch_selection_provider = (
        len(positive_actual_realization_candidates) > 0
    )

    current_repo_still_does_not_export_actual_realization_of_t213_target = (
        p759_t213_target_already_exported_only_at_future_only_strength
        and current_actual_source_topology_selector_witness_still_not_pair12_typed
        and current_actual_selector_witness_codomain_still_lacks_pair12_typed_carrier_bridge
        and current_qw2191_safe_resolution_still_lacks_pair12_typed_branch_provider
        and current_selector_witness_target_still_lacks_local_chart_sensitive_pair12_bridge
        and t213_future_only_target_symbol_frozen
        and not current_repo_has_exported_actual_pair12_source_side_branch_selection_provider
        and len(positive_actual_realization_candidates) == 0
    )

    next_honest_move_is_actual_t213_realization_attempt_or_further_provider_attack = (
        current_repo_still_does_not_export_actual_realization_of_t213_target
    )

    add_check(
        "p759_t213_target_already_exported_only_at_future_only_strength",
        {
            "t213_target_exported_on_current_repo_state": bool(
                p759.get("t213_target_exported_on_current_repo_state")
            ),
            "current_t213_target_is_source_side_observer_free": bool(
                p759.get("current_t213_target_is_source_side_observer_free")
            ),
            "current_t213_target_is_pair12_typed_and_branch_sensitive": bool(
                p759.get("current_t213_target_is_pair12_typed_and_branch_sensitive")
            ),
            "current_t213_target_is_chart_sensitive_and_residual_datum_carrier_binding": bool(
                p759.get("current_t213_target_is_chart_sensitive_and_residual_datum_carrier_binding")
            ),
            "current_t213_target_is_external_to_current_exported_p731_continuation_family": bool(
                p759.get("current_t213_target_is_external_to_current_exported_p731_continuation_family")
            ),
            "current_t213_target_is_nonconvention_nonpremise_based": bool(
                p759.get("current_t213_target_is_nonconvention_nonpremise_based")
            ),
            "current_t213_target_remains_below_t176_discharge_and_physical_sign_datum_claim": bool(
                p759.get("current_t213_target_remains_below_t176_discharge_and_physical_sign_datum_claim")
            ),
            "current_t213_target_is_future_route_only": bool(
                p759.get("current_t213_target_is_future_route_only")
            ),
        },
        {
            "t213_target_exported_on_current_repo_state": True,
            "current_t213_target_is_source_side_observer_free": True,
            "current_t213_target_is_pair12_typed_and_branch_sensitive": True,
            "current_t213_target_is_chart_sensitive_and_residual_datum_carrier_binding": True,
            "current_t213_target_is_external_to_current_exported_p731_continuation_family": True,
            "current_t213_target_is_nonconvention_nonpremise_based": True,
            "current_t213_target_remains_below_t176_discharge_and_physical_sign_datum_claim": True,
            "current_t213_target_is_future_route_only": True,
        },
        "P759 already freezes that T213 exists only as a future-only target for a genuinely new source-side pair1/pair2 branch-selection provider.",
    )
    add_check(
        "current_actual_source_topology_selector_witness_still_not_pair12_typed",
        current_actual_source_topology_selector_witness_still_not_pair12_typed,
        True,
        "The strongest current actual source-topology selector witness is real, but still remains chart-bound preLM and not pair1/pair2 typed (P741).",
    )
    add_check(
        "current_actual_selector_witness_codomain_still_lacks_pair12_typed_carrier_bridge",
        current_actual_selector_witness_codomain_still_lacks_pair12_typed_carrier_bridge,
        True,
        "The strongest current codomain continuation out of that actual witness still remains basis-free and does not export a pair1/pair2 typed residual-datum bridge (P742).",
    )
    add_check(
        "current_qw2191_safe_resolution_still_lacks_pair12_typed_branch_provider",
        current_qw2191_safe_resolution_still_lacks_pair12_typed_branch_provider,
        True,
        "The strongest current actual quotient-safe QW-2191 resolution still remains quotient-class-only and does not export a pair1/pair2 typed branch provider (P743).",
    )
    add_check(
        "current_selector_witness_target_still_lacks_local_chart_sensitive_pair12_bridge",
        current_selector_witness_target_still_lacks_local_chart_sensitive_pair12_bridge,
        True,
        "The strongest current source-side selector-witness target still remains unbridged to the local chart-sensitive pair1/pair2 atlas lane (P747).",
    )
    add_check(
        "t213_future_only_target_symbol_frozen",
        t213_future_only_target_symbol_frozen,
        True,
        "T213 already names the future-only target object and keeps it explicitly pair1/pair2-typed, branch-sensitive, and below actual realization.",
    )
    add_check(
        "positive_actual_t214_realization_candidates",
        positive_actual_realization_candidates,
        [],
        "No current exported packet realizes one actual source-side pair1/pair2 branch-selection provider corresponding to the T213 target.",
    )
    add_check(
        "current_repo_has_exported_actual_pair12_source_side_branch_selection_provider",
        current_repo_has_exported_actual_pair12_source_side_branch_selection_provider,
        False,
        "The current repo still exports no actual source-side pair1/pair2 branch-selection provider.",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t213_target",
        current_repo_still_does_not_export_actual_realization_of_t213_target,
        True,
        "Therefore the current repo still does not export one actual realization of the future-only T213 provider target.",
    )
    add_check(
        "next_honest_move_is_actual_t213_realization_attempt_or_further_provider_attack",
        next_honest_move_is_actual_t213_realization_attempt_or_further_provider_attack,
        True,
        "Hence the next honest strict move is now either one actual realization attempt of T213 or a further genuinely new provider attack if that target route stalls too.",
    )

    status = (
        "PASS_STRICT_T214_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking
        and current_repo_still_does_not_export_actual_realization_of_t213_target
        and next_honest_move_is_actual_t213_realization_attempt_or_further_provider_attack
        else "FAIL_STRICT_T214_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P760",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "theorem_result": {
            "t214_target_name": "Pair12SourceSideBranchSelectionProvider_global_C_v1_strict_v1",
            "t214_target_exported_on_current_repo_state": False,
            "current_repo_still_does_not_export_actual_realization_of_t213_target": current_repo_still_does_not_export_actual_realization_of_t213_target,
            "current_actual_source_topology_selector_witness_still_not_pair12_typed": current_actual_source_topology_selector_witness_still_not_pair12_typed,
            "current_actual_selector_witness_codomain_still_lacks_pair12_typed_carrier_bridge": current_actual_selector_witness_codomain_still_lacks_pair12_typed_carrier_bridge,
            "current_qw2191_safe_resolution_still_lacks_pair12_typed_branch_provider": current_qw2191_safe_resolution_still_lacks_pair12_typed_branch_provider,
            "current_selector_witness_target_still_lacks_local_chart_sensitive_pair12_bridge": current_selector_witness_target_still_lacks_local_chart_sensitive_pair12_bridge,
            "next_honest_move_is_actual_t213_realization_attempt_or_further_provider_attack": next_honest_move_is_actual_t213_realization_attempt_or_further_provider_attack,
            "no_false_pass": True,
        },
        "checks": checks,
        "blocking_checks": blocking,
        "positive_actual_realization_candidates": positive_actual_realization_candidates,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P760",
        "status": status,
        "as_of": AS_OF,
        "t214_target_name": artifact["theorem_result"]["t214_target_name"],
        "t214_target_exported_on_current_repo_state": artifact["theorem_result"][
            "t214_target_exported_on_current_repo_state"
        ],
        "current_repo_still_does_not_export_actual_realization_of_t213_target": artifact["theorem_result"][
            "current_repo_still_does_not_export_actual_realization_of_t213_target"
        ],
        "current_actual_source_topology_selector_witness_still_not_pair12_typed": artifact[
            "theorem_result"
        ]["current_actual_source_topology_selector_witness_still_not_pair12_typed"],
        "current_actual_selector_witness_codomain_still_lacks_pair12_typed_carrier_bridge": artifact[
            "theorem_result"
        ]["current_actual_selector_witness_codomain_still_lacks_pair12_typed_carrier_bridge"],
        "current_qw2191_safe_resolution_still_lacks_pair12_typed_branch_provider": artifact[
            "theorem_result"
        ]["current_qw2191_safe_resolution_still_lacks_pair12_typed_branch_provider"],
        "current_selector_witness_target_still_lacks_local_chart_sensitive_pair12_bridge": artifact[
            "theorem_result"
        ]["current_selector_witness_target_still_lacks_local_chart_sensitive_pair12_bridge"],
        "next_honest_move_is_actual_t213_realization_attempt_or_further_provider_attack": artifact[
            "theorem_result"
        ]["next_honest_move_is_actual_t213_realization_attempt_or_further_provider_attack"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
