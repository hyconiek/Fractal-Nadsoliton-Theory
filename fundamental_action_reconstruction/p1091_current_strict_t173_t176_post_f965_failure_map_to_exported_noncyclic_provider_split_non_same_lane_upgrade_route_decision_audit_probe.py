#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

AS_OF = "2026-03-28"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1090 = GENERATED / "p1090_current_strict_t173_t176_nadsoliton_information_process_to_orientation_supplier_family_failure_map_audit_probe_summary.json"
IN_F965 = GENERATED / "f965_current_strict_t173_t176_nadsoliton_information_process_to_orientation_supplier_failure_map_packet_summary.json"
IN_P982 = GENERATED / "p982_current_strict_qw2191_pair12_entry_point_same_lane_exhaustion_and_noncyclic_pivot_audit_probe_summary.json"
IN_F949 = GENERATED / "f949_first_current_strict_qw2191_pair12_entry_point_same_lane_exhaustion_and_noncyclic_pivot_packet_summary.json"
IN_P1064 = GENERATED / "p1064_current_nad12_sigma_pair_side_sharper_same_lane_witness_refinement_exhaustion_and_feeder_pivot_audit_probe_summary.json"
IN_F963 = GENERATED / "f963_current_nad12_sigma_pair_side_sharper_same_lane_witness_refinement_exhaustion_and_feeder_pivot_packet_summary.json"
IN_P1089 = GENERATED / "p1089_current_canonical_ontology_supported_nad12_sigma_feeder_sharper_same_lane_witness_refinement_stagnation_and_stop_audit_probe_summary.json"
IN_F964 = GENERATED / "f964_current_canonical_ontology_supported_nad12_sigma_feeder_sharper_same_lane_witness_refinement_stagnation_stop_packet_summary.json"

OUT_JSON = GENERATED / "p1091_current_strict_t173_t176_post_f965_failure_map_to_exported_noncyclic_provider_split_non_same_lane_upgrade_route_decision_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1091_current_strict_t173_t176_post_f965_failure_map_to_exported_noncyclic_provider_split_non_same_lane_upgrade_route_decision_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P1090, IN_F965, IN_P982, IN_F949, IN_P1064, IN_F963, IN_P1089, IN_F964]
    missing = [rel(p) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P1091",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1090 = load_json(IN_P1090)
    f965 = load_json(IN_F965)
    p982 = load_json(IN_P982)
    f949 = load_json(IN_F949)
    p1064 = load_json(IN_P1064)
    f963 = load_json(IN_F963)
    p1089 = load_json(IN_P1089)
    f964 = load_json(IN_F964)

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

    add_check(
        "p1090_f965_failure_map_is_real_but_not_supplier",
        (
            bool(p1090.get("process_reading_of_nadsoliton_is_repo_real"))
            and bool(p1090.get("current_repo_exports_lawful_process_to_orientation_supplier")) is False
            and f965.get("failure_map_grade") == "frontier_failure_map_only"
            and bool(f965.get("counts_as_lawful_supplier")) is False
        ),
        True,
        "P1090/F965 already freeze the informational-process reading as repo-real while keeping lawful supplier export absent.",
    )
    add_check(
        "p982_f949_pair12_entry_same_lane_reentry_is_disallowed",
        (
            bool(p982.get("same_lane_exhaustion_boundary_reached"))
            and bool(p982.get("further_same_lane_t274_style_descent_is_not_honest_primary_move"))
            and bool(p982.get("next_honest_move_is_noncyclic_pivot_to_exported_provider_split_family"))
            and p982.get("preferred_noncyclic_pivot_family") == "nad12_sigma_residual_shannon_noncyclic_provider_split_family"
            and bool(f949.get("same_lane_exhaustion_boundary_reached"))
            and f949.get("preferred_noncyclic_pivot_family") == "Xi_nad12_sigma_residual_shannon_noncyclic_provider_split_target_v1"
        ),
        True,
        "P982/F949 already disallow reopening the old pair12 entry same-lane descent as a primary move.",
    )
    add_check(
        "p1064_f963_pair_side_same_lane_reentry_is_disallowed",
        (
            bool(p1064.get("same_lane_exhaustion_boundary_reached"))
            and bool(p1064.get("further_same_lane_sharper_witness_descent_is_not_honest_primary_move"))
            and bool(p1064.get("next_honest_move_is_noncyclic_pivot_to_feeder_support_side_witness_family"))
            and p1064.get("preferred_noncyclic_pivot_family") == "Xi_nad12_sigma_residual_shannon_noncyclic_provider_split_target_v1"
            and bool(f963.get("same_lane_exhaustion_boundary_reached"))
            and f963.get("preferred_noncyclic_pivot_family") == "Xi_nad12_sigma_residual_shannon_noncyclic_provider_split_target_v1"
        ),
        True,
        "P1064/F963 already disallow reopening the downstream pair-side sharper same-lane witness family as a primary move.",
    )
    add_check(
        "p1089_f964_feeder_same_lane_reentry_is_disallowed",
        (
            bool(p1089.get("same_lane_stagnation_boundary_reached"))
            and bool(p1089.get("further_same_lane_sharper_witness_descent_is_not_honest_primary_move"))
            and bool(p1089.get("restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route"))
            and bool(f964.get("same_lane_stagnation_boundary_reached"))
            and bool(f964.get("same_lane_deeper_witness_descent_disallowed_as_primary_move"))
            and bool(f964.get("restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route"))
        ),
        True,
        "P1089/F964 already disallow reopening the downstream feeder sharper same-lane family as a primary move.",
    )
    add_check(
        "same_exported_noncyclic_provider_split_family_remains_the_honest_search_space",
        (
            f949.get("preferred_noncyclic_pivot_family") == "Xi_nad12_sigma_residual_shannon_noncyclic_provider_split_target_v1"
            and f963.get("preferred_noncyclic_pivot_family") == "Xi_nad12_sigma_residual_shannon_noncyclic_provider_split_target_v1"
            and bool(p1089.get("restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route"))
        ),
        True,
        "The repo keeps one exported noncyclic provider-split family as the lawful live search substrate, but only for non-same-lane upgrades.",
    )

    discharged = len(blocking) == 0
    status = (
        "PASS_CURRENT_STRICT_T173_T176_POST_F965_FAILURE_MAP_TO_EXPORTED_NONCYCLIC_PROVIDER_SPLIT_NON_SAME_LANE_UPGRADE_ROUTE_DECISION_AUDITED"
        if discharged
        else "P1091_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_ROUTE_STATE"
    )

    route_decision = {
        "process_failure_map_frozen": bool(p1090.get("process_reading_of_nadsoliton_is_repo_real")),
        "pair12_entry_same_lane_reentry_disallowed_as_primary_move": bool(p982.get("further_same_lane_t274_style_descent_is_not_honest_primary_move")),
        "pair_side_sharper_same_lane_reentry_disallowed_as_primary_move": bool(p1064.get("further_same_lane_sharper_witness_descent_is_not_honest_primary_move")),
        "feeder_sharper_same_lane_reentry_disallowed_as_primary_move": bool(f964.get("same_lane_deeper_witness_descent_disallowed_as_primary_move")),
        "preferred_search_family": "Xi_nad12_sigma_residual_shannon_noncyclic_provider_split_target_v1",
        "allowed_next_move_contract": "search_one_genuinely_new_non_same_lane_upgrade_route_within_exported_noncyclic_provider_split_family_or_one_genuinely_new_inversion_sensitive_source_side_provider_class",
        "active_missing_bridge": p1090.get("active_missing_bridge"),
        "reopen_old_same_lane_families": False,
    }

    artifact = {
        "stage": "P1091",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "P1090": rel(IN_P1090),
            "F965": rel(IN_F965),
            "P982": rel(IN_P982),
            "F949": rel(IN_F949),
            "P1064": rel(IN_P1064),
            "F963": rel(IN_F963),
            "P1089": rel(IN_P1089),
            "F964": rel(IN_F964),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "route_decision": route_decision,
        "audit_conclusion": {
            "reopen_old_pair12_entry_same_lane": False,
            "reopen_old_pair_side_sharper_same_lane": False,
            "reopen_old_feeder_sharper_same_lane": False,
            "search_inside_exported_noncyclic_provider_split_family": discharged,
            "search_mode": route_decision["allowed_next_move_contract"],
        },
        "hard_limits": [
            "No lawful supplier export.",
            "No strict physical orientation datum export.",
            "No T183 discharge.",
            "No T176 discharge.",
            "No QW-2191 discharge.",
            "No ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "preferred_search_family": route_decision["preferred_search_family"],
        "pair12_entry_same_lane_reentry_disallowed_as_primary_move": route_decision["pair12_entry_same_lane_reentry_disallowed_as_primary_move"],
        "pair_side_sharper_same_lane_reentry_disallowed_as_primary_move": route_decision["pair_side_sharper_same_lane_reentry_disallowed_as_primary_move"],
        "feeder_sharper_same_lane_reentry_disallowed_as_primary_move": route_decision["feeder_sharper_same_lane_reentry_disallowed_as_primary_move"],
        "allowed_next_move_contract": route_decision["allowed_next_move_contract"],
        "active_missing_bridge": route_decision["active_missing_bridge"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
