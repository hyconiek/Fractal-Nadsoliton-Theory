#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-04-25"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_S3 = ROOT / "S3_CURRENT_INFORMATION_COMPRESSION_DECOMPRESSION_ROUTE_CLASSIFICATION_NOTE.md"
IN_F966 = ROOT / "F966_CURRENT_STRICT_T173_T176_POST_F965_FAILURE_MAP_TO_EXPORTED_NONCYCLIC_PROVIDER_SPLIT_NON_SAME_LANE_UPGRADE_ROUTE_PACKET.md"
IN_F969 = ROOT / "F969_CURRENT_STRICT_T173_T176_MINIMAL_ORIENTED_NONRECIPROCAL_DEPHASING_NEW_IMPORT_BOUNDARY_PACKET.md"
IN_F970 = ROOT / "F970_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_STAGNATION_STOP_PACKET.md"
IN_F971 = ROOT / "F971_CURRENT_STRICT_T173_T176_POST_F970_MINIMAL_ONRD_SSL_RCW_STOP_TO_EXISTING_F966_NON_SAME_LANE_ROUTE_PACKET.md"
IN_F972 = ROOT / "F972_CURRENT_STRICT_T173_T176_EXISTING_F966_NON_SAME_LANE_ROUTE_CONTRACT_TO_EXPORTED_STRICT_SELECTOR_SOURCE_FRONTIER_SUBROUTE_PACKET.md"
IN_F973 = ROOT / "F973_CURRENT_STRICT_T173_T176_EXISTING_F972_SELECTOR_SOURCE_SUBROUTE_TO_EXISTING_F960_ACTIVE_BRIDGE_FRONTIER_ROUTE_PACKET.md"
IN_F975 = ROOT / "F975_CURRENT_STRICT_T173_T176_EXISTING_F974_F960_EXISTING_T178_SOURCE_TO_ATLAS_CHART_SEED_SELECTION_BRIDGE_CURRENT_BEST_NEW_NARROW_PROVIDER_CLASS_SEED_CANDIDATE_PACKET.md"
IN_N965 = ROOT / "N965_CURRENT_STRICT_T173_T176_EXISTING_F975_F960_T178_CURRENT_BEST_ACTUAL_REALIZATION_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_YET_FURTHER_LOWER_BOUNDARY_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORT_THEOREM.md"

OUT_JSON = GENERATED / "p1127_current_strict_t173_t176_s3_f966_f972_f960_source_side_oriented_memory_rule_target_admission_probe.json"
OUT_SUMMARY = GENERATED / "p1127_current_strict_t173_t176_s3_f966_f972_f960_source_side_oriented_memory_rule_target_admission_probe_summary.json"

ACTIVE_BRIDGE = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"
F975_SEED = "SourceTopologyToAtlasChartSeedSelectionBridge_global_C_v1_strict_v1"
TARGET_NAME = "SourceSideOrientedMemoryRuleTarget_against_ResidualDatumPair12OrbitDirectionSelectionBridge_v1"


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def scan_existing_target_hits() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py", "generated/*.json")
    excluded = {
        Path(__file__).name,
        "P1127_CURRENT_STRICT_T173_T176_EXISTING_S3_F966_F972_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_TARGET_ADMISSION_PROBE.md",
        "T332_CURRENT_STRICT_T173_T176_EXISTING_S3_F966_F972_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_TARGET_SPEC.md",
        "N966_CURRENT_STRICT_T173_T176_EXISTING_S3_F966_F972_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_TARGET_ADMISSION_THEOREM.md",
        OUT_JSON.name,
        OUT_SUMMARY.name,
    }
    hits: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if path in seen or path.name in excluded:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            if TARGET_NAME in text:
                hits.append(rel(path))
    return hits


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_S3, IN_F966, IN_F969, IN_F970, IN_F971, IN_F972, IN_F973, IN_F975, IN_N965]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1127",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    s3_text = load_text(IN_S3)
    f966_text = load_text(IN_F966)
    f969_text = load_text(IN_F969)
    f970_text = load_text(IN_F970)
    f971_text = load_text(IN_F971)
    f972_text = load_text(IN_F972)
    f973_text = load_text(IN_F973)
    f975_text = load_text(IN_F975)
    n965_text = load_text(IN_N965)
    existing_hits = scan_existing_target_hits()

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

    s3_source_side_oriented_memory_reading_live = all(
        needle in s3_text
        for needle in [
            "Current classification:",
            "still live",
            "candidate provider-class seed",
            "freeze one exact source-side oriented-memory-rule target",
            "strict-side provenance only",
            "no observer input",
            "explicit nonreciprocal or hysteretic term",
            "future locality and positivity audit hooks",
        ]
    )

    live_f966_f972_f973_contract_confirmed = all(
        needle in f966_text + f972_text + f973_text
        for needle in [
            ACTIVE_BRIDGE,
            "one genuinely new non-same-lane upgrade route",
            "explicit_strict_internal_selector_source_derivation_frontier",
            "build_one_narrow_probe_for_a_genuinely_new_inversion_sensitive_source_side_provider_class_against_the_existing_f960_bridge_target",
        ]
    )

    f969_oriented_memory_rule_component_explicit = all(
        needle in f969_text
        for needle in [
            "orientation_anchor_component := provenance_safe_orientation_anchor",
            "nonreciprocal_rate_asymmetry_component := nonreciprocal_rate_asymmetry_rule",
            "oriented_memory_rule_component := oriented_memory_rule",
            "admissible_as_candidate_provider_class_seed_only := yes",
            ACTIVE_BRIDGE,
        ]
    )

    onrd_same_lane_restart_disallowed = all(
        needle in f970_text + f971_text
        for needle in [
            "same_lane_deeper_route_coherence_witness_descent_disallowed_as_primary_move := yes",
            "restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route := yes",
            "current_primary_work_contract := rejoin_existing_f966_non_same_lane_route_do_not_spawn_competing_post_f970_same_lane",
        ]
    )

    f975_current_best_t178_seed_still_nonrealized = all(
        needle in f975_text + n965_text
        for needle in [
            f"current_best_new_narrow_provider_class_seed_candidate := {F975_SEED}",
            "counts_as_actual_t178_export := no",
            "lawful verdict for T329",
            "actual T178 export",
            "T176 discharge",
            "QW-2191 discharge",
        ]
    )

    no_existing_oriented_memory_rule_target_export = len(existing_hits) == 0

    strongest_honest_next_move_is_freeze_additional_target = (
        s3_source_side_oriented_memory_reading_live
        and live_f966_f972_f973_contract_confirmed
        and f969_oriented_memory_rule_component_explicit
        and onrd_same_lane_restart_disallowed
        and f975_current_best_t178_seed_still_nonrealized
        and no_existing_oriented_memory_rule_target_export
    )

    add_check(
        "s3_source_side_oriented_memory_reading_live",
        s3_source_side_oriented_memory_reading_live,
        True,
        "S3 already classifies source-side oriented memory as the only live rigorous compression/decompression reading.",
    )
    add_check(
        "live_f966_f972_f973_contract_confirmed",
        live_f966_f972_f973_contract_confirmed,
        True,
        "F966/F972/F973 keep the live non-same-lane selector-source contract explicit at the active F960 bridge frontier.",
    )
    add_check(
        "f969_oriented_memory_rule_component_explicit",
        f969_oriented_memory_rule_component_explicit,
        True,
        "F969 already keeps oriented_memory_rule explicit as one irreducible component of the minimal ONRD boundary.",
    )
    add_check(
        "onrd_same_lane_restart_disallowed",
        onrd_same_lane_restart_disallowed,
        True,
        "F970/F971 already reject restarting the old same-lane ONRD refinement as the primary move.",
    )
    add_check(
        "f975_current_best_t178_seed_still_nonrealized",
        f975_current_best_t178_seed_still_nonrealized,
        True,
        "F975 and the latest N965 chain keep the T178 lane explicit but still below actual bridge export and strict discharge.",
    )
    add_check(
        "no_existing_oriented_memory_rule_target_export",
        no_existing_oriented_memory_rule_target_export,
        True,
        "No already-exported exact target yet isolates the source-side oriented-memory-rule target itself beneath the live contract.",
    )
    add_check(
        "strongest_honest_next_move_is_freeze_additional_target",
        strongest_honest_next_move_is_freeze_additional_target,
        True,
        "Therefore the strongest honest additional exact export is one source-side oriented-memory-rule target admitted only as candidate_provider_class_seed_only and not as a fiat replacement of F975.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_EXISTING_S3_F966_F972_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_TARGET_ADMISSION_AUDITED"
        if not blocking and strongest_honest_next_move_is_freeze_additional_target
        else "FAIL_CURRENT_STRICT_T173_T176_EXISTING_S3_F966_F972_F960_SOURCE_SIDE_ORIENTED_MEMORY_RULE_TARGET_ADMISSION_AUDIT"
    )

    artifact = {
        "stage": "P1127",
        "status": status,
        "as_of": AS_OF,
        "target_name": TARGET_NAME,
        "active_missing_bridge": ACTIVE_BRIDGE,
        "current_best_f975_seed_candidate": F975_SEED,
        "s3_source_side_oriented_memory_reading_live": s3_source_side_oriented_memory_reading_live,
        "live_f966_f972_f973_contract_confirmed": live_f966_f972_f973_contract_confirmed,
        "f969_oriented_memory_rule_component_explicit": f969_oriented_memory_rule_component_explicit,
        "onrd_same_lane_restart_disallowed": onrd_same_lane_restart_disallowed,
        "f975_current_best_t178_seed_still_nonrealized": f975_current_best_t178_seed_still_nonrealized,
        "no_existing_oriented_memory_rule_target_export": no_existing_oriented_memory_rule_target_export,
        "current_repo_already_exports_same_target_hits": existing_hits,
        "strongest_honest_next_move_is_freeze_additional_target": strongest_honest_next_move_is_freeze_additional_target,
        "admission_is_additional_not_f975_replacement": True,
        "counts_as_lawful_supplier": False,
        "counts_as_strict_physical_orientation_datum": False,
        "counts_as_t183_discharge": False,
        "counts_as_t176_discharge": False,
        "counts_as_qw2191_discharge": False,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "target_name": artifact["target_name"],
        "active_missing_bridge": artifact["active_missing_bridge"],
        "current_best_f975_seed_candidate": artifact["current_best_f975_seed_candidate"],
        "strongest_honest_next_move_is_freeze_additional_target": artifact[
            "strongest_honest_next_move_is_freeze_additional_target"
        ],
        "admission_is_additional_not_f975_replacement": artifact[
            "admission_is_additional_not_f975_replacement"
        ],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
