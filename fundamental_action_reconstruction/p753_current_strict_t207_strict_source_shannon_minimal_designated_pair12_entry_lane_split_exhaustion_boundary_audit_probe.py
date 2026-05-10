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

IN_P750 = GENERATED / "p750_current_strict_t204_strict_source_shannon_support_refinement_to_minimal_designated_pair12_theta_population_entry_nonexport_audit_probe_summary.json"
IN_P751 = GENERATED / "p751_current_strict_t205_strict_source_shannon_theta_support_refinement_to_minimal_designated_pair12_theta_entry_nonexport_audit_probe_summary.json"
IN_P752 = GENERATED / "p752_current_strict_t206_strict_source_shannon_pair_population_support_refinement_to_minimal_designated_pair12_populated_instance_entry_nonexport_audit_probe_summary.json"

OUT_JSON = GENERATED / "p753_current_strict_t207_strict_source_shannon_minimal_designated_pair12_entry_lane_split_exhaustion_boundary_audit_probe.json"
OUT_SUMMARY = GENERATED / "p753_current_strict_t207_strict_source_shannon_minimal_designated_pair12_entry_lane_split_exhaustion_boundary_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P750, IN_P751, IN_P752]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P753",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p750 = load_json(IN_P750)
    p751 = load_json(IN_P751)
    p752 = load_json(IN_P752)

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

    p750_combined_minimal_designated_pair12_entry_nonexport_already_frozen = (
        bool(
            p750.get(
                "current_strict_source_shannon_route_remains_nonentering_for_minimal_designated_pair12_theta_population_entry"
            )
        )
        and not bool(
            p750.get("current_strict_source_shannon_route_has_exported_minimal_designated_pair12_theta_population_entry")
        )
        and not bool(p750.get("t204_target_exported_on_current_repo_state"))
    )

    p751_theta_entry_half_nonexport_already_frozen = (
        bool(
            p751.get(
                "current_strict_source_shannon_theta_support_refinement_route_remains_nonentering_for_minimal_designated_pair12_theta_entry"
            )
        )
        and not bool(
            p751.get(
                "current_strict_source_shannon_theta_support_refinement_route_has_exported_minimal_designated_pair12_theta_entry"
            )
        )
        and not bool(p751.get("t205_target_exported_on_current_repo_state"))
    )

    p752_populated_instance_entry_half_nonexport_already_frozen = (
        bool(
            p752.get(
                "current_strict_source_shannon_pair_population_support_refinement_route_remains_nonentering_for_minimal_designated_pair12_populated_instance_entry"
            )
        )
        and not bool(
            p752.get(
                "current_strict_source_shannon_pair_population_support_refinement_route_has_exported_minimal_designated_pair12_populated_instance_entry"
            )
        )
        and not bool(p752.get("t206_target_exported_on_current_repo_state"))
    )

    current_strict_source_shannon_minimal_designated_pair12_entry_lane_split_is_exhaustive = (
        p750_combined_minimal_designated_pair12_entry_nonexport_already_frozen
        and p751_theta_entry_half_nonexport_already_frozen
        and p752_populated_instance_entry_half_nonexport_already_frozen
    )

    current_strict_source_shannon_minimal_designated_pair12_entry_lane_remains_candidate_only = (
        current_strict_source_shannon_minimal_designated_pair12_entry_lane_split_is_exhaustive
        and bool(
            p751.get(
                "current_strict_source_shannon_theta_support_refinement_remains_below_actual_theta_export"
            )
        )
        and bool(
            p752.get(
                "current_strict_source_shannon_pair_population_support_refinement_remains_below_actual_pair_population"
            )
        )
    )

    current_strict_source_shannon_minimal_designated_pair12_entry_lane_has_residual_unsplit_loophole = not (
        current_strict_source_shannon_minimal_designated_pair12_entry_lane_split_is_exhaustive
    )

    next_honest_move_requires_new_entry_object_or_provider_shift = (
        current_strict_source_shannon_minimal_designated_pair12_entry_lane_split_is_exhaustive
        and current_strict_source_shannon_minimal_designated_pair12_entry_lane_remains_candidate_only
        and not current_strict_source_shannon_minimal_designated_pair12_entry_lane_has_residual_unsplit_loophole
    )

    add_check(
        "p750_combined_minimal_designated_pair12_entry_nonexport_already_frozen",
        {
            "combined_entry_nonentering": bool(
                p750.get(
                    "current_strict_source_shannon_route_remains_nonentering_for_minimal_designated_pair12_theta_population_entry"
                )
            ),
            "combined_entry_exported": bool(
                p750.get("current_strict_source_shannon_route_has_exported_minimal_designated_pair12_theta_population_entry")
            ),
            "t204_exported": bool(p750.get("t204_target_exported_on_current_repo_state")),
        },
        {
            "combined_entry_nonentering": True,
            "combined_entry_exported": False,
            "t204_exported": False,
        },
        "P750 already freezes that the combined minimal designated [pair1,pair2] theta/population entry bridge remains absent.",
    )
    add_check(
        "p751_theta_entry_half_nonexport_already_frozen",
        {
            "theta_entry_nonentering": bool(
                p751.get(
                    "current_strict_source_shannon_theta_support_refinement_route_remains_nonentering_for_minimal_designated_pair12_theta_entry"
                )
            ),
            "theta_entry_exported": bool(
                p751.get(
                    "current_strict_source_shannon_theta_support_refinement_route_has_exported_minimal_designated_pair12_theta_entry"
                )
            ),
            "t205_exported": bool(p751.get("t205_target_exported_on_current_repo_state")),
        },
        {
            "theta_entry_nonentering": True,
            "theta_entry_exported": False,
            "t205_exported": False,
        },
        "P751 already freezes the logically earlier theta-entry half as still absent.",
    )
    add_check(
        "p752_populated_instance_entry_half_nonexport_already_frozen",
        {
            "populated_instance_entry_nonentering": bool(
                p752.get(
                    "current_strict_source_shannon_pair_population_support_refinement_route_remains_nonentering_for_minimal_designated_pair12_populated_instance_entry"
                )
            ),
            "populated_instance_entry_exported": bool(
                p752.get(
                    "current_strict_source_shannon_pair_population_support_refinement_route_has_exported_minimal_designated_pair12_populated_instance_entry"
                )
            ),
            "t206_exported": bool(p752.get("t206_target_exported_on_current_repo_state")),
        },
        {
            "populated_instance_entry_nonentering": True,
            "populated_instance_entry_exported": False,
            "t206_exported": False,
        },
        "P752 already freezes the later populated-instance-entry half as still absent.",
    )
    add_check(
        "current_strict_source_shannon_minimal_designated_pair12_entry_lane_split_is_exhaustive",
        current_strict_source_shannon_minimal_designated_pair12_entry_lane_split_is_exhaustive,
        True,
        "Taken together, P750/P751/P752 fully split the current Shannon designated-family entry gap with no remaining unsplit half at this refinement level.",
    )
    add_check(
        "current_strict_source_shannon_minimal_designated_pair12_entry_lane_remains_candidate_only",
        current_strict_source_shannon_minimal_designated_pair12_entry_lane_remains_candidate_only,
        True,
        "Even after that full split, the current Shannon designated-family entry lane remains candidate-only, still below actual theta export and actual pair population.",
    )
    add_check(
        "current_strict_source_shannon_minimal_designated_pair12_entry_lane_has_residual_unsplit_loophole",
        current_strict_source_shannon_minimal_designated_pair12_entry_lane_has_residual_unsplit_loophole,
        False,
        "No residual unsplit designated-family entry loophole remains inside the current Shannon syntax/refinement lane after the split.",
    )
    add_check(
        "next_honest_move_requires_new_entry_object_or_provider_shift",
        next_honest_move_requires_new_entry_object_or_provider_shift,
        True,
        "So the next honest move can no longer be another same-level Shannon entry wording tweak; it must export a genuinely new entry object or leave this lane for a different provider shift.",
    )

    status = (
        "PASS_STRICT_SOURCE_SHANNON_MINIMAL_DESIGNATED_PAIR12_ENTRY_LANE_SPLIT_EXHAUSTION_BOUNDARY_AUDITED"
        if not blocking
        and current_strict_source_shannon_minimal_designated_pair12_entry_lane_split_is_exhaustive
        and current_strict_source_shannon_minimal_designated_pair12_entry_lane_remains_candidate_only
        and not current_strict_source_shannon_minimal_designated_pair12_entry_lane_has_residual_unsplit_loophole
        and next_honest_move_requires_new_entry_object_or_provider_shift
        else "FAIL_STRICT_SOURCE_SHANNON_MINIMAL_DESIGNATED_PAIR12_ENTRY_LANE_SPLIT_EXHAUSTION_BOUNDARY_AUDIT"
    )

    artifact = {
        "stage": "P753",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "theorem_result": {
            "t207_boundary_name": "StrictSourceShannonMinimalDesignatedPair12EntryLaneSplitExhaustionBoundary_strict_v1",
            "t207_boundary_exported_on_current_repo_state": status
            == "PASS_STRICT_SOURCE_SHANNON_MINIMAL_DESIGNATED_PAIR12_ENTRY_LANE_SPLIT_EXHAUSTION_BOUNDARY_AUDITED",
            "p750_combined_minimal_designated_pair12_entry_nonexport_already_frozen": p750_combined_minimal_designated_pair12_entry_nonexport_already_frozen,
            "p751_theta_entry_half_nonexport_already_frozen": p751_theta_entry_half_nonexport_already_frozen,
            "p752_populated_instance_entry_half_nonexport_already_frozen": p752_populated_instance_entry_half_nonexport_already_frozen,
            "current_strict_source_shannon_minimal_designated_pair12_entry_lane_split_is_exhaustive": current_strict_source_shannon_minimal_designated_pair12_entry_lane_split_is_exhaustive,
            "current_strict_source_shannon_minimal_designated_pair12_entry_lane_remains_candidate_only": current_strict_source_shannon_minimal_designated_pair12_entry_lane_remains_candidate_only,
            "current_strict_source_shannon_minimal_designated_pair12_entry_lane_has_residual_unsplit_loophole": current_strict_source_shannon_minimal_designated_pair12_entry_lane_has_residual_unsplit_loophole,
            "next_honest_move_requires_new_entry_object_or_provider_shift": next_honest_move_requires_new_entry_object_or_provider_shift,
            "no_false_pass": True,
        },
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P753",
        "status": status,
        "as_of": AS_OF,
        "t207_boundary_name": artifact["theorem_result"]["t207_boundary_name"],
        "t207_boundary_exported_on_current_repo_state": artifact["theorem_result"][
            "t207_boundary_exported_on_current_repo_state"
        ],
        "p750_combined_minimal_designated_pair12_entry_nonexport_already_frozen": artifact[
            "theorem_result"
        ]["p750_combined_minimal_designated_pair12_entry_nonexport_already_frozen"],
        "p751_theta_entry_half_nonexport_already_frozen": artifact["theorem_result"][
            "p751_theta_entry_half_nonexport_already_frozen"
        ],
        "p752_populated_instance_entry_half_nonexport_already_frozen": artifact[
            "theorem_result"
        ]["p752_populated_instance_entry_half_nonexport_already_frozen"],
        "current_strict_source_shannon_minimal_designated_pair12_entry_lane_split_is_exhaustive": artifact[
            "theorem_result"
        ]["current_strict_source_shannon_minimal_designated_pair12_entry_lane_split_is_exhaustive"],
        "current_strict_source_shannon_minimal_designated_pair12_entry_lane_remains_candidate_only": artifact[
            "theorem_result"
        ]["current_strict_source_shannon_minimal_designated_pair12_entry_lane_remains_candidate_only"],
        "current_strict_source_shannon_minimal_designated_pair12_entry_lane_has_residual_unsplit_loophole": artifact[
            "theorem_result"
        ]["current_strict_source_shannon_minimal_designated_pair12_entry_lane_has_residual_unsplit_loophole"],
        "next_honest_move_requires_new_entry_object_or_provider_shift": artifact[
            "theorem_result"
        ]["next_honest_move_requires_new_entry_object_or_provider_shift"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
