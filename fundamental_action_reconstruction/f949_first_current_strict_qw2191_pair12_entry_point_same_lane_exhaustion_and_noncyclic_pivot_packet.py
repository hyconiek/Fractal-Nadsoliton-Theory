#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P982 = GENERATED / "p982_current_strict_qw2191_pair12_entry_point_same_lane_exhaustion_and_noncyclic_pivot_audit_probe_summary.json"
IN_N355 = ROOT / "N355_CURRENT_FIRST_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_NONCYCLIC_PROVIDER_SPLIT_TARGET_THEOREM.md"
IN_N362 = ROOT / "N362_CURRENT_FIRST_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_PAIR_REALIZATION_SIDE_PROVIDER_SUPPORT_WITNESS_THEOREM.md"
IN_F949 = ROOT / "F949_FIRST_CURRENT_STRICT_QW2191_PAIR12_ENTRY_POINT_SAME_LANE_EXHAUSTION_AND_NONCYCLIC_PIVOT_PACKET.md"

OUT_JSON = GENERATED / "f949_first_current_strict_qw2191_pair12_entry_point_same_lane_exhaustion_and_noncyclic_pivot_packet.json"
OUT_SUMMARY = GENERATED / "f949_first_current_strict_qw2191_pair12_entry_point_same_lane_exhaustion_and_noncyclic_pivot_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P982, IN_N355, IN_N362, IN_F949]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "F949",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p982 = load_json(IN_P982)
    n355_text = load_text(IN_N355)
    n362_text = load_text(IN_N362)
    f949_text = load_text(IN_F949)

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

    p982_same_lane_exhaustion_boundary_already_reached = (
        p982.get("status")
        == "PASS_STRICT_QW2191_PAIR12_ENTRY_POINT_SAME_LANE_EXHAUSTION_AND_NONCYCLIC_PIVOT_AUDITED"
        and bool(p982.get("same_lane_exhaustion_boundary_reached"))
        and bool(p982.get("further_same_lane_t274_style_descent_is_not_honest_primary_move"))
        and bool(p982.get("next_honest_move_is_noncyclic_pivot_to_exported_provider_split_family"))
    )

    exported_noncyclic_provider_split_family_is_named = (
        "Xi_nad12_sigma_residual_shannon_noncyclic_provider_split_target_v1" in n355_text
    )

    preferred_first_pivot_branch_is_named = (
        "Lambda_nad12_sigma_residual_shannon_pair_realization_side_provider_support_witness_v1" in n362_text
    )

    f949_packet_shape_frozen = all(
        needle in f949_text
        for needle in [
            "Xi_strict_qw2191_pair12_entry_point_same_lane_exhaustion_and_noncyclic_pivot_packet_v1",
            "qw2191_pair12_entry_point_same_lane_exhaustion_boundary_reached := yes",
            "same_lane_t274_style_descent_disallowed_as_primary_move := yes",
            "preferred_noncyclic_pivot_family := Xi_nad12_sigma_residual_shannon_noncyclic_provider_split_target_v1",
            "preferred_first_pivot_branch := Lambda_nad12_sigma_residual_shannon_pair_realization_side_provider_support_witness_v1",
            "current_primary_work_contract := pivot_to_exported_noncyclic_provider_split_family_not_fake_local_pass",
        ]
    )

    packet_exported_on_current_repo_state = (
        p982_same_lane_exhaustion_boundary_already_reached
        and exported_noncyclic_provider_split_family_is_named
        and preferred_first_pivot_branch_is_named
        and f949_packet_shape_frozen
    )

    add_check(
        "p982_same_lane_exhaustion_boundary_already_reached",
        p982_same_lane_exhaustion_boundary_already_reached,
        True,
        "P982 already freezes the same-lane exhaustion boundary and the noncyclic pivot requirement.",
    )
    add_check(
        "exported_noncyclic_provider_split_family_is_named",
        exported_noncyclic_provider_split_family_is_named,
        True,
        "The noncyclic provider-split family is explicitly named in the repo.",
    )
    add_check(
        "preferred_first_pivot_branch_is_named",
        preferred_first_pivot_branch_is_named,
        True,
        "The preferred first pivot branch is explicitly named in the repo.",
    )
    add_check(
        "f949_packet_shape_frozen",
        f949_packet_shape_frozen,
        True,
        "F949 freezes the noncyclic pivot packet shape explicitly.",
    )
    add_check(
        "packet_exported_on_current_repo_state",
        packet_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one honest noncyclic pivot packet for the exhausted QW-2191 same-lane descent.",
    )

    status = (
        "PASS_CURRENT_STRICT_QW2191_PAIR12_ENTRY_POINT_SAME_LANE_EXHAUSTION_AND_NONCYCLIC_PIVOT_PACKET_EXPORTED"
        if not blocking and packet_exported_on_current_repo_state
        else "FAIL_CURRENT_STRICT_QW2191_PAIR12_ENTRY_POINT_SAME_LANE_EXHAUSTION_AND_NONCYCLIC_PIVOT_PACKET"
    )

    artifact = {
        "stage": "F949",
        "status": status,
        "as_of": AS_OF,
        "packet_exported_on_current_repo_state": packet_exported_on_current_repo_state,
        "same_lane_exhaustion_boundary_reached": bool(p982.get("same_lane_exhaustion_boundary_reached")),
        "preferred_noncyclic_pivot_family": "Xi_nad12_sigma_residual_shannon_noncyclic_provider_split_target_v1",
        "preferred_first_pivot_branch": "Lambda_nad12_sigma_residual_shannon_pair_realization_side_provider_support_witness_v1",
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "packet_exported_on_current_repo_state": artifact["packet_exported_on_current_repo_state"],
        "same_lane_exhaustion_boundary_reached": artifact["same_lane_exhaustion_boundary_reached"],
        "preferred_noncyclic_pivot_family": artifact["preferred_noncyclic_pivot_family"],
        "preferred_first_pivot_branch": artifact["preferred_first_pivot_branch"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
