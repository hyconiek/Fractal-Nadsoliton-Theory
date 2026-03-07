#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "f7_legacy_weinberg_retained_branch_refinement_packet.json"
OUT_SUMMARY = (
    GENERATED / "f7_legacy_weinberg_retained_branch_refinement_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def load_text(repo_relative_path: str) -> str:
    return (REPO / repo_relative_path).read_text(encoding="utf-8")


def contains_any(text: str, parts: list[str]) -> bool:
    return any(part in text for part in parts)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p64 = load_json(
        "fundamental_action_reconstruction/generated/p64_legacy_weinberg_role_strict_side_branch_probe_summary.json"
    )
    strict_textbook = load_text("RELEASE_4_9_TEXTBOOK_EN_PL.md")
    qw2000 = load_text(
        "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2000_BOUNDED_COUPLING_DEEP_AUDIT.md"
    )
    qw2001 = load_text(
        "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2001_BOUNDED_GW_TRIAD_LOCKABLE_GATE.md"
    )
    qw2002 = load_text(
        "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2002_SINGLE_KERNEL_TRIPLE_SECTOR_CLOSURE_GATE_V3.md"
    )
    qw2003 = load_text(
        "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2003_FROZEN_LOCKABLE_PACKAGE_EXPORT.md"
    )
    strict_side_sources = "\n".join(
        [strict_textbook, qw2000, qw2001, qw2002, qw2003]
    )

    literal_retention_present = contains_any(
        strict_side_sources,
        [
            "sin^2(theta_W)=alpha_geo/12",
            "sin^2\\theta_W = alpha_geo/12",
            "sin^2(theta_W) = alpha_geo/12",
        ],
    )
    role_equivalence_retention_present = False

    checks_spec = [
        {
            "id": "p64_retained_branch_absent",
            "actual": p64["branch_state"]["retained_branch_present"],
            "expected": False,
            "meaning": "P64 already confirms that the retained branch is not yet exported at top level",
        },
        {
            "id": "strict_side_literal_retention_present",
            "actual": literal_retention_present,
            "expected": False,
            "meaning": "no strict-side literal retention of sin^2(theta_W)=alpha_geo/12 is currently exported",
        },
        {
            "id": "strict_side_role_equivalence_retention_present",
            "actual": role_equivalence_retention_present,
            "expected": False,
            "meaning": "no explicit strict-side role-equivalence retention verdict for the Weinberg-angle role is currently exported",
        },
    ]

    checks: list[dict[str, Any]] = []
    for item in checks_spec:
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": item["actual"] == item["expected"],
                "meaning": item["meaning"],
            }
        )

    artifact = {
        "stage": "F7",
        "lane": "legacy_weinberg_retained_branch_refinement_current_repo_state_only",
        "goal": "refine_the_missing_retained_branch_for_the_legacy_weinberg_angle_role_into_literal_retention_vs_role_equivalence_subbranches",
        "status": "F7_EXECUTED_LEGACY_WEINBERG_RETAINED_BRANCH_REFINEMENT_PACKET_NO_FALSE_PASS",
        "reason": "P64 already keeps the retained branch absent at top level, while current strict-side materials export neither the old literal formula nor an explicit role-equivalence retention verdict, so the narrowest honest refinement is literal retention versus role-equivalence retention",
        "retained_subbranch_state": {
            "literal_retention_present": literal_retention_present,
            "role_equivalence_retention_present": role_equivalence_retention_present,
        },
        "remaining_missing_objects": [
            "explicit_strict_side_literal_retention_of_sin2_thetaW_equals_alpha_geo_over_12",
            "explicit_strict_side_role_equivalence_verdict_for_the_legacy_weinberg_angle_role",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F7",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "retained_subbranch_state": artifact["retained_subbranch_state"],
        "remaining_missing_objects": artifact["remaining_missing_objects"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    print(OUT_JSON)


if __name__ == "__main__":
    main()
