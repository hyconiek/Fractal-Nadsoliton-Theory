#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "f5_legacy_physical_role_partition_refinement_packet.json"
OUT_SUMMARY = (
    GENERATED / "f5_legacy_physical_role_partition_refinement_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def load_text(repo_relative_path: str) -> str:
    return (REPO / repo_relative_path).read_text(encoding="utf-8")


def contains_all(text: str, parts: list[str]) -> bool:
    return all(part in text for part in parts)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p62 = load_json(
        "fundamental_action_reconstruction/generated/p62_legacy_physical_role_transfer_to_strict_gate_kernel_probe_summary.json"
    )
    qw2005 = load_text(
        "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2005_PRE1700_TEX_REVALIDATION_MATRIX_EN_PL.md"
    )

    broad_retained_partition_present = contains_all(
        qw2005,
        [
            "Retained:",
            "Nadsoliton conceptual ontology",
            "oscillatory+damped kernel structure idea",
            "laminar-vacuum interpretation over turbulence",
        ],
    )
    broad_replaced_partition_present = contains_all(
        qw2005,
        [
            "Replaced / upgraded:",
            "strict closure stack",
            "GW robustness methodology",
            "lockability criterion and deep-audit requirements",
            "frozen package formalization",
        ],
    )
    claim_specific_strict_partition_present = False

    checks_spec = [
        {
            "id": "p62_non_transfer_probe_present",
            "actual": p62["status"],
            "expected": "CURRENT_REPO_EXPORTS_LEGACY_PHYSICAL_ROLE_CLAIMS_AND_STRICT_KERNEL_PIPELINE_BUT_NO_RIGOROUS_LEGACY_ROLE_TRANSFER_TO_STRICT_GATE_KERNEL_AFTER_P62",
            "meaning": "P62 already fixes the broad non-transfer frontier for the three legacy roles",
        },
        {
            "id": "q2005_broad_retained_partition_present",
            "actual": broad_retained_partition_present,
            "expected": True,
            "meaning": "QW-2005 exports a broad retained partition at conceptual level",
        },
        {
            "id": "q2005_broad_replaced_partition_present",
            "actual": broad_replaced_partition_present,
            "expected": True,
            "meaning": "QW-2005 exports a broad replaced/upgraded partition at methodology/closure level",
        },
        {
            "id": "claim_specific_strict_side_partition_present",
            "actual": claim_specific_strict_partition_present,
            "expected": False,
            "meaning": "the repo does not yet export a claim-specific strict-side retained-or-replaced verdict for C3/C4/C5",
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
        "stage": "F5",
        "lane": "legacy_physical_role_partition_refinement_current_repo_state_only",
        "goal": "refine_the_broad_missing_retained_vs_replaced_partition_theorem_into_claim_specific_partition_gaps",
        "status": "F5_EXECUTED_LEGACY_PHYSICAL_ROLE_PARTITION_REFINEMENT_PACKET_NO_FALSE_PASS",
        "reason": "QW-2005 already exports a broad retained-vs-replaced partition, but only at conceptual and methodology level; it does not yet export claim-specific strict-side retained-or-replaced verdicts for the legacy Weinberg-angle, fine-structure, or gravity-hierarchy roles",
        "broad_partition_state": {
            "retained": [
                "Nadsoliton conceptual ontology",
                "oscillatory+damped kernel structure idea",
                "laminar-vacuum interpretation over turbulence",
            ],
            "replaced_or_upgraded": [
                "strict closure stack",
                "GW robustness methodology",
                "lockability criterion and deep-audit requirements",
                "frozen package formalization",
            ],
        },
        "claim_specific_partition_state": {
            "legacy_weinberg_angle_role_partition_verdict_present": False,
            "legacy_fine_structure_role_partition_verdict_present": False,
            "legacy_gravity_hierarchy_role_partition_verdict_present": False,
        },
        "remaining_missing_objects": [
            "explicit_strict_side_retained_or_replaced_verdict_for_the_legacy_weinberg_angle_role",
            "explicit_strict_side_retained_or_replaced_verdict_for_the_legacy_fine_structure_role",
            "explicit_strict_side_retained_or_replaced_verdict_for_the_legacy_gravity_hierarchy_role",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F5",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "broad_partition_state": artifact["broad_partition_state"],
        "claim_specific_partition_state": artifact["claim_specific_partition_state"],
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
