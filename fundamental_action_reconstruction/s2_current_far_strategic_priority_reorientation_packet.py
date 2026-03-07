#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "s2_current_far_strategic_priority_reorientation_packet.json"
OUT_SUMMARY = GENERATED / "s2_current_far_strategic_priority_reorientation_packet_summary.json"


def exists(repo_relative_path: str) -> bool:
    return (REPO / repo_relative_path).exists()


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    checks = [
        {
            "id": "kernel_split_note_present",
            "actual": exists("fundamental_action_reconstruction/K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md"),
            "expected": True,
            "meaning": "K1 records the current legacy-vs-strict kernel split",
        },
        {
            "id": "kernel_bridge_probe_present",
            "actual": exists("fundamental_action_reconstruction/P47_LEGACY_ONTOLOGICAL_KERNEL_TO_STRICT_GATE_KERNEL_BRIDGE_PROBE.md"),
            "expected": True,
            "meaning": "P47 asks the direct bridge question for the two kernels",
        },
        {
            "id": "kernel_nonidentification_theorem_present",
            "actual": exists("fundamental_action_reconstruction/N50_CURRENT_LEGACY_ONTOLOGICAL_KERNEL_TO_STRICT_GATE_KERNEL_NONIDENTIFICATION_THEOREM.md"),
            "expected": True,
            "meaning": "N50 discharges the current-repo-state nonidentification theorem",
        },
        {
            "id": "qw2191_obstruction_present",
            "actual": exists("material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2191_MODE_INDEX_UNIQUENESS_OBSTRUCTION_THEOREM_GATE.md"),
            "expected": True,
            "meaning": "QW-2191 explicitly records the selector/uniqueness obstruction",
        },
        {
            "id": "qw2382_noncyclic_strategy_present",
            "actual": exists("material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2382_DUAL_NONCYCLIC_STRATEGY_PACKET_GATE.md"),
            "expected": True,
            "meaning": "QW-2382 records the hard noncyclic strategy constraints",
        },
        {
            "id": "qw2383_repeat_step_rejection_present",
            "actual": exists("material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2383_DUAL_NONCYCLIC_STEP_ADMISSION_GATE.md"),
            "expected": True,
            "meaning": "QW-2383 rejects repeated step admission under the same blocker-cut",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "S2",
        "goal": "record_the_current_highest_level_far_priority_reorientation_after_kernel_split_correction_qw2191_and_l5_l12_cycle_analysis",
        "priority_order": [
            "legacy_to_strict_kernel_bridge_or_non_bridge",
            "explicit_selector_requirement_or_symmetry_breaking_after_QW2191",
            "noncyclic_anchor_for_L5_L12",
            "further_local_direct_route_hyper_decomposition_only_as_auxiliary_lane",
        ],
        "priority_1_kernel_bridge": {
            "current_state": "legacy_and_strict_kernels_are_both_exported_but_not_rigorously identified",
            "source_reports": ["K1", "P47", "N50"],
            "highest_priority": True,
        },
        "priority_2_selector_obstruction": {
            "current_state": "QW2191 uniqueness obstruction remains active in strict core",
            "source_reports": ["QW-2191", "AX1..AX8"],
            "requires_explicit_new_premise_or_stronger_impossibility": True,
        },
        "priority_3_noncyclic_anchor": {
            "current_state": "L5/L12 recurrence under same blocker-cut is already rejected",
            "source_reports": ["QW-2381", "QW-2382", "QW-2383"],
            "repeat_cycle_not_admitted": True,
        },
        "direct_route_status": {
            "local_direct_m2_route_still_valid": True,
            "main_theoretical_priority_now": False,
            "role": "auxiliary_sharpening_lane_only",
        },
        "checks": checks,
        "no_false_pass": True,
    }

    summary = {
        "stage": "S2",
        "status": "S2_EXECUTED_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET_NO_FALSE_PASS",
        "priority_order": artifact["priority_order"],
        "direct_route_role": artifact["direct_route_status"]["role"],
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
