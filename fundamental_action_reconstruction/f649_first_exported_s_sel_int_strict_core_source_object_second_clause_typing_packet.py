#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

OUT_JSON = (
    GENERATED
    / "f649_first_exported_s_sel_int_strict_core_source_object_second_clause_typing_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f649_first_exported_s_sel_int_strict_core_source_object_second_clause_typing_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n540 = load_json(
        "fundamental_action_reconstruction/generated/n540_current_first_admissibility_clause_discharge_theorem_for_s_sel_int_strict_core_source_object_v1_summary.json"
    )
    f647 = load_json(
        "fundamental_action_reconstruction/generated/f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt.json"
    )

    exported_object = str(n540["theorem_result"]["exported_object"] or "S_sel_int_strict_core_source_object_v1")

    dot = float(
        f647.get("constructed_source_object", {})
        .get("exported_payload", {})
        .get("dot_w_break_with_s1", 0.0)
    )

    two_pi_over_12 = 2.0 * math.pi / 12.0
    c1 = [math.cos(two_pi_over_12 * float(x)) for x in range(12)]
    s1 = [math.sin(two_pi_over_12 * float(x)) for x in range(12)]

    carrier = {
        "typed_sum": ["L2_Z12_real_v1"],
        "index_set": "Z_12",
        "basis_frame": ["c1", "s1"],
        "basis_functions_by_x": {
            "c1": c1,
            "s1": s1,
        },
        "notes": [
            "Treat per-site arrays on Z_12 as elements of a real 12D carrier.",
            "This packet freezes only a future orientation-export target frame, not an exported E_orient object.",
        ],
    }

    checks = [
        {
            "id": "n540_first_clause_discharged",
            "actual": bool(n540["theorem_result"]["discharged"]),
            "expected": True,
            "pass": bool(n540["theorem_result"]["discharged"]) is True,
            "meaning": "N540 discharges the first admissibility clause (genuinely new strict-core source object).",
        }
    ]

    artifact = {
        "stage": "F649",
        "lane": "second_clause_typing_only",
        "goal": "freeze_minimal_typed_carrier_and_future_orientation_export_target_frame_for_S_sel_int_strict_core_source_object_v1",
        "status": "F649_EXECUTED_FIRST_EXPORTED_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_SECOND_CLAUSE_TYPING_PACKET_NO_FALSE_PASS",
        "exported_object": exported_object,
        "carrier": carrier,
        "state_support": {
            "nonzero_s1_support": dot != 0.0,
            "dot_w_break_with_s1": dot,
        },
        "future_orientation_export_target_frame": ["c1", "s1"],
        "E_orient_exported": False,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F649",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "exported_object": exported_object,
        "carrier": {
            "typed_sum": carrier["typed_sum"],
            "basis_frame": carrier["basis_frame"],
            "index_set": carrier["index_set"],
        },
        "state_support": artifact["state_support"],
        "future_orientation_export_target_frame": artifact["future_orientation_export_target_frame"],
        "E_orient_exported": False,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

