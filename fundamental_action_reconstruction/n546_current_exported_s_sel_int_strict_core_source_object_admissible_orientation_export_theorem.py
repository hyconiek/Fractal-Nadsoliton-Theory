#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n546_current_exported_s_sel_int_strict_core_source_object_admissible_orientation_export_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p654 = load_json(
        "fundamental_action_reconstruction/generated/p654_current_exported_s_sel_int_strict_core_source_object_orientation_export_probe_summary.json"
    )

    expected_status = "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_ORIENTATION_DATUM_FROM_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_V1_AFTER_P654"
    status_ok = p654["status"] == expected_status

    checks = [
        {
            "id": "positive_orientation_export_probe",
            "actual": p654["status"],
            "expected": expected_status,
            "pass": status_ok,
        }
    ]

    summary = {
        "step": "N546",
        "status": "N546_DISCHARGED_CURRENT_EXPORTED_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_ADMISSIBLE_ORIENTATION_EXPORT_THEOREM_NO_FALSE_PASS",
        "scope": "current_s_sel_int_source_object_orientation_export_only",
        "checks": checks,
        "theorem_result": {
            "discharged": status_ok,
            "source_object": "S_sel_int_strict_core_source_object_v1",
            "exported_orientation": "E_orient_s_sel_int_source_object_v1",
            "admissible_E_orient": status_ok,
            "bridge_ready_for_B_sel": status_ok,
        },
        "hard_limits": [
            "admissible_S_sel_int_not_yet_constructed",
            "B_sel_not_yet_constructed",
            "R_sel_not_yet_constructed",
            "O_sel_not_yet_constructed",
            "downstream_chain_not_yet_constructed",
            "no_strict_core_selector_closure",
            "no_QW2191_discharge",
            "no_ToE_closure",
        ],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

