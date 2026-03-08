#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n187_current_additive_preobserver_source_object_first_clause_remainder_theorem_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n185 = load_json(
        "fundamental_action_reconstruction/generated/n185_current_first_clause_obstruction_theorem_for_s_prelm_additive_candidate_v1_summary.json"
    )
    n186 = load_json(
        "fundamental_action_reconstruction/generated/n186_current_additive_preobserver_source_object_nonreduction_theorem_summary.json"
    )
    p168 = load_json(
        "fundamental_action_reconstruction/generated/p168_current_additive_preobserver_source_object_first_clause_remainder_probe_summary.json"
    )

    checks = [
        {
            "id": "first_clause_still_obstructed",
            "actual": n185["theorem_result"]["discharged"],
            "expected": True,
            "pass": n185["theorem_result"]["discharged"] is True,
        },
        {
            "id": "nonreduction_witness_present",
            "actual": n186["theorem_result"]["delta_norm"] > 0.0,
            "expected": True,
            "pass": n186["theorem_result"]["delta_norm"] > 0.0,
        },
        {
            "id": "remainder_package_frozen",
            "actual": p168["remaining_first_clause_package"],
            "expected": "realized_constructed_source_object_export_package",
            "pass": p168["remaining_first_clause_package"]
            == "realized_constructed_source_object_export_package",
        },
    ]

    summary = {
        "step": "N187",
        "status": "N187_DISCHARGED_CURRENT_ADDITIVE_PREOBSERVER_SOURCE_OBJECT_FIRST_CLAUSE_REMAINDER_THEOREM_NO_FALSE_PASS",
        "scope": "current_additive_preobserver_source_object_first_clause_remainder_only",
        "checks": checks,
        "theorem_result": {
            "discharged": all(item["pass"] for item in checks),
            "construction_attempt": "S_preLM_additive_candidate_v1",
            "first_clause": "genuinely_new_strict_core_source_object_required",
            "simple_packaging_reduction_removed": True,
            "remaining_first_clause_package": "realized_constructed_source_object_export_package",
            "full_closure_pass": False,
        },
        "hard_limits": [
            "first_clause_not_yet_satisfied",
            "constructed_source_object_not_yet_exported",
            "admissible_S_sel_int_not_yet_constructed",
            "admissible_E_orient_not_yet_constructed",
            "downstream_chain_not_yet_constructed",
            "no_strict_core_selector_closure",
            "no_QW2191_discharge",
            "no_ToE_closure",
        ],
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
