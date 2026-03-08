#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n195_current_full_admissibility_discharge_theorem_for_s_prelm_strict_core_source_object_v1_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    theorems = {
        "first": load_json("fundamental_action_reconstruction/generated/n188_current_first_admissibility_clause_discharge_theorem_for_s_prelm_strict_core_source_object_v1_summary.json"),
        "second": load_json("fundamental_action_reconstruction/generated/n189_current_second_admissibility_clause_discharge_theorem_for_s_prelm_strict_core_source_object_v1_summary.json"),
        "third": load_json("fundamental_action_reconstruction/generated/n190_current_third_admissibility_clause_discharge_theorem_for_s_prelm_strict_core_source_object_v1_summary.json"),
        "fourth": load_json("fundamental_action_reconstruction/generated/n191_current_fourth_admissibility_clause_discharge_theorem_for_s_prelm_strict_core_source_object_v1_summary.json"),
        "fifth": load_json("fundamental_action_reconstruction/generated/n192_current_fifth_admissibility_clause_discharge_theorem_for_s_prelm_strict_core_source_object_v1_summary.json"),
        "sixth": load_json("fundamental_action_reconstruction/generated/n193_current_sixth_admissibility_clause_discharge_theorem_for_s_prelm_strict_core_source_object_v1_summary.json"),
        "seventh": load_json("fundamental_action_reconstruction/generated/n194_current_seventh_admissibility_clause_discharge_theorem_for_s_prelm_strict_core_source_object_v1_summary.json"),
    }

    checks = []
    for key, data in theorems.items():
        discharged = data["theorem_result"]["discharged"]
        checks.append(
            {
                "id": f"{key}_clause_discharged",
                "actual": discharged,
                "expected": True,
                "pass": discharged is True,
            }
        )

    summary = {
        "step": "N195",
        "status": "N195_DISCHARGED_CURRENT_FULL_ADMISSIBILITY_THEOREM_FOR_S_PRELM_STRICT_CORE_SOURCE_OBJECT_V1_NO_FALSE_PASS",
        "scope": "current_full_admissibility_only",
        "checks": checks,
        "theorem_result": {
            "discharged": all(item["pass"] for item in checks),
            "exported_object": "S_preLM_strict_core_source_object_v1",
            "all_seven_clauses_discharged": all(item["pass"] for item in checks),
            "admissible_S_sel_int": all(item["pass"] for item in checks),
            "admissible_E_orient": False,
            "full_selector_closure": False,
        },
        "hard_limits": [
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
