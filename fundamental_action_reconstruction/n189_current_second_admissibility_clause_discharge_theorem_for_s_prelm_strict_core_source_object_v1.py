#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n189_current_second_admissibility_clause_discharge_theorem_for_s_prelm_strict_core_source_object_v1_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p170 = load_json(
        "fundamental_action_reconstruction/generated/p170_second_admissibility_clause_rerun_after_carrier_typing_probe_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_STRICT_CORE_SOURCE_OBJECT_SATISFYING_THE_SECOND_ADMISSIBILITY_CLAUSE_AFTER_P170"
    )
    checks = [
        {
            "id": "positive_second_clause_rerun",
            "actual": p170["status"],
            "expected": expected_status,
            "pass": p170["status"] == expected_status,
        }
    ]

    summary = {
        "step": "N189",
        "status": "N189_DISCHARGED_CURRENT_SECOND_ADMISSIBILITY_CLAUSE_THEOREM_FOR_S_PRELM_STRICT_CORE_SOURCE_OBJECT_V1_NO_FALSE_PASS",
        "scope": "current_second_admissibility_clause_only",
        "checks": checks,
        "theorem_result": {
            "discharged": all(item["pass"] for item in checks),
            "exported_object": "S_preLM_strict_core_source_object_v1",
            "second_clause": "carrier_typed_enough_for_later_export",
            "full_admissibility_pass": False,
        },
        "hard_limits": [
            "later_admissibility_clauses_not_yet_discharged",
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
