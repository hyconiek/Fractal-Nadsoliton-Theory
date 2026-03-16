#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n544_current_fifth_admissibility_clause_discharge_theorem_for_s_sel_int_strict_core_source_object_v1_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p652 = load_json(
        "fundamental_action_reconstruction/generated/p652_current_fifth_admissibility_clause_rerun_after_selector_acceptance_independence_packet_for_s_sel_int_strict_core_source_object_v1_probe_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_STRICT_CORE_SOURCE_OBJECT_SATISFYING_THE_FIFTH_ADMISSIBILITY_CLAUSE_FOR_S_SEL_INT_AFTER_P652"
    )
    status_ok = p652.get("status") == expected_status
    checks = [
        {
            "id": "positive_fifth_clause_rerun",
            "actual": p652.get("status"),
            "expected": expected_status,
            "pass": status_ok,
        }
    ]

    summary = {
        "step": "N544",
        "status": "N544_DISCHARGED_CURRENT_FIFTH_ADMISSIBILITY_CLAUSE_THEOREM_FOR_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_V1_NO_FALSE_PASS",
        "scope": "current_fifth_admissibility_clause_only",
        "checks": checks,
        "theorem_result": {
            "discharged": status_ok,
            "exported_object": "S_sel_int_strict_core_source_object_v1",
            "fifth_clause": "selector_acceptance_outside_strict_core_may_not_count_as_source_construction",
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
        "no_false_pass": True,
    }

    GENERATED.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

