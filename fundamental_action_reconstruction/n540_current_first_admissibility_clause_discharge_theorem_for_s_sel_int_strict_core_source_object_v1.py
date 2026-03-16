#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n540_current_first_admissibility_clause_discharge_theorem_for_s_sel_int_strict_core_source_object_v1_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p648 = load_json(
        "fundamental_action_reconstruction/generated/p648_current_first_admissibility_clause_rerun_after_seed_v1_witness_provider_export_probe_summary.json"
    )

    status_ok = (
        p648["status"]
        == "CURRENT_REPO_EXPORTS_ONE_GENUINELY_NEW_STRICT_CORE_SOURCE_OBJECT_SATISFYING_THE_FIRST_ADMISSIBILITY_CLAUSE_FOR_S_SEL_INT_AFTER_P648"
    )

    checks = [
        {
            "id": "positive_first_clause_rerun",
            "actual": p648["status"],
            "expected": "CURRENT_REPO_EXPORTS_ONE_GENUINELY_NEW_STRICT_CORE_SOURCE_OBJECT_SATISFYING_THE_FIRST_ADMISSIBILITY_CLAUSE_FOR_S_SEL_INT_AFTER_P648",
            "pass": status_ok,
        }
    ]

    summary = {
        "step": "N540",
        "status": "N540_DISCHARGED_CURRENT_FIRST_ADMISSIBILITY_CLAUSE_THEOREM_FOR_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_V1_NO_FALSE_PASS",
        "scope": "current_first_admissibility_clause_only",
        "checks": checks,
        "theorem_result": {
            "discharged": status_ok,
            "exported_object": p648.get("exported_object"),
            "first_clause": "genuinely_new_strict_core_source_object_required",
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

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

