#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n167_current_first_contract_compliant_additive_upstream_source_target_theorem_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p151 = load_json(
        "fundamental_action_reconstruction/generated/p151_first_contract_compliant_additive_upstream_source_target_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "p151_status",
            "actual": p151["status"],
            "expected": "CURRENT_REPO_REDUCES_THE_ONLY_HONEST_POSITIVE_WORK_TO_ONE_FIRST_CONTRACT_COMPLIANT_FUTURE_ADDITIVE_UPSTREAM_SOURCE_TARGET_AFTER_P151",
            "meaning": "P151 already reduces the only remaining honest positive work to one first contract-compliant target",
        }
    ]

    checks = []
    mismatches = []
    for item in checks_spec:
        ok = item["actual"] == item["expected"]
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": ok,
                "meaning": item["meaning"],
            }
        )
        if not ok:
            mismatches.append(item["id"])

    if mismatches:
        summary = {
            "step": "N167",
            "status": "N167_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FIRST_CONTRACT_COMPLIANT_TARGET_STATE",
            "scope": "current_first_contract_compliant_additive_upstream_source_target_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N167",
            "status": "N167_DISCHARGED_CURRENT_FIRST_CONTRACT_COMPLIANT_ADDITIVE_UPSTREAM_SOURCE_TARGET_THEOREM_NO_FALSE_PASS",
            "scope": "current_first_contract_compliant_additive_upstream_source_target_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "first_contract_compliant_additive_upstream_source_target_active": True,
                "first_target": "S_sel_int_future_additive_upstream_target_v2",
                "full_closure_pass": False,
            },
            "hard_limits": [
                "future_additive_source_object_not_yet_constructed",
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
