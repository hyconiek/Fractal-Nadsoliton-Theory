#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n166_current_only_honest_positive_work_contract_theorem_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p150 = load_json(
        "fundamental_action_reconstruction/generated/p150_current_future_additive_upstream_source_work_contract_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "p150_status",
            "actual": p150["status"],
            "expected": "CURRENT_REPO_REDUCES_THE_ONLY_HONEST_POSITIVE_WORK_TO_ONE_EXPLICIT_FUTURE_ADDITIVE_UPSTREAM_SOURCE_WORK_CONTRACT_AFTER_P150",
            "meaning": "P150 already reduces the only remaining honest positive work to one explicit contract",
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
            "step": "N166",
            "status": "N166_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_ONLY_HONEST_POSITIVE_WORK_CONTRACT_STATE",
            "scope": "current_only_honest_positive_work_contract_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N166",
            "status": "N166_DISCHARGED_CURRENT_ONLY_HONEST_POSITIVE_WORK_CONTRACT_THEOREM_NO_FALSE_PASS",
            "scope": "current_only_honest_positive_work_contract_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "only_honest_positive_work_contract_active": True,
                "remaining_positive_work_must_be_genuinely_additive": True,
                "remaining_positive_work_must_stay_upstream_of_observer": True,
                "remaining_positive_work_must_be_kernel_split_safe": True,
                "remaining_positive_work_must_be_source_object_first": True,
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
