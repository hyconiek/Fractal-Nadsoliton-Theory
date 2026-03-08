#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n165_current_selector_construction_handoff_theorem_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p149 = load_json(
        "fundamental_action_reconstruction/generated/p149_current_selector_construction_handoff_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "p149_status",
            "actual": p149["status"],
            "expected": "CURRENT_REPO_SUPPORTS_THE_CONCLUSION_THAT_THE_STOPPED_SELECTOR_CONSTRUCTION_LANE_SHOULD_BE_HANDED_OFF_TO_GENUINELY_ADDITIVE_FUTURE_UPSTREAM_SOURCE_WORK_AFTER_P149",
            "meaning": "P149 already supports the handoff conclusion",
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
            "step": "N165",
            "status": "N165_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_SELECTOR_CONSTRUCTION_HANDOFF_STATE",
            "scope": "current_selector_construction_handoff_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N165",
            "status": "N165_DISCHARGED_CURRENT_SELECTOR_CONSTRUCTION_HANDOFF_THEOREM_NO_FALSE_PASS",
            "scope": "current_selector_construction_handoff_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "stopped_selector_construction_lane_handed_off": True,
                "remaining_positive_move_class": "future_genuinely_additive_upstream_source_work_only",
                "full_closure_pass": False,
            },
            "hard_limits": [
                "no_permanent_impossibility_claim",
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
