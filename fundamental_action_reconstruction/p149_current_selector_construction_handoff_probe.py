#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p149_current_selector_construction_handoff_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f62 = load_json(
        "fundamental_action_reconstruction/generated/f62_current_selector_construction_handoff_packet_summary.json"
    )
    n164 = load_json(
        "fundamental_action_reconstruction/generated/n164_current_selector_construction_stop_condition_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "f62_handoff_active",
            "actual": f62["handoff"]["current_selector_construction_handoff_active"],
            "expected": True,
            "meaning": "F62 already activates the handoff packet",
        },
        {
            "id": "n164_lane_stopped",
            "actual": n164["theorem_result"][
                "current_selector_construction_lane_should_be_treated_as_stopped"
            ],
            "expected": True,
            "meaning": "N164 already stops the current selector-construction lane",
        },
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
            "stage": "P149",
            "lane": "current_selector_construction_handoff_probe_only",
            "status": "P149_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_HANDOFF_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P149",
            "lane": "current_selector_construction_handoff_probe_only",
            "goal": "test_whether_the_current_repo_supports_handing_off_the_stopped_selector_construction_lane_to_genuinely_additive_future_upstream_source_work",
            "status": "CURRENT_REPO_SUPPORTS_THE_CONCLUSION_THAT_THE_STOPPED_SELECTOR_CONSTRUCTION_LANE_SHOULD_BE_HANDED_OFF_TO_GENUINELY_ADDITIVE_FUTURE_UPSTREAM_SOURCE_WORK_AFTER_P149",
            "checks": checks,
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
