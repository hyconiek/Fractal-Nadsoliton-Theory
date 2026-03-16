#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n538_current_strict_witness_provider_absence_theorem_for_seed_v1_realization_attempt_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p646 = load_json(
        "fundamental_action_reconstruction/generated/p646_current_strict_witness_provider_scan_probe_for_seed_v1_realization_attempt_summary.json"
    )

    checks_spec = [
        {
            "id": "p646_candidate_count_zero",
            "actual": p646["scan_result"]["candidates_count"],
            "expected": 0,
            "meaning": "P646 finds no strict witness provider matching the F646 signature on current repo state",
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
            "step": "N538",
            "status": "N538_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_WITNESS_PROVIDER_SCAN_STATE_FOR_SEED_V1",
            "scope": "current_strict_witness_provider_scan_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
            "no_false_pass": True,
        }
    else:
        summary = {
            "step": "N538",
            "status": "N538_DISCHARGED_CURRENT_STRICT_WITNESS_PROVIDER_ABSENCE_THEOREM_FOR_SEED_V1_REALIZATION_ATTEMPT_NO_FALSE_PASS",
            "scope": "current_strict_witness_provider_scan_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "strict_witness_provider_present": False,
                "full_closure_pass": False,
            },
            "recommended_next_move": "implement_one_strict_witness_provider_export_packet_meeting_F646_contract",
            "hard_limits": [
                "no_constructed_source_object_exported_for_seed_v1_attempt",
                "no_admissible_S_sel_int",
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

