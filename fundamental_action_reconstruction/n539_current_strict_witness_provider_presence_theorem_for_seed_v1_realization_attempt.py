#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n539_current_strict_witness_provider_presence_theorem_for_seed_v1_realization_attempt_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p646 = load_json(
        "fundamental_action_reconstruction/generated/p646_current_strict_witness_provider_scan_probe_for_seed_v1_realization_attempt_summary.json"
    )

    candidates_count = int(p646["scan_result"]["candidates_count"])
    witness_provider_present = bool(p646.get("witness_provider_present"))

    checks_spec = [
        {
            "id": "p646_witness_provider_present",
            "actual": witness_provider_present,
            "expected": True,
            "meaning": "P646 reports at least one exported artifact matching the F646 signature",
        },
        {
            "id": "p646_candidates_count_ge_1",
            "actual": candidates_count >= 1,
            "expected": True,
            "meaning": "P646 candidate list is nonempty",
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
            "step": "N539",
            "status": "N539_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_WITNESS_PROVIDER_PRESENCE_STATE_FOR_SEED_V1",
            "scope": "current_strict_witness_provider_presence_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
            "no_false_pass": True,
        }
    else:
        first = p646["scan_result"]["candidates"][0]
        summary = {
            "step": "N539",
            "status": "N539_DISCHARGED_CURRENT_STRICT_WITNESS_PROVIDER_PRESENCE_THEOREM_FOR_SEED_V1_REALIZATION_ATTEMPT_NO_FALSE_PASS",
            "scope": "current_strict_witness_provider_presence_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "strict_witness_provider_present": True,
                "candidates_count": candidates_count,
                "one_witness_provider_path": first.get("path"),
                "exported_constructed_source_object": first.get("exported_constructed_source_object"),
                "full_closure_pass": False,
            },
            "recommended_next_move": "attack_admissibility_branch_and_orientation_export_on_the_exported_constructed_source_object_without_implying_selector_closure",
            "hard_limits": [
                "no_admissible_S_sel_int",
                "no_admissible_E_orient",
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

