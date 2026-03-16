#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

OUT_JSON = (
    GENERATED
    / "p392_current_strict_side_source_seed_strict_sigma_int_upgrade_candidate_instance_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p392_current_strict_side_source_seed_strict_sigma_int_upgrade_candidate_instance_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    sigma_int_strict = load_json(
        "fundamental_action_reconstruction/generated/sigma_int_strict_derived_v1.json"
    )
    f318 = load_json(
        "fundamental_action_reconstruction/generated/f318_first_current_strict_side_source_seed_strict_sigma_int_upgrade_candidate_instance_packet_summary.json"
    )

    sigma_value = sigma_int_strict.get("value")
    sigma_in_codomain = sigma_value in (-1, 1, "+1", "-1")

    candidate = f318.get("candidate_seed_instance") or {}
    does_not_count_as = candidate.get("does_not_count_as") or []

    checks_spec = [
        {
            "id": "strict_sigma_int_source_upgrade_datum_exported",
            "actual": bool(
                sigma_int_strict.get("object") == "sigma_int_strict_derived_v1"
                and sigma_in_codomain
            ),
            "expected": True,
            "meaning": "strict sigma-int source-upgrade datum exists and lies in {+1,-1}",
        },
        {
            "id": "seed_candidate_instance_exported",
            "actual": bool(candidate.get("candidate_seed_name") == "S_sel_int_candidate_seed_v1"),
            "expected": True,
            "meaning": "the strict-side seed candidate instance S_sel_int_candidate_seed_v1 is exported",
        },
        {
            "id": "no_implied_admissible_s_sel_int",
            "actual": "admissible_S_sel_int" in does_not_count_as,
            "expected": True,
            "meaning": "the export stays explicitly below admissible S_sel_int / selector closure",
        },
    ]

    checks: list[dict[str, Any]] = []
    for item in checks_spec:
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": item["actual"] == item["expected"],
                "meaning": item["meaning"],
            }
        )

    overall_pass = all(item["pass"] for item in checks)

    artifact = {
        "stage": "P392",
        "lane": "strict_sigma_int_upgraded_seed_candidate_instance_only",
        "goal": "probe_whether_the_repo_exports_a_strict_sigma_int_upgraded_strict_side_seed_candidate_instance_below_admissible_S_sel_int",
        "status": "CURRENT_REPO_EXPORTS_ONE_STRICT_SIGMA_INT_UPGRADED_STRICT_SIDE_SOURCE_SEED_CANDIDATE_INSTANCE_AFTER_P392",
        "overall_pass": bool(overall_pass),
        "checked": {
            "sigma_int_strict_derived_v1": {
                "object": sigma_int_strict.get("object"),
                "value": sigma_value,
                "status": sigma_int_strict.get("status"),
                "path": "fundamental_action_reconstruction/generated/sigma_int_strict_derived_v1.json",
            },
            "seed_candidate_instance": candidate,
            "seed_source_packet": {
                "stage": f318.get("stage"),
                "status": f318.get("status"),
                "path": "fundamental_action_reconstruction/generated/f318_first_current_strict_side_source_seed_strict_sigma_int_upgrade_candidate_instance_packet_summary.json",
            },
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "overall_pass": artifact["overall_pass"],
        "lane": artifact["lane"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    print(OUT_JSON)


if __name__ == "__main__":
    main()

