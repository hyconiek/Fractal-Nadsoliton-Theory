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
    / "p113_current_strict_core_internal_selector_source_derivation_discharge_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p113_current_strict_core_internal_selector_source_derivation_discharge_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f28 = load_json(
        "fundamental_action_reconstruction/generated/f28_current_strict_core_internal_selector_source_derivation_refinement_packet_summary.json"
    )
    p2 = load_json(
        "fundamental_action_reconstruction/generated/p2_strict_core_sigma_int_to_a1_pair1_probe_summary.json"
    )

    no_discharge_supported = (
        f28["branch_state"][
            "generic_hidden_source_branch_closed_negatively_on_current_repo_state"
        ]
        and f28["branch_state"]["psi0_branch_closed_negatively_on_current_repo_state"]
        and f28["branch_state"]["fr_branch_closed_negatively_on_current_repo_state"]
        and f28["branch_state"][
            "sigma_int_bridge_branch_closed_negatively_on_current_repo_state"
        ]
        and p2["status"] == "NOT_COMPUTABLE_FROM_CURRENT_STRICT_CORE_ROUTE"
    )

    checks_spec = [
        {
            "id": "f28_generic_branch_closed",
            "actual": f28["branch_state"][
                "generic_hidden_source_branch_closed_negatively_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "the generic hidden-source branch is already closed negatively",
        },
        {
            "id": "f28_psi0_branch_closed",
            "actual": f28["branch_state"][
                "psi0_branch_closed_negatively_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "the current psi0 branch is already closed negatively",
        },
        {
            "id": "f28_fr_branch_closed",
            "actual": f28["branch_state"][
                "fr_branch_closed_negatively_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "the current FR branch is already closed negatively",
        },
        {
            "id": "f28_sigma_int_branch_closed",
            "actual": f28["branch_state"][
                "sigma_int_bridge_branch_closed_negatively_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "the current sigma-int bridge branch is already closed negatively",
        },
        {
            "id": "p2_no_downstream_a1_pair1_reachability",
            "actual": p2["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_STRICT_CORE_ROUTE",
            "meaning": "there is still no strict-core downstream reachability to A1(pair1)",
        },
        {
            "id": "package_level_non_discharge_supported",
            "actual": no_discharge_supported,
            "expected": True,
            "meaning": "taken together, the current repo supports a package-level non-discharge conclusion for strict-core internal selector source derivation",
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

    artifact = {
        "stage": "P113",
        "lane": "current_strict_core_internal_selector_source_package_level_probe_only",
        "goal": "test_whether_the_current_repo_now_supports_an_explicit_package_level_non_discharge_conclusion_for_strict_core_internal_selector_source_derivation",
        "status": "CURRENT_REPO_DOES_NOT_EXPORT_AN_EXPLICIT_STRICT_CORE_INTERNAL_SELECTOR_SOURCE_DERIVATION_DISCHARGE_AFTER_P113",
        "reason": "B2, N4, N5, N6, N7, N8, and P2 together already close the generic hidden-source branch, the psi0 branch, the FR branch, the sigma-int bridge branch, and the downstream A1(pair1) reachability route negatively on the current repo state; therefore the current repo still exports no strict-core internal selector source derivation discharge",
        "support_state": {
            "generic_hidden_source_branch_closed_negatively": f28["branch_state"][
                "generic_hidden_source_branch_closed_negatively_on_current_repo_state"
            ],
            "psi0_branch_closed_negatively": f28["branch_state"][
                "psi0_branch_closed_negatively_on_current_repo_state"
            ],
            "fr_branch_closed_negatively": f28["branch_state"][
                "fr_branch_closed_negatively_on_current_repo_state"
            ],
            "sigma_int_bridge_branch_closed_negatively": f28["branch_state"][
                "sigma_int_bridge_branch_closed_negatively_on_current_repo_state"
            ],
            "strict_core_downstream_A1_pair1_reachability_present": False,
            "package_level_non_discharge_supported": no_discharge_supported,
        },
        "remaining_missing_objects": [
            "explicit_strict_core_internal_selector_source_derivation_discharge",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P113",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "support_state": artifact["support_state"],
        "remaining_missing_objects": artifact["remaining_missing_objects"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    print(OUT_JSON)


if __name__ == "__main__":
    main()
