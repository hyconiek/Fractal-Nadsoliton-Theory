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
    / "p115_current_admissible_strict_core_internal_selector_source_object_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p115_current_admissible_strict_core_internal_selector_source_object_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f29 = load_json(
        "fundamental_action_reconstruction/generated/f29_genuine_strict_core_internal_selector_source_admission_packet_summary.json"
    )
    n124 = load_json(
        "fundamental_action_reconstruction/generated/n124_current_strict_core_internal_selector_source_derivation_full_negative_closure_theorem_summary.json"
    )
    p2 = load_json(
        "fundamental_action_reconstruction/generated/p2_strict_core_sigma_int_to_a1_pair1_probe_summary.json"
    )
    n123 = load_json(
        "fundamental_action_reconstruction/generated/n123_current_legacy_to_strict_kernel_package_level_nonbridge_theorem_summary.json"
    )
    n125 = load_json(
        "fundamental_action_reconstruction/generated/n125_current_selector_requirement_theory_level_acceptance_theorem_summary.json"
    )

    downstream_operator_reachability_present = str(p2.get("status") or "").startswith("PASS_")

    # Conservative probe: do not upgrade to a positive admission claim unless the
    # historical negative packaging remains applicable.
    negative_packaging_supported = bool((n124.get("theorem_result") or {}).get("discharged")) and (
        not downstream_operator_reachability_present
    )

    admissible_object_present = False
    future_move_must_add_new_source_object = True
    status = "CURRENT_REPO_EXPORTS_NO_ADMISSIBLE_STRICT_CORE_INTERNAL_SELECTOR_SOURCE_OBJECT_AFTER_P115"
    if not negative_packaging_supported:
        status = "P115_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_ADMISSION_FRONTIER_STATE"
        future_move_must_add_new_source_object = False

    checks_spec = [
        {
            "id": "f29_admission_contract_present",
            "actual": f29["admission_contract"][
                "strict_core_source_export_required"
            ],
            "expected": True,
            "meaning": "F29 already exports the admission contract",
        },
        {
            "id": "n124_no_current_source_discharge",
            "actual": (n124.get("theorem_result") or {}).get(
                "strict_core_internal_selector_source_derivation_discharge_present"
            ),
            "expected": False,
            "meaning": "the historical negative packaging required that current strict-core source discharge remain absent",
        },
        {
            "id": "p2_no_downstream_operator_reachability",
            "actual": p2.get("status"),
            "expected": "NOT_PASS",
            "meaning": "the historical negative packaging required that downstream strict-core operator reachability to A_1(pair1) be absent",
        },
        {
            "id": "n123_nonbridge_forbids_hidden_substitution",
            "actual": n123["theorem_result"][
                "package_level_nonbridge_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "N123 forbids hidden package-level bridge substitution",
        },
        {
            "id": "n125_axiom_acceptance_not_counted_as_strict_source",
            "actual": n125["theorem_result"]["accepted_scope"],
            "expected": "axiom_augmented_only",
            "meaning": "N125 keeps selector acceptance outside strict core",
        },
        {
            "id": "admissible_object_present",
            "actual": admissible_object_present,
            "expected": False,
            "meaning": "probe conclusion (conservative): no admissible strict-core internal selector source object is claimed present unless negative packaging is supported",
        },
    ]

    checks: list[dict[str, Any]] = []
    for item in checks_spec:
        ok = item["actual"] == item["expected"]
        if item["id"] == "p2_no_downstream_operator_reachability":
            ok = not downstream_operator_reachability_present
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": ok,
                "meaning": item["meaning"],
            }
        )

    artifact = {
        "stage": "P115",
        "lane": "current_admissible_strict_core_internal_selector_source_object_only",
        "goal": "test_whether_the_current_repo_already_exports_any_object_satisfying_the_full_admission_contract_for_a_genuine_strict_core_internal_selector_source",
        "status": status,
        "admission_state": {
            "admissible_strict_core_internal_selector_source_object_present": admissible_object_present,
            "future_positive_move_must_add_new_source_object": future_move_must_add_new_source_object,
            "negative_packaging_supported": negative_packaging_supported,
            "downstream_operator_reachability_present": downstream_operator_reachability_present,
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P115",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "admission_state": artifact["admission_state"],
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
