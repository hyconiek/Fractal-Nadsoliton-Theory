#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p161_current_preobserver_light_matter_asymmetric_provider_packet_consistency_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f73 = load_json(
        "fundamental_action_reconstruction/generated/f73_first_preobserver_light_matter_source_provider_class_target_packet_summary.json"
    )
    f74 = load_json(
        "fundamental_action_reconstruction/generated/f74_first_preobserver_light_matter_asymmetric_provider_packet_summary.json"
    )
    n163 = load_json(
        "fundamental_action_reconstruction/generated/n163_current_observer_information_deficit_downstream_symptom_theorem_summary.json"
    )
    n166 = load_json(
        "fundamental_action_reconstruction/generated/n166_current_only_honest_positive_work_contract_theorem_summary.json"
    )
    n179 = load_json(
        "fundamental_action_reconstruction/generated/n179_current_only_honest_positive_reopening_target_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "f73_provider_class_target",
            "actual": f73["provider_class_target"],
            "expected": "preobserver_light_matter_source_provider_class_v1",
            "meaning": "F73 already fixes the provider class target",
        },
        {
            "id": "f74_kernel_split_safe",
            "actual": [
                f74["strict_kernel_control"]["kernel_used_as_operational_control_only"],
                f74["strict_kernel_control"]["legacy_bridge_claimed"],
            ],
            "expected": [True, False],
            "meaning": "F74 keeps strict-kernel use operational and non-bridged",
        },
        {
            "id": "f74_one_way_cascade",
            "actual": [
                f74["cascade"]["strictly_one_way_by_packet_construction"],
                f74["cascade"]["reverse_maps_exported"],
            ],
            "expected": [True, False],
            "meaning": "F74 keeps the cascade one-way",
        },
        {
            "id": "f74_observer_nonparticipation",
            "actual": [
                f74["observer_nonparticipation"]["dP_NL_dO_zero"],
                f74["observer_nonparticipation"]["dP_LM_dO_zero"],
                f74["observer_nonparticipation"]["observer_to_upstream_blocks_zero"],
            ],
            "expected": [True, True, True],
            "meaning": "F74 keeps observer out of upstream couplings",
        },
        {
            "id": "n163_observer_downstream",
            "actual": n163["theorem_result"]["observer_information_deficit_downstream_symptom_on_current_repo_state"],
            "expected": True,
            "meaning": "N163 already keeps observer deficit downstream",
        },
        {
            "id": "n166_upstream_contract",
            "actual": [
                n166["theorem_result"]["remaining_positive_work_must_be_genuinely_additive"],
                n166["theorem_result"]["remaining_positive_work_must_stay_upstream_of_observer"],
                n166["theorem_result"]["remaining_positive_work_must_be_source_object_first"],
            ],
            "expected": [True, True, True],
            "meaning": "N166 already fixes the additive upstream contract",
        },
        {
            "id": "n179_reopening_target",
            "actual": n179["theorem_result"]["only_honest_positive_reopening_target"],
            "expected": "preobserver_light_matter_source_provider_class_v1",
            "meaning": "N179 already reduces the reopening target to this provider class",
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
            "stage": "P161",
            "lane": "future_preobserver_light_matter_asymmetric_provider_packet_consistency_only",
            "status": "P161_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_PROVIDER_PACKET_CONSISTENCY_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P161",
            "lane": "future_preobserver_light_matter_asymmetric_provider_packet_consistency_only",
            "goal": "test_whether_the_first_preobserver_light_matter_asymmetric_provider_packet_stays_guardrail_consistent",
            "status": "CURRENT_REPO_EXPORTS_ONE_INTERNAL_GUARDRAIL_CONSISTENT_PREOBSERVER_LIGHT_MATTER_ASYMMETRIC_PROVIDER_PACKET_AFTER_P161",
            "provider_packet": "preobserver_light_matter_source_provider_packet_v1",
            "checks": checks,
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
