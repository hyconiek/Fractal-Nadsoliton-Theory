#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p215_current_source_topology_invariant_candidate_packet_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f127 = load_json(
        GENERATED / "f127_first_source_topology_invariant_candidate_packet_summary.json"
    )
    n163 = load_json(
        GENERATED / "n163_current_observer_information_deficit_downstream_symptom_theorem_summary.json"
    )
    n234 = load_json(
        GENERATED
        / "n234_current_global_selector_closure_and_qw2191_discharge_promotion_obstruction_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "future_tau_src_candidate_packet_exported",
            "actual": f127["candidate_packet"] == "tau_src_candidate_v1",
            "expected": True,
        },
        {
            "id": "candidate_upstream_of_observer",
            "actual": f127["properties"]["upstream_of_observer"],
            "expected": True,
        },
        {
            "id": "no_external_selector_import",
            "actual": not f127["properties"]["external_selector_import"],
            "expected": True,
        },
        {
            "id": "observer_downstream_only",
            "actual": n163["theorem_result"][
                "observer_information_deficit_downstream_symptom_on_current_repo_state"
            ],
            "expected": True,
        },
        {
            "id": "basis_independent_promotion_not_exported",
            "actual": not f127["properties"]["basis_independent_selector_promotion_exported"],
            "expected": True,
        },
        {
            "id": "qw2191_quotient_safe_promotion_not_exported",
            "actual": not f127["properties"]["qw2191_quotient_safe_promotion_exported"],
            "expected": True,
        },
        {
            "id": "current_selector_closure_not_claimed",
            "actual": not f127["properties"]["current_selector_closure"],
            "expected": True,
        },
        {
            "id": "current_qw2191_discharge_not_claimed",
            "actual": not f127["properties"]["current_qw2191_discharge"]
            and not n234["theorem_result"]["global_qw2191_discharge_justified_on_current_repo_state"],
            "expected": True,
        },
    ]

    checks = []
    failed = []
    for item in checks_spec:
        ok = item["actual"] == item["expected"]
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": ok,
            }
        )
        if not ok:
            failed.append(item["id"])

    if failed:
        status = "CURRENT_REPO_DOES_NOT_YET_EXPORT_A_GUARDRAIL_CONSISTENT_FUTURE_SOURCE_TOPOLOGY_INVARIANT_CANDIDATE_PACKET_AFTER_P215"
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_FUTURE_SOURCE_TOPOLOGY_INVARIANT_CANDIDATE_PACKET_BELOW_BASIS_INDEPENDENCE_AND_QW2191_QUOTIENT_SAFETY_AFTER_P215"

    summary = {
        "stage": "P215",
        "lane": "future_source_topology_selector_route_only",
        "status": status,
        "checks": checks,
        "failed_checks": failed,
        "tau_src_candidate_packet_exported": not failed,
        "basis_independent_selector_promotion_exported": False,
        "qw2191_quotient_safe_promotion_exported": False,
        "current_selector_closure": False,
        "current_qw2191_discharge": False,
        "observer_downstream_only": n163["theorem_result"][
            "observer_information_deficit_downstream_symptom_on_current_repo_state"
        ],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
