#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p169_first_admissibility_clause_rerun_after_strict_core_export_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f81 = load_json(
        "fundamental_action_reconstruction/generated/f81_first_additive_preobserver_strict_core_source_object_export_packet_summary.json"
    )
    n184 = load_json(
        "fundamental_action_reconstruction/generated/n184_current_first_additive_preobserver_source_object_first_clause_admissibility_target_theorem_summary.json"
    )

    props = f81["strict_core_export_properties"]
    checks_spec = [
        {
            "id": "first_clause_target_fixed",
            "actual": n184["theorem_result"]["first_clause"],
            "expected": "genuinely_new_strict_core_source_object_required",
        },
        {
            "id": "new_exported_identity",
            "actual": props["genuinely_new_exported_identity"],
            "expected": True,
        },
        {
            "id": "constructed_source_object_exported",
            "actual": props["constructed_source_object_exported"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": props["strict_core_only"],
            "expected": True,
        },
        {
            "id": "upstream_of_observer",
            "actual": props["upstream_of_observer"],
            "expected": True,
        },
        {
            "id": "no_external_selector_import",
            "actual": props["no_external_selector_import"],
            "expected": True,
        },
        {
            "id": "reuses_psi0_artifact",
            "actual": props["reuses_psi0_artifact"],
            "expected": False,
        },
        {
            "id": "reuses_fr_route_artifact",
            "actual": props["reuses_fr_route_artifact"],
            "expected": False,
        },
        {
            "id": "reuses_sigma_int_candidate_as_source_object",
            "actual": props["reuses_sigma_int_candidate_as_source_object"],
            "expected": False,
        },
        {
            "id": "reuses_overlay_fit_artifact",
            "actual": props["reuses_overlay_fit_artifact"],
            "expected": False,
        },
        {
            "id": "uses_axiom_lane_artifact",
            "actual": props["uses_axiom_lane_artifact"],
            "expected": False,
        },
        {
            "id": "nonreduction_witness_present",
            "actual": f81["supporting_nonreduction_witness"]["present"],
            "expected": True,
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
            }
        )
        if not ok:
            mismatches.append(item["id"])

    if mismatches:
        summary = {
            "stage": "P169",
            "lane": "first_admissibility_clause_rerun_after_strict_core_export_only",
            "status": "P169_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FIRST_CLAUSE_STATE_AFTER_STRICT_CORE_EXPORT",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P169",
            "lane": "first_admissibility_clause_rerun_after_strict_core_export_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_GENUINELY_NEW_STRICT_CORE_SOURCE_OBJECT_SATISFYING_THE_FIRST_ADMISSIBILITY_CLAUSE_AFTER_P169",
            "exported_object": "S_preLM_strict_core_source_object_v1",
            "first_clause": "genuinely_new_strict_core_source_object_required",
            "checks": checks,
            "remaining_admissibility_clauses_unresolved": [
                "carrier_typed_enough_for_later_export",
                "source_seed_only",
                "strict_core_only_already_fixed_but_full_admissibility_not_granted",
                "non_substitutive_global_compatibility",
                "selector_acceptance_independent",
                "future_bridge_compatible",
            ],
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
