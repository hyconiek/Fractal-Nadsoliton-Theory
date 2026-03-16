#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

OUT_JSON = (
    GENERATED
    / "p648_current_first_admissibility_clause_rerun_after_seed_v1_witness_provider_export_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p648_current_first_admissibility_clause_rerun_after_seed_v1_witness_provider_export_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f34 = load_json(
        "fundamental_action_reconstruction/generated/f34_minimal_admissible_strict_core_source_seed_construction_contract_packet_summary.json"
    )
    f647 = load_json(
        "fundamental_action_reconstruction/generated/f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt.json"
    )
    f647_summary = load_json(
        "fundamental_action_reconstruction/generated/f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt_summary.json"
    )

    props = f647_summary["strict_core_export_properties"]
    exported_object = str(f647_summary["exported_constructed_source_object"])

    checks_spec = [
        {
            "id": "f34_first_clause_required",
            "actual": f34["minimal_source_seed_construction_contract"][
                "genuinely_new_strict_core_source_object_required"
            ],
            "expected": True,
            "meaning": "F34 keeps the genuinely-new strict-core source-object clause active for S_sel_int",
        },
        {
            "id": "new_exported_identity",
            "actual": props["genuinely_new_exported_identity"],
            "expected": True,
            "meaning": "a new strict-core exported source-object identity exists",
        },
        {
            "id": "constructed_source_object_exported",
            "actual": props["constructed_source_object_exported"],
            "expected": True,
            "meaning": "the constructed source object is exported (seed-v1 witness provider)",
        },
        {
            "id": "strict_core_only",
            "actual": props["strict_core_only"],
            "expected": True,
            "meaning": "the object is exported on strict core only",
        },
        {
            "id": "no_external_selector_import",
            "actual": props["no_external_selector_import"],
            "expected": True,
            "meaning": "no external selector import is used as the source object",
        },
        {
            "id": "kernel_split_safe",
            "actual": props["kernel_split_safe"],
            "expected": True,
            "meaning": "no silent legacy->strict kernel substitution is used",
        },
        {
            "id": "reuses_psi0_artifact",
            "actual": props["reuses_psi0_artifact"],
            "expected": False,
            "meaning": "does not reuse psi0 artifacts as the source object",
        },
        {
            "id": "reuses_fr_route_artifact",
            "actual": props["reuses_fr_route_artifact"],
            "expected": False,
            "meaning": "does not reuse hybrid/legacy FR route artifacts as the source object",
        },
        {
            "id": "reuses_sigma_int_candidate_as_source_object",
            "actual": props["reuses_sigma_int_candidate_as_source_object"],
            "expected": False,
            "meaning": "does not identify the source object with sigma_int itself",
        },
        {
            "id": "reuses_overlay_fit_artifact",
            "actual": props["reuses_overlay_fit_artifact"],
            "expected": False,
            "meaning": "does not reuse overlay-fit artifacts as the source object",
        },
        {
            "id": "uses_axiom_lane_artifact",
            "actual": props["uses_axiom_lane_artifact"],
            "expected": False,
            "meaning": "does not rely on axiom-lane artifacts as the source object",
        },
        {
            "id": "nonreduction_witness_present",
            "actual": bool(f647.get("supporting_nonreduction_witness", {}).get("present")),
            "expected": True,
            "meaning": "a simple nonreduction witness is present (not pure packaging of sigma_int)",
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
        status = "P648_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FIRST_ADMISSIBILITY_CLAUSE_STATE_FOR_S_SEL_INT_AFTER_SEED_V1_WITNESS_PROVIDER_EXPORT"
        artifact = {
            "stage": "P648",
            "lane": "current_first_admissibility_clause_rerun_after_seed_v1_witness_provider_export_only",
            "goal": "rerun_the_first_admissibility_clause_for_S_sel_int_after_exporting_one_seed_v1_constructed_source_object_witness_provider",
            "status": status,
            "exported_object": exported_object,
            "first_clause": "genuinely_new_strict_core_source_object_required",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "no_false_pass": True,
        }
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_GENUINELY_NEW_STRICT_CORE_SOURCE_OBJECT_SATISFYING_THE_FIRST_ADMISSIBILITY_CLAUSE_FOR_S_SEL_INT_AFTER_P648"
        artifact = {
            "stage": "P648",
            "lane": "current_first_admissibility_clause_rerun_after_seed_v1_witness_provider_export_only",
            "goal": "rerun_the_first_admissibility_clause_for_S_sel_int_after_exporting_one_seed_v1_constructed_source_object_witness_provider",
            "status": status,
            "exported_object": exported_object,
            "first_clause": "genuinely_new_strict_core_source_object_required",
            "checks": checks,
            "remaining_admissibility_clauses_unresolved": [
                "carrier_typed_enough_for_later_E_orient_export_required",
                "source_seed_only_not_counted_as_E_orient_or_bridge",
                "strict_core_only_required (still not full admissibility)",
                "selector_acceptance_outside_strict_core_may_not_count_as_source_construction",
                "future_bridge_compatible_required",
            ],
            "no_false_pass": True,
        }

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "lane": artifact["lane"],
        "exported_object": artifact["exported_object"],
        "first_clause": artifact["first_clause"],
        "checks": artifact["checks"],
        "no_false_pass": True,
    }
    if "remaining_admissibility_clauses_unresolved" in artifact:
        summary["remaining_admissibility_clauses_unresolved"] = artifact[
            "remaining_admissibility_clauses_unresolved"
        ]
    if "blocking_mismatches" in artifact:
        summary["blocking_mismatches"] = artifact["blocking_mismatches"]

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

