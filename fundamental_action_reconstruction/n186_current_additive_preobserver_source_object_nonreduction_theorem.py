#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n186_current_additive_preobserver_source_object_nonreduction_theorem_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f79 = load_json(
        "fundamental_action_reconstruction/generated/f79_first_additive_preobserver_source_object_nonreduction_comparison_packet_summary.json"
    )
    p167 = load_json(
        "fundamental_action_reconstruction/generated/p167_current_additive_preobserver_source_object_nonreduction_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "packaged_realization_present",
            "actual": f79["packaged_realization"],
            "expected": "S_preLM_target_packaged_realization_v1(d_*=1)",
        },
        {
            "id": "delta_nonzero_packet",
            "actual": f79["delta_is_nonzero"],
            "expected": True,
        },
        {
            "id": "delta_nonzero_probe",
            "actual": p167["status"],
            "expected": "CURRENT_REPO_EXPORTS_ONE_EXPLICIT_NONZERO_NONREDUCTION_WITNESS_BETWEEN_S_PRELM_ADDITIVE_CANDIDATE_V1_AND_THE_SAME_BASIS_PACKAGED_F75_TARGET_REALIZATION_AFTER_P167",
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
            "step": "N186",
            "status": "N186_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_NONREDUCTION_THEOREM_STATE",
            "scope": "current_additive_preobserver_source_object_nonreduction_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N186",
            "status": "N186_DISCHARGED_CURRENT_ADDITIVE_PREOBSERVER_SOURCE_OBJECT_NONREDUCTION_THEOREM_NO_FALSE_PASS",
            "scope": "current_additive_preobserver_source_object_nonreduction_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "construction_attempt": "S_preLM_additive_candidate_v1",
                "packaged_realization": "S_preLM_target_packaged_realization_v1(d_*=1)",
                "delta": f79["delta"],
                "delta_norm": f79["delta_norm"],
                "full_closure_pass": False,
            },
            "hard_limits": [
                "first_clause_not_yet_satisfied",
                "constructed_source_object_not_yet_exported",
                "admissible_S_sel_int_not_yet_constructed",
                "admissible_E_orient_not_yet_constructed",
                "downstream_chain_not_yet_constructed",
                "no_strict_core_selector_closure",
                "no_QW2191_discharge",
                "no_ToE_closure",
            ],
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
