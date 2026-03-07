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
    / "f19_legacy_fine_structure_method_successor_semantics_refinement_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f19_legacy_fine_structure_method_successor_semantics_refinement_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def entries_have_method(payload: dict[str, Any], method_name: str) -> bool:
    return any(item.get("method") == method_name for item in payload.get("entries", []))


def updates_have_method(payload: dict[str, Any], method_name: str) -> bool:
    return any(item.get("method") == method_name for item in payload.get("updates", []))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n94 = load_json(
        "fundamental_action_reconstruction/generated/n94_current_strict_side_object_successor_branch_full_negative_closure_theorem_for_legacy_fine_structure_role_summary.json"
    )
    qw2069 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2069_full_sm_gr_derivation_package.json"
    )
    qw2098 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2098_ew_secondary_nonanchor_closure_gate.json"
    )
    qw2094 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2094_strict_rigor_defect_sweep.json"
    )

    method_name = "qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r"

    method_chain_present = (
        entries_have_method(qw2069, method_name)
        and updates_have_method(qw2098, method_name)
    )
    method_chain_consistency_present = (
        qw2094.get("checks", {}).get("QW-2069_alpha_em_inv_mz_matches_qw2098_method")
        is True
    )

    textual_method_successor_verdict_present = False
    method_lineage_upgrade_verdict_present = False

    checks_spec = [
        {
            "id": "n94_object_branch_closed",
            "actual": n94["theorem_result"]["object_successor_branch_closed_negatively_on_current_repo_state"],
            "expected": True,
            "meaning": "N94 already closes the full object-successor branch negatively on the current repo state",
        },
        {
            "id": "method_chain_present",
            "actual": method_chain_present,
            "expected": True,
            "meaning": "the repo exports the qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r method chain across QW-2069 and QW-2098",
        },
        {
            "id": "method_chain_consistency_present",
            "actual": method_chain_consistency_present,
            "expected": True,
            "meaning": "QW-2094 confirms method-side consistency of the qw2098 alpha chain",
        },
        {
            "id": "textual_method_successor_verdict_present",
            "actual": textual_method_successor_verdict_present,
            "expected": False,
            "meaning": "no explicit textual method-successor-semantics verdict is currently exported",
        },
        {
            "id": "method_lineage_upgrade_verdict_present",
            "actual": method_lineage_upgrade_verdict_present,
            "expected": False,
            "meaning": "no explicit method-lineage-upgrade verdict is currently exported",
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
        "stage": "F19",
        "lane": "legacy_fine_structure_method_successor_semantics_refinement_current_repo_state_only",
        "goal": "refine_the_missing_method_successor_semantics_verdict_into_textual_vs_lineage_upgrade_subbranches",
        "status": "F19_EXECUTED_LEGACY_FINE_STRUCTURE_METHOD_SUCCESSOR_SEMANTICS_REFINEMENT_PACKET_NO_FALSE_PASS",
        "reason": "N94 closes the object-successor branch negatively, while the repo already exports a real qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r method chain and QW-2094 method-side consistency markers; therefore the narrowest honest refinement is textual method-successor verdict vs method-lineage-upgrade verdict",
        "candidate_state": {
            "method_chain_present": method_chain_present,
            "method_chain_consistency_present": method_chain_consistency_present,
            "textual_method_successor_verdict_present": textual_method_successor_verdict_present,
            "method_lineage_upgrade_verdict_present": method_lineage_upgrade_verdict_present,
        },
        "remaining_missing_objects": [
            "explicit_textual_method_successor_semantics_verdict_identifying_qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r_as_the_strict_side_successor_semantics_replacing_the_legacy_fine_structure_role",
            "explicit_method_lineage_upgrade_verdict_elevating_the_existing_qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r_chain_into_replacement_semantics_for_the_legacy_fine_structure_role",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F19",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "candidate_state": artifact["candidate_state"],
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
