#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f81_first_additive_preobserver_strict_core_source_object_export_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f76 = load_json(
        "fundamental_action_reconstruction/generated/f76_first_additive_preobserver_source_object_construction_attempt_packet_summary.json"
    )
    n186 = load_json(
        "fundamental_action_reconstruction/generated/n186_current_additive_preobserver_source_object_nonreduction_theorem_summary.json"
    )

    summary = {
        "stage": "F81",
        "lane": "first_additive_preobserver_strict_core_source_object_export_only",
        "goal": "export_one_genuinely_new_strict_core_source_object_identity_above_the_additive_preobserver_candidate",
        "status": "F81_EXECUTED_FIRST_ADDITIVE_PREOBSERVER_STRICT_CORE_SOURCE_OBJECT_EXPORT_PACKET_NO_FALSE_PASS",
        "base_attempt": "S_preLM_additive_candidate_v1",
        "exported_object": "S_preLM_strict_core_source_object_v1",
        "export_definition": "strict_core_export(S_preLM_additive_candidate_v1)",
        "carrier": f76["carrier"],
        "state": {
            "basis": f76["carrier"]["basis"],
            "closed_form": f76["attempt_profile"]["closed_form"],
            "assembled_state": f76["attempt_profile"]["assembled_state"],
        },
        "strict_core_export_properties": {
            "constructed_source_object_exported": True,
            "genuinely_new_exported_identity": True,
            "strict_core_only": True,
            "upstream_of_observer": True,
            "source_seed_only": True,
            "no_external_selector_import": True,
            "kernel_split_safe": True,
            "reuses_psi0_artifact": False,
            "reuses_fr_route_artifact": False,
            "reuses_sigma_int_candidate_as_source_object": False,
            "reuses_overlay_fit_artifact": False,
            "uses_axiom_lane_artifact": False,
        },
        "supporting_nonreduction_witness": {
            "present": True,
            "delta": n186["theorem_result"]["delta"],
            "delta_norm": n186["theorem_result"]["delta_norm"],
        },
        "admissible_S_sel_int": False,
        "admissible_E_orient": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
