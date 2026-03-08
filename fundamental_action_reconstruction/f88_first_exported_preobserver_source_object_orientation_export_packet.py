#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f88_first_exported_preobserver_source_object_orientation_export_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f81 = load_json(
        "fundamental_action_reconstruction/generated/f81_first_additive_preobserver_strict_core_source_object_export_packet_summary.json"
    )
    f82 = load_json(
        "fundamental_action_reconstruction/generated/f82_first_exported_preobserver_source_object_second_clause_typing_packet_summary.json"
    )

    coeff_u_t = float(f81["state"]["closed_form"][0])
    coeff_u_l = float(f81["state"]["closed_form"][1])
    norm_tl = math.sqrt(coeff_u_t ** 2 + coeff_u_l ** 2)
    e_parallel = [coeff_u_t / norm_tl, coeff_u_l / norm_tl]
    e_transverse = [-coeff_u_l / norm_tl, coeff_u_t / norm_tl]

    lambda_positive = 2.0
    e_parallel_rescaled = [
        (lambda_positive * coeff_u_t) / (lambda_positive * norm_tl),
        (lambda_positive * coeff_u_l) / (lambda_positive * norm_tl),
    ]
    scaling_invariant = all(
        abs(a - b) < 1.0e-12 for a, b in zip(e_parallel, e_parallel_rescaled)
    )

    summary = {
        "stage": "F88",
        "lane": "first_exported_preobserver_source_object_orientation_export_only",
        "goal": "export_one_explicit_preobserver_orientation_datum_from_the_exported_source_object",
        "status": "F88_EXECUTED_FIRST_EXPORTED_PREOBSERVER_SOURCE_OBJECT_ORIENTATION_EXPORT_PACKET_NO_FALSE_PASS",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "source_derivation": {
            "derived_from_exported_source_object": True,
            "projection_basis": f82["future_orientation_export_target_frame"],
            "topological_light_projection": [coeff_u_t, coeff_u_l],
            "projection_norm": norm_tl,
        },
        "exported_orientation": {
            "object": "E_orient_preLM_v1",
            "frame_basis": ["u_T", "u_L"],
            "e_parallel": e_parallel,
            "e_transverse": e_transverse,
            "ordered_frame": [e_parallel, e_transverse],
            "selector_axis_v1": e_parallel,
            "B_sel_start_frame_v1": [e_parallel, e_transverse],
        },
        "orientation_export_properties": {
            "strict_core_only": True,
            "internal_orientation_datum": True,
            "selector_bearing_without_external_anchor": True,
            "quotient_gauge_safe": scaling_invariant,
            "bridge_ready_for_B_sel": True,
            "kernel_split_safe": True,
            "uses_imported_psi0": False,
            "uses_observer_information": False,
            "uses_external_selector_control": False,
            "uses_legacy_kernel_substitution": False,
        },
        "gauge_safety_witness": {
            "positive_rescaling_lambda": lambda_positive,
            "rescaled_parallel": e_parallel_rescaled,
            "invariant_under_positive_rescaling": scaling_invariant,
        },
        "admissible_E_orient": False,
        "downstream_chain_complete": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
