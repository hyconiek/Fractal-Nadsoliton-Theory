#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f76_first_additive_preobserver_source_object_construction_attempt_packet_summary.json"


def main() -> None:
    phi = 0.16250
    beta = 1.0
    cos_phi = math.cos(phi)
    d_star = 1.0
    i_mat = beta * (d_star**1.8) / (1.0 + beta * (d_star**1.8))

    summary = {
        "stage": "F76",
        "lane": "future_additive_preobserver_source_object_construction_attempt_only",
        "goal": "export_one_explicit_additive_preobserver_source_object_construction_attempt_above_the_fixed_preobserver_target",
        "status": "F76_EXECUTED_FIRST_ADDITIVE_PREOBSERVER_SOURCE_OBJECT_CONSTRUCTION_ATTEMPT_PACKET_NO_FALSE_PASS",
        "source_object_target": "preobserver_light_matter_source_object_target_v1",
        "construction_attempt": "S_preLM_additive_candidate_v1",
        "carrier": {
            "used": ["V_topo", "L_int", "M_int"],
            "basis": ["u_T", "u_L", "u_M"],
            "observer_excluded": True,
        },
        "strict_kernel_role": {
            "operational_only": True,
            "identified_with_legacy_kernel": False,
            "d_star": d_star,
            "I_mat(d_star)": i_mat,
        },
        "generator": {
            "name": "A_up",
            "matrix": [
                [0.0, 0.0, 0.0],
                [cos_phi, 0.0, 0.0],
                [0.0, i_mat, 0.0],
            ],
            "nilpotent_order_leq_3": True,
        },
        "attempt_profile": {
            "definition": "exp(A_up) u_T",
            "closed_form": [
                1.0,
                cos_phi,
                cos_phi * i_mat / 2.0,
            ],
            "assembled_state": "u_T + cos(phi) u_L + (cos(phi)/4) u_M",
        },
        "guardrails": {
            "future_only": True,
            "additive_construction_attempt_only": True,
            "upstream_of_observer": True,
            "light_before_observer": True,
            "matter_as_encoding_not_primary_selector_source": True,
            "kernel_split_safe": True,
            "no_external_selector_import": True,
            "source_object_first": True,
        },
        "constructed_source_object": False,
        "admissible_S_sel_int": False,
        "admissible_E_orient": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
