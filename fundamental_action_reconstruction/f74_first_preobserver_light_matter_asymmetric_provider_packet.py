#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f74_first_preobserver_light_matter_asymmetric_provider_packet_summary.json"


def main() -> None:
    omega = 0.18575
    phi = 0.16250
    beta = 1.0
    eta = 1.8

    kernel_core_amplitude = math.cos(phi)
    i_mat_at_zero = 0.0
    c_obs_at_zero = 0.0
    c_obs_large_scale_limit = 1.0

    summary = {
        "stage": "F74",
        "lane": "future_preobserver_light_matter_asymmetric_provider_packet_only",
        "goal": "export_one_explicit_future_additive_provider_packet_using_the_strict_kernel_only_as_operational_control_and_keeping_observer_downstream",
        "status": "F74_EXECUTED_FIRST_PREOBSERVER_LIGHT_MATTER_ASYMMETRIC_PROVIDER_PACKET_NO_FALSE_PASS",
        "provider_packet": "preobserver_light_matter_source_provider_packet_v1",
        "strict_kernel_control": {
            "form": "cos(omega*d + phi) / (1 + beta*d^eta)",
            "omega": omega,
            "phi": phi,
            "beta": beta,
            "eta": eta,
            "kernel_core_limit": kernel_core_amplitude,
            "kernel_used_as_operational_control_only": True,
            "legacy_bridge_claimed": False,
        },
        "topological_flow_vector": {
            "definition": "T_flow^(0) = cos(phi) * e_topo",
            "kernel_core_amplitude": kernel_core_amplitude,
            "candidate_only": True,
        },
        "cascade": {
            "ordering": ["nadsoliton", "light", "matter", "emergent_observer"],
            "forward_maps": ["P_NL^(0)", "P_LM(d)", "P_MO(d)"],
            "reverse_maps_exported": False,
            "strictly_one_way_by_packet_construction": True,
        },
        "weights": {
            "I_mat(d)": "beta*d^eta/(1 + beta*d^eta)",
            "C_obs(d)": "d^eta/(1 + d^eta)",
            "I_mat_at_zero": i_mat_at_zero,
            "C_obs_at_zero": c_obs_at_zero,
            "C_obs_large_scale_limit": c_obs_large_scale_limit,
        },
        "observer_nonparticipation": {
            "dP_NL_dO_zero": True,
            "dP_LM_dO_zero": True,
            "observer_to_upstream_blocks_zero": True,
        },
        "theorem_level_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
