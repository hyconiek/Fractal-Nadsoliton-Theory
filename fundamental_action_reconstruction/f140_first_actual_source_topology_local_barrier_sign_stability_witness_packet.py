#!/usr/bin/env python3

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
IN_F74 = GENERATED / "f74_first_preobserver_light_matter_asymmetric_provider_packet_summary.json"
IN_F127 = GENERATED / "f127_first_source_topology_invariant_candidate_packet_summary.json"
IN_F139 = GENERATED / "f139_first_actual_source_topology_barrier_sign_component_witness_packet_summary.json"
OUT = GENERATED / "f140_first_actual_source_topology_local_barrier_sign_stability_witness_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f74 = load_json(IN_F74)
    f127 = load_json(IN_F127)
    f139 = load_json(IN_F139)

    phi = f74["strict_kernel_control"]["phi"]
    barrier_margin = f139["barrier_margin_value"]
    local_radius = barrier_margin / 2.0
    abs_phi_plus_radius = abs(phi) + local_radius
    interval_inside_positive_cos_domain = abs_phi_plus_radius < (math.pi / 2.0)
    lower_endpoint = phi - local_radius
    upper_endpoint = phi + local_radius

    summary = {
        "packet_id": "F140",
        "status": "F140_EXECUTED_FIRST_ACTUAL_SOURCE_TOPOLOGY_LOCAL_BARRIER_SIGN_STABILITY_WITNESS_PACKET_NO_FALSE_PASS",
        "as_of": "2026-03-08",
        "input_packet": f127["candidate_packet"],
        "witness": "chi_src_local_barrier_sign_stability_witness_v1",
        "witness_definition": (
            "for all epsilon in R, if |epsilon| <= epsilon_src_local_barrier_radius_v1, "
            "then sign(cos(phi + epsilon)) = +1"
        ),
        "local_radius": "epsilon_src_local_barrier_radius_v1",
        "local_radius_definition": "delta_src_barrier_sign_margin_v1 / 2",
        "local_radius_value": local_radius,
        "local_radius_positive": local_radius > 0.0,
        "phi_value": phi,
        "barrier_margin_value": barrier_margin,
        "scalar_sign_component_from_f139": f139["witness_value"],
        "abs_phi_plus_radius": abs_phi_plus_radius,
        "interval_inside_positive_cos_domain": interval_inside_positive_cos_domain,
        "local_interval_lower_endpoint": lower_endpoint,
        "local_interval_upper_endpoint": upper_endpoint,
        "source_limit_tag": f127["components"]["source_limit_tag_v1"],
        "phi_barrier_tag": f127["components"]["phi_barrier_tag_v1"],
        "observer_role": "downstream_only",
        "actual_local_sign_stability_witness_exported": (
            f139["witness_value"] == 1
            and local_radius > 0.0
            and interval_inside_positive_cos_domain is True
        ),
        "barrier_protected_sign_discharged": False,
        "full_source_topology_nontriviality_discharged": False,
        "basis_independence_discharged": False,
        "qw2191_quotient_safe_discharged": False,
        "current_selector_closure": False,
        "current_global_qw2191_discharge": False,
        "kernel_split_safe": True,
        "legacy_to_strict_bridge_claimed": False,
        "no_false_pass": True,
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
