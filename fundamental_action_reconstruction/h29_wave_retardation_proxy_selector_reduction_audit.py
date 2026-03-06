#!/usr/bin/env python3
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
GEN.mkdir(exist_ok=True)
OUT = GEN / "h29_wave_retardation_proxy_selector_reduction_audit.json"
OUT_SUM = GEN / "h29_wave_retardation_proxy_selector_reduction_audit_summary.json"


def main() -> None:
    data = {
        "step": "H29",
        "title": "wave_retardation_proxy_selector_reduction_audit",
        "inputs": {
            "anisotropic_psf_term": "cos(2*(th-psi0)+phase)",
            "retardation_proxy": "retard_phase",
            "memory_proxy": "observer_tau",
            "feedback_gain_proxy": "observer_feedback_gain",
            "orientation_anchor": "orientation_psi0",
        },
        "reduction": {
            "expanded_angular_term": "cos(2*(th-psi0))*cos(phase) - sin(2*(th-psi0))*sin(phase)",
            "proxy_roles": {
                "retard_phase": "reweights_existing_angular_harmonics",
                "observer_tau": "sets_response_time_or_memory_scale",
                "observer_feedback_gain": "sets_feedback_amplitude",
                "orientation_psi0": "supplies_preoriented_directional_anchor",
            },
        },
        "selector_sector_result": {
            "wave_and_memory_proxies_alone_generate_anchor": False,
            "wave_and_memory_proxies_only_modulate_preanchored_channel": True,
            "strict_core_internal_theta_source_recovered": False,
        },
        "frontier": "H29_B1",
        "status": "PASS_PARTIAL_PREORIENTED_PROXY_REDUCTION_ONLY",
    }
    summary = {
        "step": "H29",
        "status": data["status"],
        "frontier": "H29_B1 := existing wave/retardation/memory proxies only modulate a pre-oriented anisotropic channel and do not by themselves generate an internal strict-core orientation anchor for theta_i",
        "no_theorem_level_pass": True,
        "no_full_closure_pass": True,
    }
    OUT.write_text(json.dumps(data, indent=2) + "\n")
    OUT_SUM.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
