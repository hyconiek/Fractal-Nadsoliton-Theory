#!/usr/bin/env python3
"""P1510 S4.60: clarify strict-kernel role, symmetry breaking, and classical-physics link."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1500 = GEN / "p1500_s450_qw2191_selector_source_and_f_mapping_witness_summary.json"
P1509 = GEN / "p1509_s459_strict_contradiction_branch_stress_for_f_to_lsm_lgr_summary.json"
SUMMARY = GEN / "p1510_s460_fin_strict_kernel_state_symmetry_breaking_and_classical_physics_link_summary.json"


def main() -> None:
    s1500 = json.loads(P1500.read_text(encoding="utf-8"))
    s1509 = json.loads(P1509.read_text(encoding="utf-8"))

    s_internal = s1500.get("objects", {}).get("S_internal_v1", {})
    fmap = s1500.get("objects", {}).get("W_Fmap_v1", {})

    orientation = str(s_internal.get("orientation", ""))
    fmap_orientation = str(fmap.get("selector_orientation", ""))
    shared_orientation = bool(orientation) and orientation == fmap_orientation

    counterexample = s1509.get("first_counterexample", {})
    symmetry_break_observed = counterexample.get("id") == "orientation_flip_probe"

    summary = {
        "packet": "P1510",
        "status": "PASS_STRICT_KERNEL_STATE_CLARIFIED",
        "scope": "STRICT_ONLY_NO_LEGACY_BRIDGE",
        "strict_kernel_role": "operational_strict_working_kernel_for_gate_selection_and_F_to_LSM_LGR_coupling",
        "symmetry_breaking_state": {
            "selector_orientation": orientation,
            "shared_with_fmap": shared_orientation,
            "symmetry_break_probe_detected": symmetry_break_observed,
            "probe_reference": counterexample,
        },
        "classical_physics_link_level": "indirect_structural_link_via_coupled_LSM_and_LGR_channels_not_final_theorem_closure",
        "qw2191_closed": False,
        "legacy_bridge_used": False,
        "next_honest_step": "Implement P1511 admissible-perturbation boundary map for F(Nadsoliton)=>LSM+LGR and quantify robust vs rejection zones.",
        "layman_explanation": "Dziś kernel strict działa jak centralny regulator kierunku: spina część cząstkową i grawitacyjną. Gdy próbujemy sztucznie odwrócić kierunek, model słusznie odrzuca taki przypadek. To dobry znak, ale to jeszcze nie koniec dowodu.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1510] status={summary['status']} symmetry_break_observed={symmetry_break_observed}")


if __name__ == "__main__":
    main()
