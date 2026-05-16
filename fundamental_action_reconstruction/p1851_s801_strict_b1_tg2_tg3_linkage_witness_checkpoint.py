#!/usr/bin/env python3
"""P1851 S801 strict B1 TG2/TG3 linkage witness checkpoint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    p1850 = load("p1850_s800_strict_gravity_background_b1_symbolic_coefficients_checkpoint.json")
    p1839 = load("p1839_s789_strict_full_lagrangian_term_delivery_contract_checkpoint.json")
    p1844 = load("p1844_s794_strict_toe_qg_closure_blocker_matrix_checkpoint.json")

    tg2_contract = {
        "gate": "TG2_BRST",
        "required_witnesses": [
            "brst_nilpotency_nonproxy_proof",
            "cohomology_class_stability_under_b1_renormalization",
            "ward_identity_compatibility_with_delta_c_gr_i",
        ],
        "b1_link": "counterterm_cancellation_identity_b1 must be BRST-exact compatible (no anomaly term on declared B1 family)",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    tg3_contract = {
        "gate": "TG3_CUT",
        "required_witnesses": [
            "cutkosky_discontinuity_full_sector_theorem",
            "positive_spectral_density_for_dressed_graviton_channels",
            "background_family_B1_to_background_independence_lift_theorem",
        ],
        "b1_link": "B1 renormalized couplings and residues must preserve optical-theorem discontinuity relations",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    update_matrix = {
        "renormalization_counterterm_closure": "PARTIAL_TRACE_ONLY",
        "brst_nilpotency": "OPEN_OBSTRUCTION_WITH_TRACE",
        "cutkosky_unitarity": "OPEN_OBSTRUCTION_WITH_TRACE",
        "background_independence": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    out = {
        "packet_id": "P1851",
        "stage_id": "S801",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1850_present": "counterterm_cancellation_identity_b1" in p1850,
            "p1839_present": "theorem_gate_dependencies" in p1839,
            "p1844_present": "toe_qg_blocker_matrix" in p1844,
        },
        "b1_forward_reverse_chain": "K_strict -> symbolic_a_i(B1) -> L_total counterterm dressing -> EOM residual discipline -> TG2_BRST/TG3_CUT linkage",
        "tg2_brst_link_contract": tg2_contract,
        "tg3_cut_background_link_contract": tg3_contract,
        "qg_blocker_matrix_update_proposal": update_matrix,
        "proven": "B1 renormalization trace is now explicitly linked to TG2/TG3 witness obligations in a strict forward-reverse contract.",
        "open": "TG2/TG3 theorem witnesses remain missing; background-independence lift beyond B1 is still open.",
        "false_pass_risk": "Having B1 symbolic counterterm identities does not imply BRST nilpotency, Cutkosky unitarity, or global background-independence closure.",
        "next_honest_step": "Export first BRST anomaly-check witness on B1 and first Cutkosky discontinuity witness for dressed graviton-matter channels in the same renormalization scheme.",
        "lay_explanation": "Połączyliśmy już matematykę kasowania nieskończoności z kolejnymi testami fizycznymi (unitarność i symetrie), ale dalej brakuje twardych dowodów, więc ToE nie jest jeszcze domknięte.",
    }

    path = GEN / "p1851_s801_strict_b1_tg2_tg3_linkage_witness_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
