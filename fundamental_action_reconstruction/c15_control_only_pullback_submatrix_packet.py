from __future__ import annotations

import json
from pathlib import Path

root = Path(__file__).resolve().parent

payload = {
    "stage": "C15",
    "status": "C15_EXECUTED_CONTROL_ONLY_PULLBACK_PACKET_NO_FALSE_PASS",
    "as_of": "2026-03-06",
    "goal": "Reduce the missing assembled submatrix blocker by showing that the strict core already admits a formal control-only pullback packet M_control = T_control^T H_PsiPsi T_control, even if no coefficient-filled canonical Psi-Psi block is exported yet.",
    "inputs": {
        "strict_admissible": [
            "QW-2190",
            "QW-2164",
            "QW-2166",
            "QW-2180",
            "C12",
            "C14",
            "A10",
        ]
    },
    "formal_objects": {
        "T_control_shape": [12, 4],
        "H_PsiPsi_shape": [12, 12],
        "M_control_shape": [4, 4],
        "control_basis": ["c1", "s1", "c2", "s2"],
        "assembly_formula": "M_control = T_control^T H_PsiPsi T_control",
    },
    "result": {
        "control_only_pullback_formula_present": "yes",
        "coefficient_filled_H_PsiPsi_present": "not_shown",
        "coefficient_filled_M_control_present": "not_shown",
        "orientation_slice_restriction_present": "not_shown",
    },
    "residual_blockers": {
        "C15_B1": "no_explicit_coefficient_filled_canonical_Psi_x_Psi_block_H_PsiPsi_for_evaluating_the_control_pullback",
        "C15_B2": "no_explicit_restriction_from_M_control_to_the_candidate_orientation_slice",
        "C9_B2": "no_explicit_restriction_from_that_Psi_sector_quadratic_carrier_to_the_candidate_orientation_slice",
    },
    "hard_limits": [
        "no_theorem_level_pass",
        "no_full_closure_pass",
        "no_coefficient_filled_h_psi_psi_claim",
        "no_coefficient_filled_m_control_claim",
        "no_orientation_slice_restriction_claim",
        "no_qw2191_discharge_claim",
    ],
    "next_step": "C16",
}

out = root / "generated" / "c15_control_only_pullback_submatrix_packet_summary.json"
out.write_text(json.dumps(payload, indent=2) + "\n", encoding="ascii")
print(out)
