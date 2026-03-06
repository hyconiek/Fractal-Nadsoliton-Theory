from __future__ import annotations

import json
import math
from pathlib import Path

root = Path(__file__).resolve().parent


def local_projectors(theta: float) -> dict[str, object]:
    c = math.cos(theta)
    s = math.sin(theta)
    p_tan = [[s * s, -s * c], [-s * c, c * c]]
    p_red = [[c * c, c * s], [c * s, s * s]]
    return {
        "theta": theta,
        "P_tan": p_tan,
        "P_red": p_red,
        "trace_P_tan": p_tan[0][0] + p_tan[1][1],
        "trace_P_red": p_red[0][0] + p_red[1][1],
    }


samples = [
    local_projectors(0.0),
    local_projectors(math.pi / 4.0),
    local_projectors(math.pi / 2.0),
]

payload = {
    "stage": "C29",
    "status": "C29_EXECUTED_LOCAL_PROJECTOR_FORMULA_PACKET_NO_FALSE_PASS",
    "as_of": "2026-03-06",
    "goal": "Reduce C28_B1 by checking whether the strict core already supports an explicit serialized local projector formula on each mode-pair plane, even if no global pair-to-pair gluing rule is exported yet.",
    "inputs": {
        "strict_admissible": [
            "C4",
            "C28",
            "C14",
            "C15",
            "A10",
        ]
    },
    "local_frame": {
        "e_theta": "(cos(theta), sin(theta))",
        "tau_theta": "(-sin(theta), cos(theta))",
    },
    "serialized_formulas": {
        "P_tan(theta)": "[[sin^2(theta), -sin(theta)cos(theta)], [-sin(theta)cos(theta), cos^2(theta)]]",
        "P_red(theta)": "[[cos^2(theta), cos(theta)sin(theta)], [cos(theta)sin(theta), sin^2(theta)]]",
        "relations": [
            "P_tan + P_red = I_2",
            "P_tan^2 = P_tan",
            "P_red^2 = P_red",
            "P_tan P_red = 0",
        ],
    },
    "numeric_samples": samples,
    "result": {
        "explicit_local_projector_formula_present": "yes_partial",
        "global_pair_to_pair_gluing_rule_present": "not_shown",
        "final_basis_level_slice_extraction_present": "not_shown",
    },
    "residual_blockers": {
        "C29_B1": "no_explicit_pair_to_pair_global_gluing_rule_assembling_the_local_reduced_lines_into_a_single_reduced_control_plane",
        "C26_B2": "no_explicit_basis_level_embedding_or_extraction_of_the_candidate_two_dimensional_orientation_slice_inside_that_reduced_plane",
    },
    "hard_limits": [
        "no_theorem_level_pass",
        "no_full_closure_pass",
        "no_claim_that_global_gluing_is_resolved",
        "no_claim_that_final_slice_extraction_is_resolved",
        "no_claim_that_selector_track_is_closed",
    ],
    "next_step": "C30",
}

out = root / "generated" / "c29_local_projector_formula_packet_summary.json"
out.write_text(json.dumps(payload, indent=2) + "\n", encoding="ascii")
print(out)
