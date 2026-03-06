#!/usr/bin/env python3
"""C33: local phase export formula class audit without false PASS."""

from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "c33_local_phase_export_class_audit_summary.json"


def main() -> None:
    samples = []
    for theta in [0.0, math.pi / 6, math.pi / 2, 7 * math.pi / 10]:
        x = math.cos(theta)
        y = math.sin(theta)
        recovered = math.atan2(y, x)
        samples.append(
            {
                "theta": theta,
                "x=<c_i,u_i>": round(x, 12),
                "y=<s_i,u_i>": round(y, 12),
                "atan2(y,x)": round(recovered, 12),
            }
        )

    summary = {
        "stage": "C33",
        "status": "C33_EXECUTED_LOCAL_PHASE_EXPORT_CLASS_AUDIT_NO_FALSE_PASS",
        "as_of": "2026-03-06",
        "goal": "Reduce C32_B1 by checking whether strict core already supports a packet-ready formula class for exporting the local phase theta_i on a single mode pair, even if no explicit representative u_i is exported for the actual pair frames.",
        "inputs": {
            "strict_admissible": ["C4", "C28", "C29", "C31", "A10"]
        },
        "formula_class": {
            "representative": "u_i = x_i c_i + y_i s_i with ||u_i||=1",
            "coordinates": "x_i=<c_i,u_i>, y_i=<s_i,u_i>",
            "phase_export": "theta_i = atan2(<s_i,u_i>, <c_i,u_i>)"
        },
        "numeric_samples": samples,
        "result": {
            "local_phase_export_formula_class_present": "yes_partial",
            "explicit_representative_u_i_for_actual_pair_frames": "not_shown",
            "explicit_exported_theta_1_theta_2_for_actual_pair_frames": "not_shown",
            "raw_cross_pair_overlap_scalar_route": "blocked_by_degeneracy",
            "final_basis_level_slice_extraction_present": "not_shown"
        },
        "residual_blockers": {
            "C33_B1": "no_explicit_export_of_normalized_local_reduced_representatives_u_1_u_2_for_the_actual_pair_frames_from_which_theta_1_theta_2_could_be_serialized_via_atan2",
            "C32_B2": "raw_cross_pair_overlap_scalar_route_is_formally_degenerate_under_the_strict_orthonormal_disjoint_mode_scaffold_and_thus_does_not_export_alpha_12",
            "C26_B2": "no_explicit_basis_level_embedding_or_extraction_of_the_candidate_two_dimensional_orientation_slice_inside_that_reduced_plane"
        },
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_u_1_u_2_are_exported",
            "no_claim_that_theta_1_theta_2_are_exported",
            "no_claim_that_alpha_12_is_exported",
            "no_claim_that_final_slice_extraction_is_resolved"
        ],
        "next_step": "C34"
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
