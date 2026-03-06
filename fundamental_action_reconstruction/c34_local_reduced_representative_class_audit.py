#!/usr/bin/env python3
"""C34: local reduced representative class audit without false PASS."""

from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "c34_local_reduced_representative_class_audit_summary.json"


def mv(m: list[list[float]], v: list[float]) -> list[float]:
    return [m[0][0] * v[0] + m[0][1] * v[1], m[1][0] * v[0] + m[1][1] * v[1]]


def main() -> None:
    samples = []
    for theta in [0.0, math.pi / 6, math.pi / 3, 7 * math.pi / 10]:
        c = math.cos(theta)
        s = math.sin(theta)
        u = [c, s]
        p_red = [[c * c, c * s], [c * s, s * s]]
        p_tan = [[s * s, -s * c], [-s * c, c * c]]
        red_u = mv(p_red, u)
        tan_u = mv(p_tan, u)
        samples.append(
            {
                "theta": theta,
                "u_i": [round(u[0], 12), round(u[1], 12)],
                "norm_sq": round(u[0] * u[0] + u[1] * u[1], 12),
                "P_red_u": [round(red_u[0], 12), round(red_u[1], 12)],
                "P_tan_u": [round(tan_u[0], 12), round(tan_u[1], 12)],
            }
        )

    summary = {
        "stage": "C34",
        "status": "C34_EXECUTED_LOCAL_REDUCED_REPRESENTATIVE_CLASS_AUDIT_NO_FALSE_PASS",
        "as_of": "2026-03-06",
        "goal": "Reduce C33_B1 by checking whether strict core already supports a packet-ready explicit normalized representative class on a single local reduced line, even if actual theta_1, theta_2 are not exported for the current pair frames.",
        "inputs": {"strict_admissible": ["C4", "C28", "C29", "C33", "A10"]},
        "representative_class": {
            "formula": "u_i(theta_i) = cos(theta_i) c_i + sin(theta_i) s_i",
            "normalization": "||u_i|| = 1",
            "projector_relations": [
                "P_red(theta_i) u_i = u_i",
                "P_tan(theta_i) u_i = 0"
            ]
        },
        "numeric_samples": samples,
        "result": {
            "local_reduced_representative_class_present": "yes_partial",
            "projector_compatibility_present": "yes_partial",
            "explicit_exported_u_1_u_2_for_actual_pair_frames": "not_shown",
            "explicit_exported_theta_1_theta_2_for_actual_pair_frames": "not_shown",
            "raw_cross_pair_overlap_scalar_route": "blocked_by_degeneracy",
            "final_basis_level_slice_extraction_present": "not_shown"
        },
        "residual_blockers": {
            "C34_B1": "no_explicit_export_of_actual_local_phase_coordinates_theta_1_theta_2_needed_to_materialize_the_normalized_local_reduced_representatives_u_1_u_2_for_the_actual_pair_frames",
            "C32_B2": "raw_cross_pair_overlap_scalar_route_is_formally_degenerate_under_the_strict_orthonormal_disjoint_mode_scaffold_and_thus_does_not_export_alpha_12",
            "C26_B2": "no_explicit_basis_level_embedding_or_extraction_of_the_candidate_two_dimensional_orientation_slice_inside_that_reduced_plane"
        },
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_theta_1_theta_2_are_exported",
            "no_claim_that_u_1_u_2_are_materialized_for_actual_pair_frames",
            "no_claim_that_alpha_12_is_exported",
            "no_claim_that_final_slice_extraction_is_resolved"
        ],
        "next_step": "C35"
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
