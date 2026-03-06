#!/usr/bin/env python3
"""C31: transition-angle source class audit without false PASS."""

from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "c31_transition_angle_source_candidate_audit_summary.json"


def rot(alpha: float) -> list[list[float]]:
    c = math.cos(alpha)
    s = math.sin(alpha)
    return [[c, -s], [s, c]]


def mv(a: list[list[float]], v: list[float]) -> list[float]:
    return [sum(a[i][k] * v[k] for k in range(2)) for i in range(2)]


def e(theta: float) -> list[float]:
    return [math.cos(theta), math.sin(theta)]


def p_red(theta: float) -> list[list[float]]:
    c = math.cos(theta)
    s = math.sin(theta)
    return [[c * c, c * s], [c * s, s * s]]


def mm(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [
        [sum(a[i][k] * b[k][j] for k in range(2)) for j in range(2)]
        for i in range(2)
    ]


def mt(a: list[list[float]]) -> list[list[float]]:
    return [[a[j][i] for j in range(2)] for i in range(2)]


def rounded_vec(v: list[float]) -> list[float]:
    return [round(x, 12) for x in v]


def rounded_mat(m: list[list[float]]) -> list[list[float]]:
    return [[round(x, 12) for x in row] for row in m]


def main() -> None:
    samples = []
    for theta1, theta2 in [(0.0, 0.0), (math.pi / 6, math.pi / 2), (math.pi / 5, 7 * math.pi / 10)]:
        alpha12 = theta2 - theta1
        g12 = rot(alpha12)
        lhs_vec = mv(g12, e(theta1))
        lhs_proj = mm(mm(g12, p_red(theta1)), mt(g12))
        samples.append(
            {
                "theta_1": theta1,
                "theta_2": theta2,
                "alpha_12": alpha12,
                "G(alpha_12) e(theta_1)": rounded_vec(lhs_vec),
                "e(theta_2)": rounded_vec(e(theta2)),
                "G(alpha_12) P_red(theta_1) G(alpha_12)^T": rounded_mat(lhs_proj),
                "P_red(theta_2)": rounded_mat(p_red(theta2)),
            }
        )

    summary = {
        "stage": "C31",
        "status": "C31_EXECUTED_TRANSITION_ANGLE_SOURCE_CANDIDATE_AUDIT_NO_FALSE_PASS",
        "as_of": "2026-03-06",
        "goal": "Reduce C30_B1 by checking whether strict core already supports a packet-ready source class for the pair-to-pair transition angle, even if no explicit exported value is present for the actual pair frames.",
        "inputs": {
            "strict_admissible": ["C4", "C29", "C30", "C3", "C13", "A10"]
        },
        "source_class": {
            "candidate": "alpha_12 = theta_2 - theta_1 (mod 2pi)",
            "interpretation": "relative phase difference between the two local orbit coordinates"
        },
        "numeric_samples": samples,
        "result": {
            "transition_angle_source_class_present": "yes_partial",
            "explicit_theta_1_theta_2_export_for_actual_pair_frames": "not_shown",
            "explicit_overlap_scalar_for_actual_pair_frames": "not_shown",
            "explicit_serialized_alpha_12_for_actual_pair_frames": "not_shown",
            "final_basis_level_slice_extraction_present": "not_shown"
        },
        "residual_blockers": {
            "C31_B1": "no_explicit_export_of_local_phase_coordinates_theta_1_theta_2_or_equivalent_pair_overlap_scalar_for_serializing_alpha_12_between_the_two_local_pair_frames",
            "C26_B2": "no_explicit_basis_level_embedding_or_extraction_of_the_candidate_two_dimensional_orientation_slice_inside_that_reduced_plane"
        },
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_alpha_12_is_exported_for_actual_pair_frames",
            "no_claim_that_global_gluing_is_resolved",
            "no_claim_that_final_slice_extraction_is_resolved"
        ],
        "next_step": "C32"
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
