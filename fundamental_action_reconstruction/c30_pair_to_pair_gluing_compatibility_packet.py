#!/usr/bin/env python3
"""C30: local pair-to-pair gluing compatibility packet without false PASS."""

from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "c30_pair_to_pair_gluing_compatibility_packet_summary.json"

IN_ALPHA12 = (
    GENERATED
    / "alpha12_pair1_pair2_transition_angle_strict_derived_from_sigma_int_slot_free_theta_pair_v1.json"
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def rot(alpha: float) -> list[list[float]]:
    c = math.cos(alpha)
    s = math.sin(alpha)
    return [[c, -s], [s, c]]


def e(theta: float) -> list[float]:
    return [math.cos(theta), math.sin(theta)]


def p_red(theta: float) -> list[list[float]]:
    c = math.cos(theta)
    s = math.sin(theta)
    return [[c * c, c * s], [c * s, s * s]]


def p_tan(theta: float) -> list[list[float]]:
    c = math.cos(theta)
    s = math.sin(theta)
    return [[s * s, -s * c], [-s * c, c * c]]


def mm(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [
        [sum(a[i][k] * b[k][j] for k in range(2)) for j in range(2)]
        for i in range(2)
    ]


def mt(a: list[list[float]]) -> list[list[float]]:
    return [[a[j][i] for j in range(2)] for i in range(2)]


def rounded(m: list[list[float]]) -> list[list[float]]:
    return [[round(x, 12) for x in row] for row in m]


def main() -> None:
    samples = []
    for theta, alpha in [(0.0, 0.0), (math.pi / 6, math.pi / 4), (math.pi / 5, math.pi / 3)]:
        g = rot(alpha)
        pred_push = mm(mm(g, p_red(theta)), mt(g))
        ptan_push = mm(mm(g, p_tan(theta)), mt(g))
        samples.append(
            {
                "theta": theta,
                "alpha": alpha,
                "P_red(theta)": rounded(p_red(theta)),
                "G(alpha)": rounded(g),
                "G P_red(theta) G^T": rounded(pred_push),
                "P_red(theta+alpha)": rounded(p_red(theta + alpha)),
                "G P_tan(theta) G^T": rounded(ptan_push),
                "P_tan(theta+alpha)": rounded(p_tan(theta + alpha)),
            }
        )

    alpha12_exported = IN_ALPHA12.exists()
    alpha12_mod_2pi = None
    if alpha12_exported:
        try:
            a = load_json(IN_ALPHA12)
            alpha12_mod_2pi = float(((a.get("outputs") or {}).get("alpha_12_mod_2pi")))
        except Exception:
            alpha12_exported = False
            alpha12_mod_2pi = None

    summary = {
        "stage": "C30",
        "status": "C30_EXECUTED_PAIR_TO_PAIR_GLUING_COMPATIBILITY_PACKET_NO_FALSE_PASS",
        "as_of": "2026-03-15",
        "goal": "Reduce C29_B1 by checking whether strict core already supports a packet-ready pair-to-pair overlap compatibility law for local reduced lines under orthogonal transition, even if no explicit exported transition matrix is present.",
        "inputs": {
            "strict_admissible": ["C4", "C28", "C29", "C14", "A10"]
        },
        "compatibility_law": {
            "transition": "G(alpha) in O(2)",
            "vector_rule": "G(alpha)e(theta) = e(theta+alpha)",
            "reduced_projector_rule": "G(alpha) P_red(theta) G(alpha)^T = P_red(theta+alpha)",
            "tangent_projector_rule": "G(alpha) P_tan(theta) G(alpha)^T = P_tan(theta+alpha)",
        },
        "numeric_samples": samples,
        "result": {
            "local_pair_to_pair_overlap_compatibility_present": "yes_partial",
            "explicit_serialized_transition_matrix_between_actual_pair1_and_pair2": "not_shown",
            "explicit_transition_angle_between_actual_pair1_and_pair2": (
                "exported_in_declared_sigma_int_scope (F457)"
                if alpha12_exported
                else "not_shown"
            ),
            "alpha_12_mod_2pi_if_exported": alpha12_mod_2pi,
            "final_basis_level_slice_extraction_present": "not_shown",
        },
        "residual_blockers": {
            "C30_B1": (
                "no_global_pair_to_pair_gluing_object_lifting_local_compatibility_to_global_transport; "
                "lane_scoped_transition_angle_export_exists_via_F457_but_no_global_transport_between_disjoint_pair_frames_is_exported"
            ),
            "C26_B2": "no_explicit_basis_level_embedding_or_extraction_of_the_candidate_two_dimensional_orientation_slice_inside_that_reduced_plane",
        },
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_global_gluing_is_resolved",
            "no_claim_that_transition_matrix_G12_is_exported",
            "no_claim_that_final_slice_extraction_is_resolved",
        ],
        "next_step": "C31",
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
