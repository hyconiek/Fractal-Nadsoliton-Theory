#!/usr/bin/env python3
"""C31: transition-angle source class audit without false PASS."""

from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "c31_transition_angle_source_candidate_audit_summary.json"

IN_THETA_PAIR = GENERATED / "theta_pair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1.json"
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

    theta_exported = IN_THETA_PAIR.exists()
    theta_1_actual = None
    theta_2_actual = None
    if theta_exported:
        try:
            obj = load_json(IN_THETA_PAIR)
            theta_1_actual = float(((obj.get("outputs") or {}).get("pair1") or {}).get("theta_1"))
            theta_2_actual = float(((obj.get("outputs") or {}).get("pair2") or {}).get("theta_2"))
        except Exception:
            theta_exported = False
            theta_1_actual = None
            theta_2_actual = None

    alpha12_exported = IN_ALPHA12.exists()
    alpha12_mod_2pi = None
    alpha12_mod_pi = None
    if alpha12_exported:
        try:
            a = load_json(IN_ALPHA12)
            alpha12_mod_2pi = float(((a.get("outputs") or {}).get("alpha_12_mod_2pi")))
            alpha12_mod_pi = float(((a.get("outputs") or {}).get("alpha_12_mod_pi")))
        except Exception:
            alpha12_exported = False
            alpha12_mod_2pi = None
            alpha12_mod_pi = None

    summary = {
        "stage": "C31",
        "status": "C31_EXECUTED_TRANSITION_ANGLE_SOURCE_CANDIDATE_AUDIT_NO_FALSE_PASS",
        "as_of": "2026-03-15",
        "goal": "Reduce C30_B1 by checking whether strict core already supports a packet-ready source class for the pair-to-pair transition angle, even if no explicit exported value is present for the actual pair frames.",
        "inputs": {
            "strict_admissible": ["C4", "C29", "C30", "C3", "C13", "A10"]
        },
        "source_class": {
            "candidate": "alpha_12 = theta_2 - theta_1 (mod 2pi)",
            "interpretation": "relative phase difference between the two local orbit coordinates"
        },
        "numeric_samples": samples,
        "current_repo_state_if_present": {
            "theta_pair_export_ref": str(IN_THETA_PAIR.relative_to(ROOT)) if theta_exported else None,
            "theta_1": theta_1_actual,
            "theta_2": theta_2_actual,
            "alpha12_export_ref": str(IN_ALPHA12.relative_to(ROOT)) if alpha12_exported else None,
            "alpha_12_mod_2pi": alpha12_mod_2pi,
            "alpha_12_mod_pi": alpha12_mod_pi,
        },
        "result": {
            "transition_angle_source_class_present": "yes_partial",
            "explicit_theta_1_theta_2_export_for_actual_pair_frames": (
                "exported_in_declared_sigma_int_scope (F451)" if theta_exported else "not_shown"
            ),
            "explicit_overlap_scalar_for_actual_pair_frames": "not_shown",
            "explicit_serialized_alpha_12_for_actual_pair_frames": (
                "exported_in_declared_sigma_int_scope (F457)"
                if alpha12_exported
                else "not_shown"
            ),
            "final_basis_level_slice_extraction_present": "not_shown"
        },
        "residual_blockers": {
            "C31_B1": (
                "no_overlap_scalar_route_for_alpha12_under_disjoint_orthonormal_mode_scaffold; "
                "alpha12_is_exported_only_via_theta_supply_in_declared_sigma_int_scope"
            ),
            "C26_B2": "no_explicit_basis_level_embedding_or_extraction_of_the_candidate_two_dimensional_orientation_slice_inside_that_reduced_plane"
        },
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_alpha_12_is_exported_globally_beyond_declared_sigma_int_scope",
            "no_claim_that_global_gluing_is_resolved",
            "no_claim_that_final_slice_extraction_is_resolved"
        ],
        "next_step": "C32"
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
