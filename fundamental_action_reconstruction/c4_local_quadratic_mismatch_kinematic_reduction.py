from __future__ import annotations

import json
import math
from pathlib import Path

root = Path(__file__).resolve().parent


def sample(theta: float, a: float, b: float) -> dict[str, float]:
    cos_t = math.cos(theta)
    du_norm_sq = 2.0 * (1.0 - cos_t)
    dv_norm_sq = 2.0 * (1.0 - cos_t)
    du_dv_inner = 0.0
    q_loc = a * du_norm_sq + b * dv_norm_sq
    closed = 2.0 * (a + b) * (1.0 - cos_t)
    return {
        "theta": theta,
        "a": a,
        "b": b,
        "du_norm_sq": du_norm_sq,
        "dv_norm_sq": dv_norm_sq,
        "du_dv_inner": du_dv_inner,
        "q_local": q_loc,
        "closed_form": closed,
    }


rows = [
    sample(0.0, 1.0, 1.0),
    sample(math.pi / 3.0, 1.0, 1.0),
    sample(math.pi / 2.0, 2.0, 1.0),
    sample(math.pi, 1.2, 2.8),
]

payload = {
    "stage": "C4",
    "status": "C4_EXECUTED_KINEMATIC_REDUCTION_OF_C2_B2_LOCAL_METRIC_IDENTIFICATION_PENDING_NO_FALSE_PASS",
    "as_of": "2026-03-06",
    "goal": "Reduce C2_B2 by showing that on the O(2) rotation orbit the J_ab family follows kinematically from any diagonal positive local quadratic mismatch metric.",
    "inputs": {
        "strict_admissible": [
            "QW-2190",
            "QW-2191",
            "QW-2192",
            "QW-2193",
            "A3",
            "A7",
            "C2",
            "C3",
            "A10",
        ]
    },
    "kinematic_identities": {
        "du_norm_sq": "2(1-cos(theta))",
        "dv_norm_sq": "2(1-cos(theta))",
        "du_dv_inner": "0 on the rotation orbit",
    },
    "conditional_result": {
        "assumption": "a positive diagonal local mismatch metric exists on the candidate orientation plane",
        "result": "Q_loc(theta)=a||Delta u||^2 + b||Delta v||^2 = 2(a+b)(1-cos(theta))",
        "status": "derived_conditionally"
    },
    "numeric_samples": rows,
    "frontier_reduction": {
        "C3_B1": "no_physical_elevation_of_deterministic_mode_pair_to_internal_orientation_datum",
        "C4_B1": "no_internal_identification_of_the_physical_positive_local_metric_on_candidate_orientation_plane"
    },
    "hard_limits": [
        "no_theorem_level_pass",
        "no_full_closure_pass",
        "no_c2_b2_pass",
        "no_internal_origin_of_metric_claim",
        "no_qw2191_discharge_claim"
    ],
    "next_step": "C5"
}

out = root / "generated" / "c4_local_quadratic_mismatch_kinematic_reduction_summary.json"
out.write_text(json.dumps(payload, indent=2) + "\n", encoding="ascii")
print(out)
