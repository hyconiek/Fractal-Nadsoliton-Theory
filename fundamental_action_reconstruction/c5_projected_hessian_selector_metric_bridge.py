from __future__ import annotations

import json
import math
from pathlib import Path

root = Path(__file__).resolve().parent


def sample(theta: float, a: float, b: float, c: float) -> dict[str, float]:
    cos_t = math.cos(theta)
    du_du = 2.0 * (1.0 - cos_t)
    dv_dv = 2.0 * (1.0 - cos_t)
    du_dv = 0.0
    q_h = a * du_du + 2.0 * c * du_dv + b * dv_dv
    closed = 2.0 * (a + b) * (1.0 - cos_t)
    return {
        "theta": theta,
        "a": a,
        "b": b,
        "c": c,
        "du_du": du_du,
        "dv_dv": dv_dv,
        "du_dv": du_dv,
        "q_h": q_h,
        "closed_form": closed,
    }


rows = [
    sample(0.0, 1.0, 1.0, 0.2),
    sample(math.pi / 3.0, 1.0, 1.0, 0.2),
    sample(math.pi / 2.0, 2.0, 1.0, 0.4),
    sample(math.pi, 1.2, 2.8, 0.5),
]

payload = {
    "stage": "C5",
    "status": "C5_EXECUTED_PROJECTED_HESSIAN_BRIDGE_CONDITIONAL_LOCAL_IDENTIFICATION_PENDING_NO_FALSE_PASS",
    "as_of": "2026-03-06",
    "goal": "Test whether the selector-family mismatch metric can be interpreted as the orbital form of a projected local second variation / Hessian.",
    "inputs": {
        "strict_admissible": [
            "A3",
            "A7",
            "QW-2190",
            "QW-2191",
            "C3",
            "C4",
            "A10",
        ]
    },
    "conditional_bridge": {
        "assumption": "an explicit projected second variation exists on the candidate orientation plane and has a local symmetric quadratic form with a strict-scope positivity certificate",
        "quadratic_form": "Q_H(theta)=a<Delta u,Delta u> + 2c<Delta u,Delta v> + b<Delta v,Delta v>",
        "orbital_reduction": "Q_H(theta)=2(a+b)(1-cos(theta))",
        "status": "derived_conditionally"
    },
    "numeric_samples": rows,
    "frontier_reduction": {
        "C3_B1": "no_physical_elevation_of_deterministic_mode_pair_to_internal_orientation_datum",
        "C5_B1": "no_explicit_projected_second_variation_with_strict_scope_positivity_certificate_on_candidate_orientation_plane"
    },
    "hard_limits": [
        "no_theorem_level_pass",
        "no_full_closure_pass",
        "no_c2_b2_pass",
        "no_projected_hessian_found_claim",
        "no_qw2191_discharge_claim"
    ],
    "next_step": "C6"
}

out = root / "generated" / "c5_projected_hessian_selector_metric_bridge_summary.json"
out.write_text(json.dumps(payload, indent=2) + "\n", encoding="ascii")
print(out)
