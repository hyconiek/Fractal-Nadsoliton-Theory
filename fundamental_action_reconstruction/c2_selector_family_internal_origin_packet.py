from __future__ import annotations

import json
import math
from pathlib import Path

root = Path(__file__).resolve().parent

samples = []
for a, b in [(1.0, 1.0), (2.0, 1.0), (1.2, 2.8)]:
    theta0 = 0.0
    thetapi = math.pi
    j0 = 2.0 * (a + b) * (1.0 - math.cos(theta0))
    jpi = 2.0 * (a + b) * (1.0 - math.cos(thetapi))
    second0 = 2.0 * (a + b)
    samples.append(
        {
            "a": a,
            "b": b,
            "J_theta0": j0,
            "J_theta_pi": jpi,
            "second_derivative_at_0": second0,
            "argmin_mod_2pi": 0.0,
        }
    )

payload = {
    "stage": "C2",
    "status": "C2_EXECUTED_CONDITIONAL_ORIGIN_REDUCTION_PACKET_READY_NO_FALSE_PASS",
    "as_of": "2026-03-06",
    "goal": "Reduce the selector-family origin blocker by deriving the J_ab form conditionally from a reference pair and a positive quadratic mismatch principle.",
    "assumptions": {
        "C2_A1": "internal reference pair exists in the degenerate two-mode plane",
        "C2_A2": "selector cost is local, positive, and quadratic: J_ab(theta)=a||u-c_ref||^2 + b||v-s_ref||^2 with a>0, b>0"
    },
    "conditional_result": {
        "closed_form": "J_ab(theta)=2(a+b)(1-cos(theta))",
        "theta_star": "0 mod 2pi",
        "status": "derived_conditionally"
    },
    "numeric_samples": samples,
    "reduced_blockers": {
        "C2_B1": "no_derived_internal_reference_pair_for_degenerate_mode_plane",
        "C2_B2": "no_derived_positive_local_quadratic_mismatch_principle"
    },
    "hard_limits": [
        "no_theorem_level_pass",
        "no_full_closure_pass",
        "no_internal_origin_of_Jab_claim",
        "no_axiom_free_uniqueness_claim"
    ],
    "next_step": "C3"
}

out = root / "generated" / "c2_selector_family_internal_origin_packet_summary.json"
out.write_text(json.dumps(payload, indent=2) + "\n", encoding="ascii")
print(out)
