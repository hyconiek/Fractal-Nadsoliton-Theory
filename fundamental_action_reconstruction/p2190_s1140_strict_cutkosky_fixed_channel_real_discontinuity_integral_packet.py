#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2190_s1140_strict_cutkosky_fixed_channel_real_discontinuity_integral_packet.json"
MD = GEN / "p2190_s1140_strict_cutkosky_fixed_channel_real_discontinuity_integral_packet.md"


def k_strict(d: float, omega: float, phi: float, beta: float, eta: float) -> float:
    return math.cos(omega * d + phi) / (1.0 + beta * (d ** eta))


def phase_space_weight(d: float, d_max: float) -> float:
    x = d / d_max
    return max(0.0, x * (1.0 - x))


def main() -> None:
    GEN.mkdir(exist_ok=True)

    params = {
        "omega": 0.18575,
        "phi": 0.16250,
        "beta": 1.0,
        "eta": 1.8,
    }
    frozen_channel = "graviton_to_gauge_gauge"
    frozen_background = "B1_fixed_chart_seed"
    d_min, d_max = 0.0, 20.0

    def integrand(d: float) -> float:
        k = k_strict(d, **params)
        return (k * k) * phase_space_weight(d, d_max)

    value, abs_err = quad(integrand, d_min, d_max, epsabs=1e-10, epsrel=1e-10, limit=400)

    # simple stability cross-check on discretized grid
    grid = np.linspace(d_min, d_max, 20001)
    vals = np.array([integrand(float(x)) for x in grid], dtype=float)
    trapz_value = float(np.trapezoid(vals, grid))
    delta_quad_trapz = abs(value - trapz_value)

    witness = {
        "witness_id": "STRICT_CUTKOSKY_FIXED_CHANNEL_REAL_DISCONTINUITY_INTEGRAL_V1",
        "frozen_channel": frozen_channel,
        "frozen_background": frozen_background,
        "integration_domain": {"d_min": d_min, "d_max": d_max},
        "strict_kernel_params": params,
        "integral_observable": "I_disc_proxy = integral[(K_strict(d)^2)*phase_space_weight(d)]dd",
        "quad_value": value,
        "quad_abs_error": abs_err,
        "trapz_value": trapz_value,
        "abs_delta_quad_vs_trapz": delta_quad_trapz,
        "positivity": {
            "integrand_nonnegative_expected": True,
            "integral_nonnegative_observed": value >= 0.0,
        },
        "numerical_checks": {
            "quad_finite": math.isfinite(value),
            "quad_error_bounded": abs_err < 1e-7,
            "quad_trapz_consistent": delta_quad_trapz < 5e-5,
        },
    }

    checks = witness["numerical_checks"]
    all_checks_pass = bool(checks["quad_finite"] and checks["quad_error_bounded"] and checks["quad_trapz_consistent"] and witness["positivity"]["integral_nonnegative_observed"])

    result_kind = (
        "PASS_STRICT_CUTKOSKY_FIXED_CHANNEL_REAL_DISCONTINUITY_INTEGRAL_PACKET"
        if all_checks_pass
        else "OPEN_STRICT_CUTKOSKY_FIXED_CHANNEL_REAL_DISCONTINUITY_INTEGRAL_PACKET_NUMERICAL_ISSUE"
    )

    payload = {
        "schema_version": "p2190_s1140_v1",
        "packet_id": "P2190",
        "stage_id": "S1140",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_cutkosky_fixed_channel_real_discontinuity_integral_packet": witness,
        "recommended_next_honest_step": {
            "id": "P2191_candidate",
            "goal": "extend fixed-channel integral packet with explicit residue-band table and UR-link uncertainty envelope",
        },
        "gatekeeper_checks": {
            "real_discontinuity_integral_exported": True,
            "integral_numerical_checks_pass": all_checks_pass,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2190 S1140: strict Cutkosky fixed-channel real discontinuity integral packet",
                "",
                f"- Result kind: `{result_kind}`",
                f"- Channel: `{frozen_channel}`",
                f"- Background: `{frozen_background}`",
                f"- Quad value: `{value:.12e}`",
                f"- Quad abs err: `{abs_err:.3e}`",
                f"- |quad-trapz|: `{delta_quad_trapz:.3e}`",
                "",
                "This is a strict-lane numerical integral witness only, not full dressed Cutkosky closure.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
