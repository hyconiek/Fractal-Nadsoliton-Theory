#!/usr/bin/env python3
"""Scratch probe: monotone pointwise inverse branch for normalized kernels.

The cutoff Puiseux branch failed global admissibility.  This probe asks a more
controlled question: ignoring the truncated local series, does the normalized
legacy profile itself admit a real monotone inverse branch that hits the strict
normalized values on d=0..11?

Yes, on the first decreasing legacy interval [0,2].  This is a useful existence
witness, but it is not a bridge theorem: the branch is defined by pointwise
output matching, not by a strict-side transport/RG law, and it compresses the
whole Z12 distance range into a short legacy coordinate interval.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
GLOBAL_STRESS = HERE / "bridge_phase_normalized_global_admissibility_report.json"
OUT_JSON = HERE / "bridge_phase_normalized_monotone_inverse_branch_report.json"
OUT_MD = HERE / "bridge_phase_normalized_monotone_inverse_branch_report.md"

LEGACY = {"omega": math.pi / 4.0, "phi": math.pi / 6.0, "beta_tors": 0.01}
STRICT = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 9.0 / 5.0}
BRACKET = (0.0, 2.0)
Z12_D = np.arange(0, 12, dtype=float)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def legacy_norm(x: float) -> float:
    return math.cos(LEGACY["omega"] * x + LEGACY["phi"]) / math.cos(LEGACY["phi"]) / (1.0 + LEGACY["beta_tors"] * x)


def legacy_derivative(x: float) -> float:
    omega = LEGACY["omega"]
    phi = LEGACY["phi"]
    beta = LEGACY["beta_tors"]
    theta = omega * x + phi
    numerator = -omega * math.sin(theta) * (1.0 + beta * x) - beta * math.cos(theta)
    denominator = math.cos(phi) * (1.0 + beta * x) ** 2
    return numerator / denominator


def strict_norm(d: float) -> float:
    return math.cos(STRICT["omega"] * d + STRICT["phi"]) / math.cos(STRICT["phi"]) / (1.0 + d ** STRICT["eta"])


def solve_first_branch(y: float, lo: float = BRACKET[0], hi: float = BRACKET[1]) -> float:
    f_lo = legacy_norm(lo) - y
    f_hi = legacy_norm(hi) - y
    if abs(f_lo) < 1e-15:
        return lo
    if f_lo * f_hi > 0:
        raise ValueError(f"target {y} not bracketed on [{lo}, {hi}]: {f_lo}, {f_hi}")
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        f_mid = legacy_norm(mid) - y
        if f_lo * f_mid <= 0:
            hi = mid
            f_hi = f_mid
        else:
            lo = mid
            f_lo = f_mid
    return 0.5 * (lo + hi)


def derivative_certificate() -> dict[str, Any]:
    grid = np.linspace(BRACKET[0], BRACKET[1], 2001)
    vals = np.array([legacy_derivative(float(x)) for x in grid])
    # Analytic coarse bound on [0,2]: theta in [pi/6, 2pi/3], sin(theta)>=1/2 and cos(theta)>=-1/2.
    omega = LEGACY["omega"]
    beta = LEGACY["beta_tors"]
    phi = LEGACY["phi"]
    numerator_upper_bound = -omega * 0.5 + beta * 0.5
    denominator_lower_bound = math.cos(phi)
    derivative_upper_bound = numerator_upper_bound / denominator_lower_bound
    return {
        "interval": list(BRACKET),
        "grid_max_derivative": float(np.max(vals)),
        "grid_min_derivative": float(np.min(vals)),
        "all_grid_derivatives_negative": bool(np.all(vals < 0.0)),
        "coarse_analytic_derivative_upper_bound": derivative_upper_bound,
        "coarse_bound_is_negative": derivative_upper_bound < 0.0,
        "bound_note": "On x in [0,2], theta=omega*x+phi lies in [pi/6,2pi/3], so sin(theta)>=1/2 and cos(theta)>=-1/2; this makes the derivative numerator strictly negative.",
    }


def branch_rows() -> list[dict[str, float | bool]]:
    rows = []
    for d in Z12_D:
        target = strict_norm(float(d))
        x = solve_first_branch(target)
        residual = legacy_norm(x) - target
        rows.append(
            {
                "d": float(d),
                "strict_norm": target,
                "x_first_branch": x,
                "legacy_norm_at_x": legacy_norm(x),
                "residual": residual,
                "abs_residual": abs(residual),
                "in_bracket": BRACKET[0] <= x <= BRACKET[1],
            }
        )
    return rows


def main() -> None:
    global_stress = load_json(GLOBAL_STRESS)
    cert = derivative_certificate()
    rows = branch_rows()
    x_values = [float(row["x_first_branch"]) for row in rows]
    dx = np.diff(np.array(x_values))
    max_residual = max(float(row["abs_residual"]) for row in rows)
    compression_ratio = x_values[-1] / float(Z12_D[-1])

    report = {
        "status": "OPEN_MONOTONE_INVERSE_BRANCH_EXISTENCE_WITNESS_NO_BRIDGE_THEOREM",
        "result_kind": "SCRATCH_PHASE_NORMALIZED_MONOTONE_INVERSE_BRANCH_PROBE__NOT_A_THEOREM",
        "source_reports": {"global_admissibility_stress": str(GLOBAL_STRESS.relative_to(HERE.parents[1]))},
        "branch_interval": list(BRACKET),
        "legacy_monotonicity_certificate": cert,
        "target_range": {
            "strict_min_on_z12": min(float(row["strict_norm"]) for row in rows),
            "strict_max_on_z12": max(float(row["strict_norm"]) for row in rows),
            "legacy_left_value": legacy_norm(BRACKET[0]),
            "legacy_right_value": legacy_norm(BRACKET[1]),
            "targets_bracketed": all(legacy_norm(BRACKET[0]) >= float(row["strict_norm"]) >= legacy_norm(BRACKET[1]) for row in rows),
        },
        "branch_rows": rows,
        "branch_metrics": {
            "max_abs_residual": max_residual,
            "all_roots_in_bracket": all(bool(row["in_bracket"]) for row in rows),
            "is_monotone_non_decreasing_on_z12": bool(np.all(dx >= -1e-12)),
            "x_min": min(x_values),
            "x_max": max(x_values),
            "z12_to_legacy_coordinate_compression_x11_over_11": compression_ratio,
            "cutoff_formal_inverse_global_admissible": global_stress["global_admissible_candidate"],
        },
        "honest_interpretation": [
            "The cutoff Puiseux branch failed, but a separate pointwise monotone inverse branch exists on the first legacy decreasing interval.",
            "This branch exactly matches normalized output values on d=0..11 to numerical precision and remains nonnegative/bounded.",
            "It is not a bridge theorem because it is defined by inverse matching, not derived from strict-side transport/RG, and it compresses Z12 distances into x<=~1.35.",
        ],
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No strict-side law derives this inverse branch.",
            "No proof that the inverse branch preserves physical roles or distance ontology is exported.",
            "No legacy physical-role transfer is licensed.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Treat the inverse branch as a candidate boundary condition: test whether a strict-side monotone transport/RG ODE can derive x'(d) and the observed compression without importing legacy roles.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch phase-normalized monotone inverse branch probe\n\n"
        "Status: pointwise monotone inverse branch exists; no bridge theorem.\n\n"
        f"- Legacy branch derivative negative on `{BRACKET}`: grid `{cert['all_grid_derivatives_negative']}`, coarse bound `{cert['coarse_analytic_derivative_upper_bound']:.6f}`.\n"
        f"- Z12 branch residual max `{max_residual:.3e}`, monotone `{bool(np.all(dx >= -1e-12))}`, roots in bracket `{all(bool(row['in_bracket']) for row in rows)}`.\n"
        f"- Compression: `x(11)={x_values[-1]:.12f}`, so `x(11)/11={compression_ratio:.12f}`.\n"
        "- Honest read: this rescues pointwise global admissibility but only as an inverse-defined branch, not a derived physical bridge.\n"
        "- No false pass: no kernel identity, no strict-side branch law, no physical-role transfer, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
