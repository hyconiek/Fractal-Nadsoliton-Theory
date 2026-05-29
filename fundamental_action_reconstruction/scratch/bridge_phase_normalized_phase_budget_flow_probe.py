#!/usr/bin/env python3
"""Scratch probe: normalized phase-budget compression flow.

The horizon and curvature probes made the ontological candidate sharper: strict
Z12 distance is compressed into a finite legacy-side horizon x_h=4/3 by a
nonlinear inverse-branch flow.  This packet rewrites the branch in normalized
phase-budget coordinates

    v = d / d_h,      u = x(d) / x_h,

where d_h is the strict carrier-zero distance and x_h is the legacy carrier-zero
horizon.  It checks endpoint identities, the normalized compression density
du/dv, the inflection point of the flow, and how badly the naive affine phase
budget u=v fails.  It is still a bridge-prep object only: the flow is induced by
output matching unless strict-side information geometry derives it.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
HORIZON = HERE / "bridge_phase_normalized_compression_horizon_report.json"
CURVATURE = HERE / "bridge_phase_normalized_compression_curvature_report.json"
OUT_JSON = HERE / "bridge_phase_normalized_phase_budget_flow_report.json"
OUT_MD = HERE / "bridge_phase_normalized_phase_budget_flow_report.md"

LEGACY = {"omega": math.pi / 4.0, "phi": math.pi / 6.0, "beta_tors": 0.01}
STRICT = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 9.0 / 5.0}
BRACKET = (0.0, 2.0)
GRID_N = 1201
SAMPLE_V = np.array([0.0, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0], dtype=float)


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


def legacy_second_derivative(x: float) -> float:
    omega = LEGACY["omega"]
    phi = LEGACY["phi"]
    beta = LEGACY["beta_tors"]
    theta = omega * x + phi
    cos_phi = math.cos(phi)
    n = math.cos(theta) / cos_phi
    np1 = -omega * math.sin(theta) / cos_phi
    np2 = -(omega**2) * math.cos(theta) / cos_phi
    den = 1.0 + beta * x
    return np2 / den - 2.0 * np1 * beta / den**2 + 2.0 * n * beta**2 / den**3


def strict_norm(d: float) -> float:
    return math.cos(STRICT["omega"] * d + STRICT["phi"]) / math.cos(STRICT["phi"]) / (1.0 + d ** STRICT["eta"])


def strict_derivative(d: float) -> float:
    omega = STRICT["omega"]
    phi = STRICT["phi"]
    eta = STRICT["eta"]
    theta = omega * d + phi
    carrier = math.cos(theta) / math.cos(phi)
    carrier_prime = -omega * math.sin(theta) / math.cos(phi)
    denom = 1.0 + d**eta
    denom_prime = 0.0 if d == 0.0 else eta * d ** (eta - 1.0)
    return (carrier_prime * denom - carrier * denom_prime) / denom**2


def strict_second_derivative(d: float) -> float:
    omega = STRICT["omega"]
    phi = STRICT["phi"]
    eta = STRICT["eta"]
    theta = omega * d + phi
    carrier = math.cos(theta) / math.cos(phi)
    carrier_prime = -omega * math.sin(theta) / math.cos(phi)
    carrier_second = -(omega**2) * math.cos(theta) / math.cos(phi)
    denom = 1.0 + d**eta
    denom_prime = eta * d ** (eta - 1.0)
    denom_second = eta * (eta - 1.0) * d ** (eta - 2.0)
    return (carrier_second * denom - carrier * denom_second) / denom**2 - 2.0 * (carrier_prime * denom - carrier * denom_prime) * denom_prime / denom**3


def solve_first_branch(y: float, lo: float = BRACKET[0], hi: float = BRACKET[1]) -> float:
    f_lo = legacy_norm(lo) - y
    if abs(f_lo) < 1e-15:
        return lo
    for _ in range(90):
        mid = 0.5 * (lo + hi)
        f_mid = legacy_norm(mid) - y
        if f_lo * f_mid <= 0:
            hi = mid
        else:
            lo = mid
            f_lo = f_mid
    return 0.5 * (lo + hi)


def branch_x(d: float) -> float:
    return solve_first_branch(strict_norm(d))


def x_prime(d: float, x: float) -> float:
    return strict_derivative(d) / legacy_derivative(x)


def x_second(d: float, x: float) -> float:
    xp = x_prime(d, x)
    return (strict_second_derivative(d) - legacy_second_derivative(x) * xp**2) / legacy_derivative(x)


def inflection_point(d_h: float) -> dict[str, float]:
    lo, hi = 0.1, 0.5
    f_lo = x_second(lo, branch_x(lo))
    f_hi = x_second(hi, branch_x(hi))
    if f_lo * f_hi > 0:
        raise ValueError("inflection not bracketed")
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        f_mid = x_second(mid, branch_x(mid))
        if f_lo * f_mid <= 0:
            hi = mid
            f_hi = f_mid
        else:
            lo = mid
            f_lo = f_mid
    d_star = 0.5 * (lo + hi)
    x_star = branch_x(d_star)
    return {"d_inflection": d_star, "v_inflection": d_star / d_h, "x_inflection": x_star}


def sample_rows(d_h: float, x_h: float) -> list[dict[str, float]]:
    rows = []
    for v in SAMPLE_V:
        d = float(v * d_h)
        x = branch_x(d)
        xp = x_prime(d, x)
        u = x / x_h
        density = xp * d_h / x_h
        rows.append(
            {
                "v": float(v),
                "d": d,
                "u": u,
                "u_minus_v": u - float(v),
                "du_dv": density,
                "affine_abs_error": abs(u - float(v)),
            }
        )
    return rows


def main() -> None:
    horizon_report = load_json(HORIZON)
    curvature_report = load_json(CURVATURE)
    x_h = float(horizon_report["horizon"]["x_horizon"])
    d_h = float(horizon_report["horizon"]["strict_zero_d"])
    grid_v = np.linspace(0.0, 1.0, GRID_N)
    grid_d = grid_v * d_h
    grid_x = np.array([branch_x(float(d)) for d in grid_d])
    grid_u = grid_x / x_h
    affine_error = grid_u - grid_v
    densities = np.array([x_prime(float(d), float(x)) * d_h / x_h for d, x in zip(grid_d, grid_x)])
    integral_density = float(np.trapezoid(densities, grid_v))
    inflection = inflection_point(d_h)
    rows = sample_rows(d_h, x_h)

    report = {
        "status": "OPEN_PHASE_BUDGET_COMPRESSION_FLOW_TRACE_NO_ONTOLOGICAL_BRIDGE_THEOREM",
        "result_kind": "SCRATCH_PHASE_NORMALIZED_PHASE_BUDGET_FLOW_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "compression_horizon": str(HORIZON.relative_to(HERE.parents[1])),
            "compression_curvature": str(CURVATURE.relative_to(HERE.parents[1])),
        },
        "normalized_coordinates": {
            "u": "x(d)/x_h",
            "v": "d/d_h",
            "x_h": x_h,
            "d_h": d_h,
            "endpoint_residual_u0": float(grid_u[0]),
            "endpoint_residual_u1_minus_1": float(grid_u[-1] - 1.0),
        },
        "flow_density": {
            "definition": "du/dv=(dx/dd)*(d_h/x_h)",
            "density_min": float(np.min(densities)),
            "density_max": float(np.max(densities)),
            "density_at_v0": float(densities[0]),
            "density_at_v1": float(densities[-1]),
            "trapezoid_integral_density": integral_density,
            "integral_residual_vs_1": integral_density - 1.0,
        },
        "affine_budget_failure": {
            "max_abs_u_minus_v": float(np.max(np.abs(affine_error))),
            "mean_abs_u_minus_v": float(np.mean(np.abs(affine_error))),
            "v_at_max_abs_error": float(grid_v[int(np.argmax(np.abs(affine_error)))]),
            "affine_phase_budget_acceptable": bool(np.max(np.abs(affine_error)) < 1e-2),
        },
        "inflection": inflection,
        "sample_rows": rows,
        "upstream_replay": {
            "horizon_candidate_supported": horizon_report["candidate_ontological_reading"]["supported_by_this_probe"],
            "curvature_affine_ruled_out": curvature_report["curvature_certificate"]["affine_bridge_ruled_out_by_curvature"],
            "curvature_positive_near_origin": curvature_report["curvature_certificate"]["has_positive_curvature_near_origin"],
            "curvature_negative_on_tail": curvature_report["curvature_certificate"]["has_negative_curvature_on_tail"],
        },
        "candidate_ontological_reading": {
            "name": "normalized_phase_budget_compression_flow_candidate",
            "supported_by_this_probe": bool(abs(grid_u[0]) < 1e-14 and abs(grid_u[-1] - 1.0) < 1e-12 and np.max(np.abs(affine_error)) > 1e-1),
            "content": "The bridge candidate can be written as a nonlinear normalized flow from strict phase budget v to finite information-horizon budget u.",
            "not_yet_proven_because": "The normalized flow is still induced by output matching, not derived from strict-side information geometry.",
        },
        "honest_interpretation": [
            "The normalized compression flow exactly maps endpoints u(0)=0 and u(1)=1 and has density integrating to one.",
            "The flow is strongly non-affine, with a large u-v deviation and a resolved inflection point.",
            "This is a useful ontological coordinate candidate, but a strict-side derivation of the flow remains missing.",
        ],
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No strict-side derivation of the normalized phase-budget flow is exported.",
            "No proof that u is the nadsoliton information-horizon coordinate is exported.",
            "No legacy physical-role transfer is licensed.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Try to derive the normalized density du/dv from strict-side information geometry; otherwise classify the phase-budget flow as an output-matching coordinate artifact.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch phase-normalized phase-budget flow probe\n\n"
        "Status: normalized nonlinear compression-flow trace; no bridge theorem.\n\n"
        f"- Endpoints: `u(0)={grid_u[0]:.3e}`, `u(1)-1={grid_u[-1]-1.0:.3e}`; density integral residual `{integral_density - 1.0:.3e}`.\n"
        f"- Affine budget failure max `|u-v|={np.max(np.abs(affine_error)):.3e}` at `v={grid_v[int(np.argmax(np.abs(affine_error)))]:.6f}`.\n"
        f"- Inflection at `d={inflection['d_inflection']:.12f}`, `v={inflection['v_inflection']:.12f}`, `x={inflection['x_inflection']:.12f}`.\n"
        "- Honest read: a finite normalized information-budget flow exists as a candidate, but remains output-matching until strict-side-derived.\n"
        "- No false pass: no kernel identity, no strict flow derivation, no information-coordinate theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
