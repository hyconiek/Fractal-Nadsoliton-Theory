#!/usr/bin/env python3
"""Scratch probe: curvature certificate for the nonlinear compression flow.

The affine phase-alignment packet blocked simple affine bridges.  This probe
makes that obstruction differential: for the monotone inverse branch defined by

    L_norm(x(d)) = S_strict_norm(d),

implicit differentiation gives x' and x''.  A genuinely affine bridge would have
x''=0.  The computed branch has large positive curvature near the origin and
negative curvature across the Z12 tail, so the bridge candidate is necessarily a
nonlinear compression flow.  This still does not prove an ontological bridge;
it only records a sharper constraint that any strict-side derivation must meet.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
PHASE_OBSTRUCTION = HERE / "bridge_phase_normalized_phase_alignment_obstruction_report.json"
INVERSE_ODE = HERE / "bridge_phase_normalized_inverse_branch_ode_report.json"
OUT_JSON = HERE / "bridge_phase_normalized_compression_curvature_report.json"
OUT_MD = HERE / "bridge_phase_normalized_compression_curvature_report.md"

LEGACY = {"omega": math.pi / 4.0, "phi": math.pi / 6.0, "beta_tors": 0.01}
STRICT = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 9.0 / 5.0}
BRACKET = (0.0, 2.0)
SAMPLE_D = np.array([1e-4, 1e-3, 1e-2, 0.1, 0.5, 1.0, 2.0, 4.0, 8.0, 11.0], dtype=float)
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


def branch_x(d: float) -> float:
    return solve_first_branch(strict_norm(d))


def derivative_row(d: float) -> dict[str, float | bool]:
    x = branch_x(d)
    lp = legacy_derivative(x)
    sp = strict_derivative(d)
    xp = sp / lp
    lpp = legacy_second_derivative(x)
    spp = strict_second_derivative(d)
    xpp = (spp - lpp * xp**2) / lp
    residual = lpp * xp**2 + lp * xpp - spp
    return {
        "d": d,
        "x": x,
        "x_prime": xp,
        "x_second": xpp,
        "legacy_prime": lp,
        "legacy_second": lpp,
        "strict_prime": sp,
        "strict_second": spp,
        "second_derivative_identity_residual": residual,
        "curvature_nonzero": abs(xpp) > 1e-6,
        "curvature_sign": 1.0 if xpp > 0 else -1.0 if xpp < 0 else 0.0,
    }


def z12_discrete_curvature() -> dict[str, Any]:
    xs = np.array([branch_x(float(d)) for d in Z12_D])
    second_diff = xs[2:] - 2.0 * xs[1:-1] + xs[:-2]
    return {
        "x_values": [float(v) for v in xs],
        "second_differences_d1_to_d10": [float(v) for v in second_diff],
        "all_tail_second_differences_negative_from_d1": bool(np.all(second_diff[0:] < 0.0)),
        "max_abs_second_difference": float(np.max(np.abs(second_diff))),
    }


def main() -> None:
    phase_obstruction = load_json(PHASE_OBSTRUCTION)
    inverse_ode = load_json(INVERSE_ODE)
    rows = [derivative_row(float(d)) for d in SAMPLE_D]
    max_identity_residual = max(abs(float(row["second_derivative_identity_residual"])) for row in rows)
    nonzero_count = sum(bool(row["curvature_nonzero"]) for row in rows)
    has_positive = any(float(row["x_second"]) > 0.0 for row in rows)
    has_negative = any(float(row["x_second"]) < 0.0 for row in rows)
    discrete = z12_discrete_curvature()

    report = {
        "status": "OPEN_NONLINEAR_COMPRESSION_CURVATURE_CERTIFICATE_NO_BRIDGE_THEOREM",
        "result_kind": "SCRATCH_PHASE_NORMALIZED_COMPRESSION_CURVATURE_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "phase_alignment_obstruction": str(PHASE_OBSTRUCTION.relative_to(HERE.parents[1])),
            "inverse_branch_ode": str(INVERSE_ODE.relative_to(HERE.parents[1])),
        },
        "implicit_derivative_formulas": {
            "first": "x'=S'(d)/L'(x)",
            "second": "x''=(S''(d)-L''(x)*(x')**2)/L'(x)",
            "affine_bridge_test": "affine bridge would require x''=0",
        },
        "sample_rows": rows,
        "discrete_z12_curvature": discrete,
        "curvature_certificate": {
            "max_second_derivative_identity_residual": max_identity_residual,
            "nonzero_curvature_sample_count": nonzero_count,
            "sample_count": len(rows),
            "has_positive_curvature_near_origin": has_positive,
            "has_negative_curvature_on_tail": has_negative,
            "affine_bridge_ruled_out_by_curvature": nonzero_count == len(rows) and has_positive and has_negative,
        },
        "upstream_replay": {
            "affine_routes_blocked": phase_obstruction["restricted_obstruction"]["all_simple_affine_routes_blocked"],
            "ode_candidate_supported": inverse_ode["candidate_ontological_reading"]["supported_by_this_probe"],
            "ode_derivation_status": inverse_ode["transport_equation"]["derivation_status"],
        },
        "honest_interpretation": [
            "The inverse-branch compression flow has nonzero curvature; it cannot be reduced to any affine phase/distance map.",
            "The curvature is positive near the origin and negative on the Z12 tail, so the flow bends from initial acceleration into horizon approach.",
            "This is a constraint on a future ontological bridge, not a strict-side derivation of that bridge.",
        ],
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No strict-side derivation of the curvature law is exported.",
            "No proof that the curved compression flow is the nadsoliton information metric is exported.",
            "No legacy physical-role transfer is licensed.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "A real ontological bridge must derive both the first-order compression ODE and this nonzero curvature profile from strict-side information geometry, not from output matching.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch phase-normalized compression curvature probe\n\n"
        "Status: nonlinear curvature certificate; no bridge theorem.\n\n"
        f"- Curvature nonzero samples `{nonzero_count}/{len(rows)}`; second-derivative identity residual max `{max_identity_residual:.3e}`.\n"
        f"- Positive near-origin curvature `{has_positive}` and negative tail curvature `{has_negative}`; affine bridge ruled out `{report['curvature_certificate']['affine_bridge_ruled_out_by_curvature']}`.\n"
        f"- Z12 max abs second difference `{discrete['max_abs_second_difference']:.3e}`; tail second differences negative `{discrete['all_tail_second_differences_negative_from_d1']}`.\n"
        "- Honest read: any ontological bridge must derive a nonlinear curved compression flow, not only a monotone branch.\n"
        "- No false pass: no kernel identity, no strict curvature derivation, no information-metric theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
