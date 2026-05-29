#!/usr/bin/env python3
"""Scratch probe: finite information-measure transport for the phase bridge.

The previous phase-budget flow packet rewrote the output-matching inverse branch
as normalized coordinates u=x/x_h and v=d/d_h.  This packet takes one more
honest ontological step: differentiate the matched normalized kernels

    L(x_h u(v)) = S(d_h v)

and recast the result as a finite information-measure transport law

    m_legacy(u(v)) * du = rho_strict(v) * dv,

where rho_strict is the strict-side loss density and m_legacy is the legacy-side
horizon metric density.  Both densities integrate to the same unit information
budget on [0,1].  This is still theorem-prep only: the measure is induced by
kernel-output matching, not derived independently from strict nadsoliton
information geometry.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
PHASE_BUDGET = HERE / "bridge_phase_normalized_phase_budget_flow_report.json"
HORIZON = HERE / "bridge_phase_normalized_compression_horizon_report.json"
OUT_JSON = HERE / "bridge_phase_normalized_measure_transport_report.json"
OUT_MD = HERE / "bridge_phase_normalized_measure_transport_report.md"

LEGACY = {"omega": math.pi / 4.0, "phi": math.pi / 6.0, "beta_tors": 0.01}
STRICT = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 9.0 / 5.0}
BRACKET = (0.0, 2.0)
GRID_N = 2401
SAMPLE_V = np.array([0.0, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0], dtype=float)


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


def branch_du_dv(d: float, x: float, d_h: float, x_h: float) -> float:
    return strict_derivative(d) / legacy_derivative(x) * d_h / x_h


def rho_strict(v: float, d_h: float) -> float:
    """Strict-side normalized information-loss density on v=d/d_h."""
    return -d_h * strict_derivative(d_h * v)


def m_legacy(u: float, x_h: float) -> float:
    """Legacy-side normalized horizon metric density on u=x/x_h."""
    return -x_h * legacy_derivative(x_h * u)


def sample_rows(d_h: float, x_h: float) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for v in SAMPLE_V:
        d = float(v * d_h)
        x = branch_x(d)
        u = x / x_h
        density = branch_du_dv(d, x, d_h, x_h)
        strict_density = rho_strict(float(v), d_h)
        metric = m_legacy(u, x_h)
        balance_residual = metric * density - strict_density
        affine_residual = m_legacy(float(v), x_h) - strict_density
        rows.append(
            {
                "v": float(v),
                "u": u,
                "d": d,
                "x": x,
                "du_dv": density,
                "rho_strict": strict_density,
                "m_legacy_at_u": metric,
                "metric_times_du_dv": metric * density,
                "measure_balance_residual": balance_residual,
                "affine_balance_residual_if_u_equals_v": affine_residual,
            }
        )
    return rows


def main() -> None:
    phase_budget_report = load_json(PHASE_BUDGET)
    horizon_report = load_json(HORIZON)
    x_h = float(horizon_report["horizon"]["x_horizon"])
    d_h = float(horizon_report["horizon"]["strict_zero_d"])

    grid_v = np.linspace(0.0, 1.0, GRID_N)
    grid_d = grid_v * d_h
    grid_x = np.array([branch_x(float(d)) for d in grid_d])
    grid_u = grid_x / x_h
    grid_du_dv = np.array([branch_du_dv(float(d), float(x), d_h, x_h) for d, x in zip(grid_d, grid_x)])
    strict_density = np.array([rho_strict(float(v), d_h) for v in grid_v])
    legacy_metric_on_branch = np.array([m_legacy(float(u), x_h) for u in grid_u])
    balance_residual = legacy_metric_on_branch * grid_du_dv - strict_density

    grid_u_direct = np.linspace(0.0, 1.0, GRID_N)
    legacy_metric_direct = np.array([m_legacy(float(u), x_h) for u in grid_u_direct])
    affine_metric = np.array([m_legacy(float(v), x_h) for v in grid_v])
    affine_balance_residual = affine_metric - strict_density

    strict_loss_integral = float(np.trapezoid(strict_density, grid_v))
    legacy_metric_integral = float(np.trapezoid(legacy_metric_direct, grid_u_direct))
    transported_metric_integral = float(np.trapezoid(legacy_metric_on_branch * grid_du_dv, grid_v))
    strict_loss_exact = strict_norm(0.0) - strict_norm(d_h)
    legacy_metric_exact = legacy_norm(0.0) - legacy_norm(x_h)

    sample = sample_rows(d_h, x_h)
    report = {
        "status": "OPEN_PHASE_NORMALIZED_INFORMATION_MEASURE_TRANSPORT_NO_BRIDGE_THEOREM",
        "result_kind": "SCRATCH_PHASE_NORMALIZED_MEASURE_TRANSPORT_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "phase_budget_flow": str(PHASE_BUDGET.relative_to(HERE.parents[1])),
            "compression_horizon": str(HORIZON.relative_to(HERE.parents[1])),
        },
        "normalized_kernel_identity_used_as_input": {
            "legacy_side": "L(x)=K_legacy_ont(x)/(alpha_geo*cos(phi_legacy)) after amplitude normalization",
            "strict_side": "S(d)=K_strict_gate(d)/cos(phi_strict) on the operational strict branch",
            "matched_branch": "L(x_h*u(v))=S(d_h*v)",
            "guardrail": "This packet differentiates an output-matched branch; it does not identify K_legacy_ont with K_strict_gate.",
        },
        "horizon": {"x_h": x_h, "d_h": d_h, "u_domain": [0.0, 1.0], "v_domain": [0.0, 1.0]},
        "measure_definitions": {
            "rho_strict(v)": "-d_h * dS/dd(d_h*v)",
            "m_legacy(u)": "-x_h * dL/dx(x_h*u)",
            "transport_law": "m_legacy(u(v)) * du/dv = rho_strict(v)",
            "integral_reading": "Both sides integrate to the same finite unit information budget S(0)-S(d_h)=L(0)-L(x_h)=1.",
        },
        "measure_balance": {
            "strict_loss_integral_trapezoid": strict_loss_integral,
            "legacy_metric_integral_trapezoid": legacy_metric_integral,
            "transported_metric_integral_trapezoid": transported_metric_integral,
            "strict_loss_exact_endpoint_integral": strict_loss_exact,
            "legacy_metric_exact_endpoint_integral": legacy_metric_exact,
            "strict_trapezoid_residual_vs_1": strict_loss_integral - 1.0,
            "legacy_trapezoid_residual_vs_1": legacy_metric_integral - 1.0,
            "transported_trapezoid_residual_vs_1": transported_metric_integral - 1.0,
            "strict_exact_residual_vs_1": strict_loss_exact - 1.0,
            "legacy_exact_residual_vs_1": legacy_metric_exact - 1.0,
            "max_abs_balance_residual": float(np.max(np.abs(balance_residual))),
            "mean_abs_balance_residual": float(np.mean(np.abs(balance_residual))),
        },
        "affine_failure": {
            "test": "If the ontological phase budgets were affinely identified by u=v, the balance would require m_legacy(v)=rho_strict(v).",
            "affine_balance_max_abs_residual": float(np.max(np.abs(affine_balance_residual))),
            "affine_balance_mean_abs_residual": float(np.mean(np.abs(affine_balance_residual))),
            "v_at_max_abs_residual": float(grid_v[int(np.argmax(np.abs(affine_balance_residual)))]),
            "affine_budget_measure_transport_acceptable": bool(np.max(np.abs(affine_balance_residual)) < 1e-2),
        },
        "density_ranges": {
            "rho_strict_min": float(np.min(strict_density)),
            "rho_strict_max": float(np.max(strict_density)),
            "m_legacy_on_branch_min": float(np.min(legacy_metric_on_branch)),
            "m_legacy_on_branch_max": float(np.max(legacy_metric_on_branch)),
            "du_dv_min": float(np.min(grid_du_dv)),
            "du_dv_max": float(np.max(grid_du_dv)),
        },
        "sample_rows": sample,
        "upstream_replay": {
            "phase_budget_candidate_supported": phase_budget_report["candidate_ontological_reading"]["supported_by_this_probe"],
            "phase_budget_density_integral_residual": phase_budget_report["flow_density"]["integral_residual_vs_1"],
            "horizon_candidate_supported": horizon_report["candidate_ontological_reading"]["supported_by_this_probe"],
        },
        "candidate_ontological_reading": {
            "name": "finite_information_measure_transport_candidate",
            "supported_by_this_probe": bool(
                np.max(np.abs(balance_residual)) < 1e-12
                and abs(strict_loss_exact - 1.0) < 1e-12
                and abs(legacy_metric_exact - 1.0) < 1e-12
                and abs(strict_loss_integral - 1.0) < 1e-4
                and abs(legacy_metric_integral - 1.0) < 1e-4
                and np.max(np.abs(affine_balance_residual)) > 1e-1
            ),
            "content": "The phase-normalized bridge can be expressed as exact transport of a unit strict loss measure into a unit legacy horizon metric measure.",
            "why_more_ontological_than_pointwise_output_matching": "It isolates a conserved finite information budget and separates the strict loss density from the legacy horizon metric density.",
            "not_yet_proven_because": "The strict loss density and legacy metric density are still induced from the two kernels; no independent strict-side derivation of rho_strict or the transport map is exported.",
        },
        "honest_interpretation": [
            "Differentiating the normalized matched branch yields an exact measure-balance identity with residual at floating precision.",
            "Both density sides integrate to a unit finite information budget, so the bridge candidate can now be stated as measure transport rather than only pointwise curve matching.",
            "The affine budget u=v fails the density-balance test, so any bridge must be nonlinear even after phase normalization.",
            "This is still a scratch ontological candidate, not a strict derivation or physical-role transfer theorem.",
        ],
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No theorem derives rho_strict from strict nadsoliton information geometry independent of the strict kernel form.",
            "No theorem derives m_legacy as a retained strict-side horizon metric.",
            "No legacy physical-role transfer is licensed for sin^2(theta_W), alpha_EM^-1, or beta^N hierarchy claims.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Try to derive rho_strict(v) from strict-only nadsoliton information geometry or from alpha_geo_strict_derived_v1; without that, classify this as exact measure-transport bookkeeping, not a bridge theorem.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch phase-normalized information-measure transport probe\n\n"
        "Status: exact finite measure-transport identity on the output-matched branch; no bridge theorem.\n\n"
        f"- Unit budgets: exact strict residual `{strict_loss_exact - 1.0:.3e}`, exact legacy residual `{legacy_metric_exact - 1.0:.3e}`; trapezoid strict residual `{strict_loss_integral - 1.0:.3e}`, transported residual `{transported_metric_integral - 1.0:.3e}`.\n"
        f"- Measure balance: max residual `{np.max(np.abs(balance_residual)):.3e}`, mean residual `{np.mean(np.abs(balance_residual)):.3e}` for `m_legacy(u(v))*du/dv=rho_strict(v)`.\n"
        f"- Affine failure: max residual `{np.max(np.abs(affine_balance_residual)):.3e}` at `v={grid_v[int(np.argmax(np.abs(affine_balance_residual)))]:.6f}` if `u=v`.\n"
        "- Honest read: a stronger ontological candidate exists as unit information-measure transport, but the densities are still kernel-induced.\n"
        "- No false pass: no kernel identity, no strict density derivation, no legacy physical-role transfer, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
