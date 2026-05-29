#!/usr/bin/env python3
"""Scratch probe: ODE form of the monotone inverse branch.

The monotone inverse-branch packet showed that normalized legacy outputs can be
matched pointwise to strict normalized outputs on d=0..11.  This packet takes a
more proof-oriented step: differentiating the output-matching condition

    L_norm(x(d)) = S_norm(d)

gives the transport equation

    x'(d) = S_norm'(d) / L_norm'(x(d)).

We integrate that ODE from x(0)=0 and compare it with the bisection-defined
inverse branch.  This is the closest thing in this scratch lane to an
"ontological bridge" candidate: it treats x(d) as an information-distance
compression flow.  It is still not a bridge theorem, because the ODE is obtained
from output matching rather than derived from strict-side nadsoliton dynamics.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
MONOTONE_BRANCH = HERE / "bridge_phase_normalized_monotone_inverse_branch_report.json"
OUT_JSON = HERE / "bridge_phase_normalized_inverse_branch_ode_report.json"
OUT_MD = HERE / "bridge_phase_normalized_inverse_branch_ode_report.md"

LEGACY = {"omega": math.pi / 4.0, "phi": math.pi / 6.0, "beta_tors": 0.01}
STRICT = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 9.0 / 5.0}
BRACKET = (0.0, 2.0)
D_MAX = 11.0
H = 1e-3
SAMPLE_D = np.array([0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 11.0], dtype=float)


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
    return (carrier_prime * denom - carrier * denom_prime) / (denom * denom)


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


def ode_rhs(d: float, x: float) -> float:
    return strict_derivative(d) / legacy_derivative(x)


def rk4_integrate() -> dict[str, Any]:
    n_steps = int(round(D_MAX / H))
    ds = np.linspace(0.0, D_MAX, n_steps + 1)
    xs = np.zeros(n_steps + 1, dtype=float)
    for idx in range(n_steps):
        d = float(ds[idx])
        x = float(xs[idx])
        k1 = ode_rhs(d, x)
        k2 = ode_rhs(d + 0.5 * H, x + 0.5 * H * k1)
        k3 = ode_rhs(d + 0.5 * H, x + 0.5 * H * k2)
        k4 = ode_rhs(d + H, x + H * k3)
        xs[idx + 1] = x + (H / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    branch = np.array([branch_x(float(d)) for d in ds])
    residual = xs - branch
    output_residual = np.array([legacy_norm(float(x)) - strict_norm(float(d)) for d, x in zip(ds, xs)])
    dx = np.diff(xs)
    return {
        "step": H,
        "n_steps": n_steps,
        "max_abs_x_residual_vs_bisection_branch": float(np.max(np.abs(residual))),
        "max_abs_output_residual": float(np.max(np.abs(output_residual))),
        "is_ode_solution_nonnegative": bool(np.all(xs >= -1e-12)),
        "is_ode_solution_monotone_non_decreasing": bool(np.all(dx >= -1e-10)),
        "x_end": float(xs[-1]),
        "branch_end": float(branch[-1]),
        "compression_x11_over_11": float(xs[-1] / D_MAX),
    }


def sample_rows() -> list[dict[str, float | bool]]:
    rows = []
    for d in SAMPLE_D:
        x = branch_x(float(d))
        xp = ode_rhs(float(d), x)
        lhs = legacy_derivative(x) * xp
        rhs = strict_derivative(float(d))
        rows.append(
            {
                "d": float(d),
                "x_branch": x,
                "x_prime_from_ode": xp,
                "legacy_derivative": legacy_derivative(x),
                "strict_derivative": rhs,
                "differentiated_identity_residual": lhs - rhs,
                "positive_flow": bool(xp >= 0.0),
                "compression_x_over_d": None if d == 0.0 else x / float(d),
            }
        )
    return rows


def main() -> None:
    previous = load_json(MONOTONE_BRANCH)
    integration = rk4_integrate()
    rows = sample_rows()
    positive_samples = all(bool(row["positive_flow"]) for row in rows)

    report = {
        "status": "OPEN_INVERSE_BRANCH_ODE_CANDIDATE_NO_ONTOLOGICAL_BRIDGE_THEOREM",
        "result_kind": "SCRATCH_PHASE_NORMALIZED_INVERSE_BRANCH_ODE_PROBE__NOT_A_THEOREM",
        "source_reports": {"monotone_inverse_branch": str(MONOTONE_BRANCH.relative_to(HERE.parents[1]))},
        "transport_equation": {
            "output_matching_identity": "L_norm(x(d)) = S_strict_norm(d)",
            "derived_ode": "x'(d)=S_strict_norm'(d)/L_norm'(x(d))",
            "initial_condition": "x(0)=0",
            "derivation_status": "DIFFERENTIATED_OUTPUT_MATCHING_NOT_STRICT_SIDE_DYNAMICS",
        },
        "integration_check": integration,
        "sample_rows": rows,
        "candidate_ontological_reading": {
            "name": "information_distance_compression_flow_candidate",
            "supported_by_this_probe": bool(
                integration["is_ode_solution_nonnegative"]
                and integration["is_ode_solution_monotone_non_decreasing"]
                and integration["max_abs_output_residual"] < 1e-6
                and positive_samples
            ),
            "content": "If a bridge exists, x(d) may be interpreted as a compressed information-distance coordinate generated by a monotone flow.",
            "not_yet_proven_because": "The flow is derived from matching two normalized outputs, not from nadsoliton ontology or strict-side transport/RG equations.",
        },
        "upstream_replay": {
            "previous_branch_max_abs_residual": previous["branch_metrics"]["max_abs_residual"],
            "previous_branch_monotone": previous["branch_metrics"]["is_monotone_non_decreasing_on_z12"],
            "previous_cutoff_formal_inverse_global_admissible": previous["branch_metrics"]["cutoff_formal_inverse_global_admissible"],
        },
        "honest_interpretation": [
            "The monotone inverse branch can be encoded as a first-order compression-flow ODE and numerically integrated back to the same branch.",
            "This is a stronger bridge candidate than an unconstrained local Puiseux fit because it gives a global flow with an initial condition.",
            "It is still not an ontological bridge theorem: the ODE is tautologically induced by output matching unless strict-side dynamics derives it independently.",
        ],
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No strict-side derivation of x'(d) is exported.",
            "No proof that x(d) is the nadsoliton information metric is exported.",
            "No legacy physical-role transfer is licensed.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Try to derive the same compression-flow ODE from strict-side information geometry; if that fails, classify the monotone branch as a reparameterization artifact rather than an ontological bridge.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch phase-normalized inverse-branch ODE probe\n\n"
        "Status: compression-flow ODE candidate; no ontological bridge theorem.\n\n"
        f"- RK4 integration max output residual `{integration['max_abs_output_residual']:.3e}` and max x residual vs bisection `{integration['max_abs_x_residual_vs_bisection_branch']:.3e}`.\n"
        f"- ODE solution nonnegative `{integration['is_ode_solution_nonnegative']}`, monotone `{integration['is_ode_solution_monotone_non_decreasing']}`, x(11)/11 `{integration['compression_x11_over_11']:.12f}`.\n"
        f"- Candidate ontology support flag `{report['candidate_ontological_reading']['supported_by_this_probe']}`: compressed information-distance flow, not strict-side-derived.\n"
        "- Honest read: stronger than local Puiseux fitting, but still output-matching tautology until derived from strict information geometry.\n"
        "- No false pass: no kernel identity, no strict-side x'(d) derivation, no information-metric theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
