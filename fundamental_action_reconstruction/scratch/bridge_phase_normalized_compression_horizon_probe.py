#!/usr/bin/env python3
"""Scratch probe: compression horizon of the monotone inverse branch.

The inverse-branch ODE suggests an ontological candidate: strict Z12 distance d
is compressed into a legacy-side information coordinate x(d).  This probe asks
what global structure that coordinate has on the controlled first branch.

Because the first zero of the normalized legacy carrier occurs at

    omega_legacy*x + phi_legacy = pi/2,

there is a closed-form compression horizon x_h = 4/3.  The monotone inverse
branch crosses this horizon exactly when the strict normalized kernel changes
sign.  This gives a sharper, computable bridge-prep object: a finite information
horizon with a local linearization around x_h.  It is still not an ontological
bridge theorem unless strict-side information geometry derives this horizon and
flow independently.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
INVERSE_ODE = HERE / "bridge_phase_normalized_inverse_branch_ode_report.json"
OUT_JSON = HERE / "bridge_phase_normalized_compression_horizon_report.json"
OUT_MD = HERE / "bridge_phase_normalized_compression_horizon_report.md"

LEGACY = {"omega": math.pi / 4.0, "phi": math.pi / 6.0, "beta_tors": 0.01}
STRICT = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 9.0 / 5.0}
BRACKET = (0.0, 2.0)
Z12_D = np.arange(0, 12, dtype=float)
TAIL_D = np.array([6.0, 7.0, 8.0, 9.0, 10.0, 11.0], dtype=float)


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


def branch_x(d: float) -> float:
    return solve_first_branch(strict_norm(d))


def strict_zero() -> float:
    lo, hi = 7.0, 8.0
    f_lo, f_hi = strict_norm(lo), strict_norm(hi)
    if f_lo * f_hi > 0:
        raise ValueError("strict zero not bracketed in [7,8]")
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        f_mid = strict_norm(mid)
        if f_lo * f_mid <= 0:
            hi = mid
            f_hi = f_mid
        else:
            lo = mid
            f_lo = f_mid
    return 0.5 * (lo + hi)


def horizon_data() -> dict[str, Any]:
    x_h = (math.pi / 2.0 - LEGACY["phi"]) / LEGACY["omega"]
    derivative = legacy_derivative(x_h)
    d_zero = strict_zero()
    return {
        "x_horizon_formula": "(pi/2 - phi_legacy)/omega_legacy",
        "x_horizon": x_h,
        "x_horizon_expected_exact": 4.0 / 3.0,
        "horizon_formula_residual_vs_4_over_3": x_h - 4.0 / 3.0,
        "legacy_norm_at_horizon": legacy_norm(x_h),
        "legacy_derivative_at_horizon": derivative,
        "strict_zero_d": d_zero,
        "branch_x_at_strict_zero": branch_x(d_zero),
        "branch_horizon_residual_at_strict_zero": branch_x(d_zero) - x_h,
    }


def z12_rows(horizon: dict[str, Any]) -> list[dict[str, float | bool]]:
    x_h = float(horizon["x_horizon"])
    lprime = float(horizon["legacy_derivative_at_horizon"])
    rows = []
    for d in Z12_D:
        s = strict_norm(float(d))
        x = branch_x(float(d))
        linearized = x_h + s / lprime
        rows.append(
            {
                "d": float(d),
                "strict_norm": s,
                "x_branch": x,
                "x_minus_horizon": x - x_h,
                "same_sign_as_strict_norm_after_lprime": bool((x - x_h) * lprime * s >= -1e-14),
                "linearized_horizon_x": linearized,
                "linearized_abs_error": abs(linearized - x),
            }
        )
    return rows


def tail_summary(rows: list[dict[str, float | bool]]) -> dict[str, float]:
    tail = [row for row in rows if row["d"] in set(float(v) for v in TAIL_D)]
    return {
        "max_tail_linearized_abs_error": float(max(float(row["linearized_abs_error"]) for row in tail)),
        "mean_tail_linearized_abs_error": float(np.mean([float(row["linearized_abs_error"]) for row in tail])),
        "max_abs_x_minus_horizon_tail": float(max(abs(float(row["x_minus_horizon"])) for row in tail)),
    }


def main() -> None:
    inverse_ode = load_json(INVERSE_ODE)
    horizon = horizon_data()
    rows = z12_rows(horizon)
    tail = tail_summary(rows)
    all_signs_consistent = all(bool(row["same_sign_as_strict_norm_after_lprime"]) for row in rows)

    report = {
        "status": "OPEN_COMPRESSION_HORIZON_TRACE_NO_ONTOLOGICAL_BRIDGE_THEOREM",
        "result_kind": "SCRATCH_PHASE_NORMALIZED_COMPRESSION_HORIZON_PROBE__NOT_A_THEOREM",
        "source_reports": {"inverse_branch_ode": str(INVERSE_ODE.relative_to(HERE.parents[1]))},
        "horizon": horizon,
        "z12_rows": rows,
        "tail_linearization_summary": tail,
        "sign_consistency_all_z12": all_signs_consistent,
        "candidate_ontological_reading": {
            "name": "finite_information_horizon_candidate",
            "supported_by_this_probe": bool(
                abs(float(horizon["horizon_formula_residual_vs_4_over_3"])) < 1e-14
                and abs(float(horizon["branch_horizon_residual_at_strict_zero"])) < 1e-12
                and all_signs_consistent
            ),
            "content": "The inverse-branch coordinate has a closed-form horizon x_h=4/3 where strict normalized output crosses zero; strict distances d=0..11 are compressed into a neighborhood of this finite horizon.",
            "not_yet_proven_because": "The horizon is derived from the legacy carrier zero and output matching, not from a strict-side information metric theorem.",
        },
        "upstream_replay": {
            "ode_candidate_supported": inverse_ode["candidate_ontological_reading"]["supported_by_this_probe"],
            "ode_derivation_status": inverse_ode["transport_equation"]["derivation_status"],
            "ode_compression_x11_over_11": inverse_ode["integration_check"]["compression_x11_over_11"],
        },
        "honest_interpretation": [
            "The monotone inverse branch has a computable finite horizon x_h=4/3 tied to the first zero of the normalized legacy carrier.",
            "The strict sign change maps to crossing this horizon, and high-d Z12 points stay close enough for first-order horizon linearization to be informative.",
            "This is an ontological-bridge candidate only if strict-side information geometry independently derives the same finite horizon and compression flow.",
        ],
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No strict-side derivation of the horizon x_h=4/3 is exported.",
            "No proof that the horizon is the nadsoliton information horizon is exported.",
            "No legacy physical-role transfer is licensed.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Search for a strict-side derivation of the finite horizon x_h=4/3 or falsify it as a legacy-carrier artifact; do not treat the horizon as ontological until that derivation exists.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch phase-normalized compression horizon probe\n\n"
        "Status: finite-horizon candidate; no ontological bridge theorem.\n\n"
        f"- Horizon `x_h={horizon['x_horizon']:.12f}` with residual vs `4/3` `{horizon['horizon_formula_residual_vs_4_over_3']:.3e}`.\n"
        f"- Strict zero at `d={horizon['strict_zero_d']:.12f}` maps to branch residual `{horizon['branch_horizon_residual_at_strict_zero']:.3e}` from the horizon.\n"
        f"- Tail linearization max abs error `{tail['max_tail_linearized_abs_error']:.3e}`; sign consistency on Z12 `{all_signs_consistent}`.\n"
        f"- Candidate ontology support flag `{report['candidate_ontological_reading']['supported_by_this_probe']}`: finite information horizon, not strict-side-derived.\n"
        "- No false pass: no kernel identity, no strict-side horizon derivation, no information-horizon theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
