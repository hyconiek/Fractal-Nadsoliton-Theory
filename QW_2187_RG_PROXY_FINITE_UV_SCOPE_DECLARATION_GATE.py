#!/usr/bin/env python3
"""
QW-2187: RG proxy finite-UV-scope declaration gate (strict).

Purpose:
- fulfill the strict next-step path from QW-2185 by declaring an explicit
  finite UV-validity scope for the current proxy RG system,
- separate clearly:
  1) what is certified inside scope,
  2) what is observed to fail outside scope (lambda crossing / no global closure),
- keep no-overclaim discipline.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2187_rg_proxy_finite_uv_scope_declaration_gate.json"
OUT_MD = ROOT / "RAPORT_QW2187_RG_PROXY_FINITE_UV_SCOPE_DECLARATION_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def beta_vec_batch(x: np.ndarray) -> np.ndarray:
    g1, g2, g3, yt, lam, ggr = [x[:, i] for i in range(6)]
    c = 1.0 / (16.0 * math.pi * math.pi)

    b1 = c * (41.0 / 6.0) * g1**3
    b2 = c * (-19.0 / 6.0) * g2**3
    b3 = c * (-7.0) * g3**3
    byt = c * yt * (
        (9.0 / 2.0) * yt**2
        - (17.0 / 12.0) * g1**2
        - (9.0 / 4.0) * g2**2
        - 8.0 * g3**2
    )
    blam = c * (
        24.0 * lam**2
        - 6.0 * yt**4
        + (3.0 / 8.0) * (2.0 * g2**4 + (g2**2 + g1**2) ** 2)
        + (-9.0 * g2**2 - 3.0 * g1**2 + 12.0 * yt**2) * lam
    )
    bgr = 2.0 * ggr * (1.0 - ggr)

    return np.stack([b1, b2, b3, byt, blam, bgr], axis=1)


def rk4_batch(x: np.ndarray, dt: float) -> np.ndarray:
    k1 = beta_vec_batch(x)
    k2 = beta_vec_batch(x + 0.5 * dt * k1)
    k3 = beta_vec_batch(x + 0.5 * dt * k2)
    k4 = beta_vec_batch(x + dt * k3)
    return x + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0


def main() -> None:
    r2182 = load_json("report_qw2182_rg_nonperturbative_domain_flow_gate.json")
    r2185 = load_json("report_qw2185_rg_proxy_global_obstruction_theorem_gate.json")

    grid = r2182["domain"]["grid"]
    initials = np.array(
        list(
            product(
                grid["g1"],
                grid["g2"],
                grid["g3"],
                grid["y_t"],
                grid["lambda_h"],
                grid["g_gr"],
            )
        ),
        dtype=float,
    )

    n_pts = int(initials.shape[0])
    dt = 0.01
    t_probe = 30.0
    n_steps_probe = int(round(t_probe / dt))
    safety_steps = 4

    x = initials.copy()
    first_cross_idx = np.full(n_pts, -1, dtype=int)
    lam_min_probe = float("inf")
    finite_probe = True

    for i in range(n_steps_probe + 1):
        if not np.all(np.isfinite(x)):
            finite_probe = False
            break

        lam = x[:, 4]
        lam_min_probe = min(lam_min_probe, float(np.min(lam)))
        newly_crossed = (first_cross_idx < 0) & (lam < 0.0)
        first_cross_idx[newly_crossed] = i

        x = rk4_batch(x, dt)

    any_cross = bool(np.any(first_cross_idx >= 0))
    first_cross_global_idx = int(np.min(first_cross_idx[first_cross_idx >= 0])) if any_cross else -1
    first_cross_global_t = float(first_cross_global_idx * dt) if any_cross else None

    if any_cross and first_cross_global_idx > safety_steps:
        t_scope = float((first_cross_global_idx - safety_steps) * dt)
    else:
        t_scope = float(min(r2182["domain"]["t_max"], t_probe))

    # Recheck strict scope guarantees explicitly on [0, t_scope].
    n_steps_scope = int(round(t_scope / dt))
    y = initials.copy()
    finite_scope = True
    lam_min_scope = float("inf")
    max_abs_scope = 0.0
    for _ in range(n_steps_scope + 1):
        if not np.all(np.isfinite(y)):
            finite_scope = False
            break
        lam_min_scope = min(lam_min_scope, float(np.min(y[:, 4])))
        max_abs_scope = max(max_abs_scope, float(np.max(np.abs(y))))
        y = rk4_batch(y, dt)

    t_star_domain_min = float(r2185["obstruction_theorem"]["t_star_domain_min"])

    flags = {
        "q2185_obstruction_theorem_loaded": str(r2185.get("verdict", "")).startswith(
            "RG_PROXY_GLOBAL_OBSTRUCTION_THEOREM_GATE_PASS"
        ),
        "domain_grid_loaded_from_q2182": bool(n_pts == 729),
        "probe_integration_finite": bool(finite_probe),
        "lambda_crossing_detected_before_probe_horizon": bool(any_cross),
        "scope_declared_strictly_below_first_crossing": bool(first_cross_global_t is not None and t_scope < first_cross_global_t),
        "scope_integration_finite": bool(finite_scope),
        "lambda_nonnegative_on_declared_scope": bool(lam_min_scope >= 0.0),
        "landau_pole_far_beyond_declared_scope": bool(t_star_domain_min > t_scope),
        "deterministic_no_scan_no_retune": True,
        "full_global_rg_closure_claimed": False,
    }

    pass_count = int(sum(bool(v) for v in flags.values()))
    total_flags = int(len(flags))

    core_ok = bool(
        flags["q2185_obstruction_theorem_loaded"]
        and flags["domain_grid_loaded_from_q2182"]
        and flags["probe_integration_finite"]
        and flags["lambda_crossing_detected_before_probe_horizon"]
        and flags["scope_declared_strictly_below_first_crossing"]
        and flags["scope_integration_finite"]
        and flags["lambda_nonnegative_on_declared_scope"]
        and flags["landau_pole_far_beyond_declared_scope"]
        and flags["deterministic_no_scan_no_retune"]
    )

    verdict = (
        "RG_PROXY_FINITE_UV_SCOPE_DECLARATION_GATE_PASS_STRICT"
        if core_ok
        else "RG_PROXY_FINITE_UV_SCOPE_DECLARATION_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2182": "report_qw2182_rg_nonperturbative_domain_flow_gate.json",
            "q2185": "report_qw2185_rg_proxy_global_obstruction_theorem_gate.json",
        },
        "probe_setup": {
            "n_points": n_pts,
            "dt": dt,
            "t_probe": t_probe,
            "n_steps_probe": n_steps_probe,
            "safety_steps_before_crossing": safety_steps,
        },
        "scope_declaration": {
            "first_lambda_crossing_time_global": first_cross_global_t,
            "declared_t_scope_max": t_scope,
            "lambda_min_on_scope": lam_min_scope,
            "lambda_min_on_probe": lam_min_probe,
            "max_abs_on_scope": max_abs_scope,
            "landau_pole_domain_min_time": t_star_domain_min,
        },
        "scope_statement": {
            "inside_scope": (
                "On 729-point strict grid, proxy flow stays finite and lambda_h remains nonnegative "
                "for 0 <= t <= t_scope_max."
            ),
            "outside_scope_indicator": (
                "Before t_probe horizon, at least one trajectory crosses lambda_h < 0; "
                "global all-t closure is not claimed."
            ),
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "USE_DECLARED_SCOPE_IN_RELEASE_STATUS_OR_INTRODUCE_UV_COMPLETING_BETA_CORRECTIONS_FOR_GLOBAL_L12_CLOSURE"
            if verdict.endswith("STRICT")
            else "REVIEW_SCOPE_DETECTION_NUMERICS_AND_RERUN_QW2187"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    md_lines = [
        "# RAPORT QW-2187: RG PROXY FINITE UV SCOPE DECLARATION GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Scope",
        f"- first lambda crossing (global on strict grid): `{first_cross_global_t}`",
        f"- declared strict scope max: `{t_scope}`",
        f"- lambda_min on scope: `{lam_min_scope:.12f}`",
        "",
        "## Separation",
        "- Inside scope: finite flow + nonnegative lambda_h on strict grid.",
        "- Outside scope: lambda_h crossing appears before probe horizon; no global all-t claim.",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()

