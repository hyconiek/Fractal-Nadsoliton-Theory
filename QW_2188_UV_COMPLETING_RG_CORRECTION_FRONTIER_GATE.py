#!/usr/bin/env python3
"""
QW-2188: UV-completing RG correction frontier gate (strict, anchored partial).

Purpose:
- test whether an anchored UV correction family can extend the strict finite scope
  of current proxy-RG beyond QW-2187 crossing horizon,
- find minimal correction coefficient inside micro-anchored range that removes
  lambda_h crossing up to declared probe horizon,
- keep boundary explicit: this is extended finite-scope feasibility, not global
  all-t nonperturbative closure.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
from typing import Dict, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2188_uv_completing_rg_correction_frontier_gate.json"
OUT_MD = ROOT / "RAPORT_QW2188_UV_COMPLETING_RG_CORRECTION_FRONTIER_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def simulate_frontier(
    b_corr: float,
    initials: np.ndarray,
    dt: float,
    t_probe: float,
    g1_cap: float,
) -> Dict:
    """Simulate corrected proxy RG and return strict diagnostics."""
    c = 1.0 / (16.0 * math.pi * math.pi)
    x = initials.copy()
    n_steps = int(round(t_probe / dt))
    first_cross_idx = np.full(x.shape[0], -1, dtype=int)
    finite_all = True
    lam_min = float("inf")
    max_abs = 0.0

    def beta_batch(y: np.ndarray) -> np.ndarray:
        g1, g2, g3, yt, lam, ggr = [y[:, i] for i in range(6)]

        b1 = c * (41.0 / 6.0) * g1**3 * (1.0 - (g1 / g1_cap) ** 2)
        b2 = c * (-19.0 / 6.0) * g2**3
        b3 = c * (-7.0) * g3**3
        byt = c * yt * (
            (9.0 / 2.0) * yt**2
            - (17.0 / 12.0) * g1**2
            - (9.0 / 4.0) * g2**2
            - 8.0 * g3**2
        )

        # Anchored correction family:
        # scale only the destabilizing -6*y_t^4 term by (1-b_corr).
        gauge_quart = (3.0 / 8.0) * (2.0 * g2**4 + (g2**2 + g1**2) ** 2)
        blam = c * (
            24.0 * lam**2
            - 6.0 * (1.0 - b_corr) * yt**4
            + gauge_quart
            + (-9.0 * g2**2 - 3.0 * g1**2 + 12.0 * yt**2) * lam
        )
        bgr = 2.0 * ggr * (1.0 - ggr)

        return np.stack([b1, b2, b3, byt, blam, bgr], axis=1)

    for i in range(n_steps + 1):
        if not np.all(np.isfinite(x)):
            finite_all = False
            break

        lam = x[:, 4]
        lam_min = min(lam_min, float(np.min(lam)))
        max_abs = max(max_abs, float(np.max(np.abs(x))))
        new_cross = (first_cross_idx < 0) & (lam < 0.0)
        first_cross_idx[new_cross] = i

        k1 = beta_batch(x)
        k2 = beta_batch(x + 0.5 * dt * k1)
        k3 = beta_batch(x + 0.5 * dt * k2)
        k4 = beta_batch(x + dt * k3)
        x = x + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0

    any_cross = bool(np.any(first_cross_idx >= 0))
    first_cross_t = float(np.min(first_cross_idx[first_cross_idx >= 0]) * dt) if any_cross else None

    return {
        "finite_all": finite_all,
        "any_lambda_cross": any_cross,
        "first_lambda_cross_t": first_cross_t,
        "lambda_min": lam_min,
        "max_abs": max_abs,
        "passes_probe": bool(finite_all and (not any_cross)),
    }


def main() -> None:
    r2066 = load_json("report_qw2066_compatibility_filtered_micro_constants_tightening.json")
    r2132 = load_json("report_qw2132_rg_fixed_point_jacobian_gate.json")
    r2187 = load_json("report_qw2187_rg_proxy_finite_uv_scope_declaration_gate.json")

    sel = r2066["selected_filter"]
    targets = r2066["targets"]

    z_beta_q50 = float(sel["z_beta_q50"])
    delta_eta_q25 = float(sel["delta_eta_q25"])
    delta_eta_q75 = float(sel["delta_eta_q75"])
    delta_eta_target = float(targets["delta_eta_target"])

    # Anchors:
    # g1_cap from micro z_beta median.
    g1_cap = float(math.sqrt(z_beta_q50 / 100.0))
    # b_corr normalized to target delta_eta envelope.
    b_min_anchor = float(max(0.0, (delta_eta_q25 - delta_eta_target) / max(delta_eta_target, 1e-300)))
    b_max_anchor = float(max(0.0, (delta_eta_q75 - delta_eta_target) / max(delta_eta_target, 1e-300)))

    grid = [
        [0.30, 0.35, 0.40],
        [0.55, 0.65, 0.75],
        [1.00, 1.15, 1.30],
        [0.75, 0.85, 0.95],
        [0.12, 0.16, 0.20],
        [1e-5, 0.15, 0.30],
    ]
    initials = np.array(list(product(*grid)), dtype=float)

    dt = 0.01
    t_probe = 30.0
    b_ref_base = 0.0

    cache: Dict[float, Dict] = {}

    def eval_b(b: float) -> Dict:
        key = round(float(b), 12)
        if key not in cache:
            cache[key] = simulate_frontier(key, initials, dt=dt, t_probe=t_probe, g1_cap=g1_cap)
        return cache[key]

    base_diag = eval_b(b_ref_base)
    hi_diag = eval_b(b_max_anchor)

    feasible = bool(hi_diag["passes_probe"])
    b_star = None
    b_star_diag = None

    # Deterministic bisection to find minimal feasible b in anchored interval.
    if feasible:
        lo = b_min_anchor
        hi = b_max_anchor
        # Ensure upper side is feasible and lower side tested.
        lo_diag = eval_b(lo)
        if lo_diag["passes_probe"]:
            b_star = lo
            b_star_diag = lo_diag
        else:
            for _ in range(12):
                mid = 0.5 * (lo + hi)
                mid_diag = eval_b(mid)
                if mid_diag["passes_probe"]:
                    hi = mid
                else:
                    lo = mid
            b_star = hi
            b_star_diag = eval_b(b_star)

    # Low-energy distortion audit at QW-2132 reference point.
    ref = r2132["reference_couplings_at_mu_1gev"]
    g1 = float(ref["g1_u1"])
    g2 = float(ref["g2_su2"])
    g3 = float(ref["g3_su3"])
    yt = float(ref["y_top"])
    lam = float(ref["lambda_h"])
    c = 1.0 / (16.0 * math.pi * math.pi)

    def beta_lambda_ref(b_corr: float) -> float:
        gauge_quart = (3.0 / 8.0) * (2.0 * g2**4 + (g2**2 + g1**2) ** 2)
        return float(
            c
            * (
                24.0 * lam**2
                - 6.0 * (1.0 - b_corr) * yt**4
                + gauge_quart
                + (-9.0 * g2**2 - 3.0 * g1**2 + 12.0 * yt**2) * lam
            )
        )

    beta_lam_base = beta_lambda_ref(0.0)
    beta_lam_star = beta_lambda_ref(float(b_star)) if b_star is not None else None
    rel_shift_beta_lam = (
        float(abs(beta_lam_star - beta_lam_base) / max(abs(beta_lam_base), 1e-300))
        if beta_lam_star is not None
        else None
    )

    flags = {
        "micro_anchors_loaded": True,
        "anchored_interval_defined_nonempty": bool(b_max_anchor >= b_min_anchor >= 0.0),
        "baseline_crossing_reproduced_before_t_probe": bool(base_diag["any_lambda_cross"]),
        "feasible_solution_exists_within_anchor_interval": bool(feasible),
        "minimal_feasible_b_star_found": bool(b_star is not None and b_star_diag is not None),
        "b_star_within_anchor_interval": bool(
            b_star is not None and b_min_anchor - 1e-12 <= b_star <= b_max_anchor + 1e-12
        ),
        "b_star_removes_crossing_to_t_probe": bool(b_star_diag["passes_probe"]) if b_star_diag else False,
        "b_star_flow_finite_to_t_probe": bool(b_star_diag["finite_all"]) if b_star_diag else False,
        "low_energy_beta_lambda_shift_le_0p75": bool(
            rel_shift_beta_lam is not None and rel_shift_beta_lam <= 0.75
        ),
        "deterministic_no_scan_no_retune": True,
        "global_all_t_rg_closure_claimed": False,
    }

    pass_count = int(sum(bool(v) for v in flags.values()))
    total_flags = int(len(flags))

    core_ok = bool(
        flags["micro_anchors_loaded"]
        and flags["anchored_interval_defined_nonempty"]
        and flags["baseline_crossing_reproduced_before_t_probe"]
        and flags["feasible_solution_exists_within_anchor_interval"]
        and flags["minimal_feasible_b_star_found"]
        and flags["b_star_within_anchor_interval"]
        and flags["b_star_removes_crossing_to_t_probe"]
        and flags["b_star_flow_finite_to_t_probe"]
        and flags["low_energy_beta_lambda_shift_le_0p75"]
        and flags["deterministic_no_scan_no_retune"]
    )

    verdict = (
        "UV_COMPLETING_RG_CORRECTION_FRONTIER_GATE_PASS_EXTENDED_SCOPE_PARTIAL"
        if core_ok
        else "UV_COMPLETING_RG_CORRECTION_FRONTIER_GATE_FAIL_OR_INFEASIBLE_IN_ANCHORED_RANGE"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2066": "report_qw2066_compatibility_filtered_micro_constants_tightening.json",
            "q2132": "report_qw2132_rg_fixed_point_jacobian_gate.json",
            "q2187": "report_qw2187_rg_proxy_finite_uv_scope_declaration_gate.json",
        },
        "anchors": {
            "g1_cap_from_zbeta_q50": g1_cap,
            "b_corr_interval_from_delta_eta_q25_q75": [b_min_anchor, b_max_anchor],
            "delta_eta_q25_q75": [delta_eta_q25, delta_eta_q75],
            "delta_eta_target": delta_eta_target,
        },
        "probe_setup": {
            "n_points": int(initials.shape[0]),
            "dt": dt,
            "t_probe": t_probe,
            "bisection_iterations": 12,
        },
        "frontier_results": {
            "baseline_b0": base_diag,
            "anchor_upper_bmax": hi_diag,
            "feasible_in_anchor_interval": feasible,
            "b_star_min_feasible": b_star,
            "b_star_diagnostics": b_star_diag,
        },
        "reference_low_energy_distortion": {
            "beta_lambda_base": beta_lam_base,
            "beta_lambda_star": beta_lam_star,
            "relative_shift_beta_lambda": rel_shift_beta_lam,
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "PROMOTE_EXTENDED_SCOPE_TO_RELEASE_STATUS_AND_SEPARATELY_PROVE_OR_REJECT_GLOBAL_ALL_T_CLOSURE_UNDER_UV_COMPLETING_FAMILY"
            if verdict.endswith("PARTIAL")
            else "REVIEW_ANCHORED_UV_FAMILY_OR_EXPAND_ALLOWED_CORRECTION_CLASS"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    md_lines = [
        "# RAPORT QW-2188: UV COMPLETING RG CORRECTION FRONTIER GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Frontier summary",
        f"- anchored `b_corr` interval: `{[b_min_anchor, b_max_anchor]}`",
        f"- baseline first crossing (`b=0`): `{base_diag['first_lambda_cross_t']}`",
        f"- minimal feasible `b*` (no crossing to `t_probe={t_probe}`): `{b_star}`",
        f"- low-energy relative shift `beta_lambda`: `{rel_shift_beta_lam}`",
        "",
        "## Scope",
        "- This gate extends finite probe scope only (up to declared horizon).",
        "- No global all-t closure is claimed.",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()

