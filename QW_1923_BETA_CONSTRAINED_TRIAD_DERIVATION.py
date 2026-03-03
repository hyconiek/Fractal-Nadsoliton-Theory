#!/usr/bin/env python3
"""
QW-1923: Integrate selected beta observable as explicit constraint in triad derivation.

Uses QW-1922 selected observable signal (B7 intervention effect) to build an external
beta prior and injects it into no-ansatz triad fitting on extended micromodel profiles.

Outputs:
- unconstrained vs constrained triad fit,
- fit degradation check,
- boundary-pressure comparison.
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1923_beta_constrained_triad_derivation.json"
OUT_MD = ROOT / "RAPORT_QW1923_BETA_CONSTRAINED_TRIAD_DERIVATION.md"


def wrap_phi(phi: float) -> float:
    return float((phi + math.pi) % (2.0 * math.pi) - math.pi)


def load_qw1739_module():
    path = ROOT / "QW_1739_SIGNED_DYNAMIC_MICROMODEL_DERIVATION.py"
    spec = importlib.util.spec_from_file_location("qw1739_mod_1923", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["qw1739_mod_1923"] = mod
    spec.loader.exec_module(mod)
    return mod


def build_profiles(n_grid: List[int], seeds_per_n: int, dmax: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    mod = load_qw1739_module()
    rows = []
    for n in n_grid:
        for sidx in range(seeds_per_n):
            seed = 192000 + 1000 * n + sidx
            cfg = mod.MicroCfg(n_nodes=n, seed=seed)
            rng = np.random.default_rng(seed)
            theta, q = mod.build_micro_state(cfg, rng)
            w = mod.build_signed_dynamic_matrix(cfg, theta, q, rng)
            g = mod.effective_path_sum(w, lam=cfg.dyn_lambda, kmax=cfg.dyn_kmax)
            y = mod.profile_from_matrix(g, dmax=min(int(dmax), n // 2))
            rows.append(np.array(y, dtype=float))

    Y = np.array(rows, dtype=float)
    d = np.arange(1, Y.shape[1] + 1, dtype=float)
    W = np.array([1.0 / max(float(np.var(y)), 1e-5) for y in Y], dtype=float)
    return d, Y, W


def moving_average(x: np.ndarray, win: int) -> np.ndarray:
    w = max(3, int(win) | 1)
    pad = w // 2
    xp = np.pad(x, (pad, pad), mode="edge")
    ker = np.ones(w, dtype=float) / float(w)
    y = np.convolve(xp, ker, mode="same")
    return y[pad : pad + len(x)]


def robust_scale(x: np.ndarray) -> float:
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    return max(1e-9, 1.4826 * mad)


def b7_effect_for_vector(theta: np.ndarray, hxy: np.ndarray) -> float:
    order = np.argsort(theta)
    y = hxy[order]

    trend = moving_average(y, win=41)
    resid = y - trend
    abs_resid = np.abs(resid)
    b0 = moving_average(abs_resid, win=31)

    t = (theta - float(np.min(theta))) / max(float(np.max(theta) - np.min(theta)), 1e-12)
    hxy_beta = np.clip(hxy / (1.0 + 0.80 * t), 0.0, 1.0)
    yb = hxy_beta[order]
    trend_b = moving_average(yb, win=41)
    resid_b = yb - trend_b
    bb = moving_average(np.abs(resid_b), win=31)

    s = robust_scale(b0)
    return float(abs(np.median(bb) - np.median(b0)) / s)


def calibrate_beta_prior_from_qw1922(n_rep: int = 1000) -> Dict[str, float]:
    d1922 = json.loads((ROOT / "report_qw1922_beta_observable_blind_external_intervention.json").read_text(encoding="utf-8"))
    e_obs_primary = float(d1922["datasets"]["primary"]["holdout"]["effect_beta"])
    e_obs_stress = float(d1922["datasets"]["stress"]["holdout"]["effect_beta"])
    e_obs = float(np.median([e_obs_primary, e_obs_stress]))

    # Theta support from primary external dataset.
    pta = pd.read_csv(ROOT / "external_confirmatory_v2/confirmatory_dataset_external_source_rebuild_v2_1831cfg/pta_v2_pairs.csv")
    theta = pta["theta_deg"].to_numpy(dtype=float)
    t = (theta - float(np.min(theta))) / max(float(np.max(theta) - np.min(theta)), 1e-12)

    rng = np.random.default_rng(192301)
    beta_grid = np.linspace(0.005, 0.60, 80)
    eff_med = []

    for beta in beta_grid:
        vals = []
        for _ in range(max(80, n_rep // 10)):
            omega = float(np.clip(rng.normal(math.pi / 4.0, 0.08), 0.35, 1.10))
            phi = float(np.clip(rng.normal(math.pi / 6.0 if rng.random() < 0.5 else -math.pi / 6.0, 0.20), -1.4, 1.4))
            amp = float(np.clip(rng.normal(1.0, 0.15), 0.4, 1.8))
            noise = rng.normal(0.0, 0.03, size=len(theta))

            hxy = amp * np.cos(omega * (1.0 + 11.0 * t) + phi) / (1.0 + float(beta) * (1.0 + 11.0 * t))
            hxy = np.clip(0.40 + 0.20 * hxy + noise, 0.0, 1.0)

            vals.append(b7_effect_for_vector(theta, hxy))
        eff_med.append(float(np.median(vals)))

    eff_med_arr = np.array(eff_med, dtype=float)
    idx = int(np.argmin(np.abs(eff_med_arr - e_obs)))
    beta_target = float(beta_grid[idx])

    # Local scale from near-target region.
    close = np.abs(eff_med_arr - e_obs) <= 0.15
    if np.sum(close) >= 4:
        b_lo = float(np.quantile(beta_grid[close], 0.20))
        b_hi = float(np.quantile(beta_grid[close], 0.80))
    else:
        b_lo = max(0.005, beta_target - 0.08)
        b_hi = min(0.60, beta_target + 0.08)

    beta_scale = max(0.03, 0.5 * (b_hi - b_lo))

    return {
        "effect_obs_primary": e_obs_primary,
        "effect_obs_stress": e_obs_stress,
        "effect_obs_median": e_obs,
        "beta_target": beta_target,
        "beta_interval": [b_lo, b_hi],
        "beta_scale": beta_scale,
    }


def objective(
    theta: Tuple[float, float, float],
    d: np.ndarray,
    y: np.ndarray,
    weights: np.ndarray,
    beta_target: float | None = None,
    beta_scale: float | None = None,
    lambda_beta: float = 0.0,
) -> float:
    omega, phi, beta = theta
    if omega < 0.02 or omega > 1.90 or beta <= 1e-4 or beta > 2.0:
        return float("inf")

    b = np.cos(omega * d + phi) / (1.0 + beta * d)
    bb = float(np.dot(b, b))
    if bb <= 1e-12:
        return float("inf")

    a = (y @ b) / bb
    pred = a[:, None] * b[None, :]
    res = y - pred
    sse = float(np.sum(weights[:, None] * (res ** 2)))

    if beta_target is not None and beta_scale is not None and lambda_beta > 0.0:
        z = (beta - beta_target) / max(beta_scale, 1e-9)
        sse += float(lambda_beta) * float(z * z)

    return float(sse)


def coordinate_refine(
    start: Tuple[float, float, float],
    d: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    beta_target: float | None,
    beta_scale: float | None,
    lambda_beta: float,
) -> Tuple[Tuple[float, float, float], float]:
    omega, phi, beta = start
    omega = float(np.clip(omega, 0.02, 1.90))
    phi = wrap_phi(phi)
    beta = float(np.clip(beta, 1e-4, 2.0))

    cur = (omega, phi, beta)
    fcur = objective(cur, d, y, w, beta_target=beta_target, beta_scale=beta_scale, lambda_beta=lambda_beta)

    step_plan = [
        (0.30, 1.00, 0.50),
        (0.12, 0.45, 0.22),
        (0.045, 0.16, 0.08),
        (0.015, 0.07, 0.03),
    ]

    for so, sp, sb in step_plan:
        improved = True
        while improved:
            improved = False

            best = (cur, fcur)
            for om in np.linspace(max(0.02, cur[0] - so), min(1.90, cur[0] + so), 9):
                th = (float(om), cur[1], cur[2])
                f = objective(th, d, y, w, beta_target=beta_target, beta_scale=beta_scale, lambda_beta=lambda_beta)
                if f < best[1]:
                    best = (th, f)
            if best[1] < fcur:
                cur, fcur = best
                improved = True

            best = (cur, fcur)
            for ph in np.linspace(cur[1] - sp, cur[1] + sp, 9):
                th = (cur[0], wrap_phi(float(ph)), cur[2])
                f = objective(th, d, y, w, beta_target=beta_target, beta_scale=beta_scale, lambda_beta=lambda_beta)
                if f < best[1]:
                    best = (th, f)
            if best[1] < fcur:
                cur, fcur = best
                improved = True

            best = (cur, fcur)
            for be in np.linspace(max(1e-4, cur[2] - sb), min(2.0, cur[2] + sb), 9):
                th = (cur[0], cur[1], float(be))
                f = objective(th, d, y, w, beta_target=beta_target, beta_scale=beta_scale, lambda_beta=lambda_beta)
                if f < best[1]:
                    best = (th, f)
            if best[1] < fcur:
                cur, fcur = best
                improved = True

    return cur, fcur


def fit_global(
    d: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    beta_target: float | None,
    beta_scale: float | None,
    lambda_beta: float,
) -> Dict[str, object]:
    starts = [
        (0.118, 0.475, 0.300),
        (0.088, 0.890, 0.400),
        (0.75, 0.50, 0.02),
        (1.10, -0.30, 0.20),
        (0.20, -0.90, 0.80),
    ]
    rng = np.random.default_rng(192302)
    for _ in range(22):
        starts.append(
            (
                float(rng.uniform(0.02, 1.90)),
                float(rng.uniform(-math.pi, math.pi)),
                float(10 ** rng.uniform(-3.0, math.log10(2.0))),
            )
        )

    sols = []
    for st in starts:
        th, f = coordinate_refine(st, d, y, w, beta_target=beta_target, beta_scale=beta_scale, lambda_beta=lambda_beta)
        sols.append((th, f))
    sols = sorted(sols, key=lambda x: x[1])

    th, f = sols[0]
    return {
        "optimum": {"omega": float(th[0]), "phi": float(th[1]), "beta": float(th[2]), "objective": float(f)},
        "top_solutions": [{"theta": [float(x) for x in a], "objective": float(b)} for a, b in sols[:10]],
    }


def main() -> None:
    d, Y, W = build_profiles(n_grid=[64, 96, 128], seeds_per_n=14, dmax=24)

    prior = calibrate_beta_prior_from_qw1922(n_rep=1000)

    fit_un = fit_global(d, Y, W, beta_target=None, beta_scale=None, lambda_beta=0.0)
    fit_con = fit_global(
        d,
        Y,
        W,
        beta_target=float(prior["beta_target"]),
        beta_scale=float(prior["beta_scale"]),
        lambda_beta=12.0,
    )

    obj_un = float(fit_un["optimum"]["objective"])
    obj_con = float(fit_con["optimum"]["objective"])
    rel_loss = float((obj_con - obj_un) / max(abs(obj_un), 1e-12))

    beta_un = float(fit_un["optimum"]["beta"])
    beta_con = float(fit_con["optimum"]["beta"])

    # Pragmatic boundary relief criterion.
    beta_boundary_relief = bool(beta_un >= 1.80 and beta_con <= 1.20)
    fit_preserved = bool(rel_loss <= 0.25)

    if beta_boundary_relief and fit_preserved:
        verdict = "BETA_CONSTRAINED_TRIAD_DERIVATION_PARTIAL_SUCCESS"
        next_step = "REPEAT_STAGE_B_GATE_WITH_CONSTRAINED_TRIAD_AND_BLIND_EXTERNAL_RECHECK"
    elif beta_boundary_relief:
        verdict = "BETA_CONSTRAINED_TRIAD_DERIVATION_TRADEOFF_HIGH"
        next_step = "TUNE_CONSTRAINT_WEIGHT_AND_RETEST_TRANSFER"
    else:
        verdict = "BETA_CONSTRAINED_TRIAD_DERIVATION_NO_BOUNDARY_RELIEF"
        next_step = "DESIGN_TRUE_EXTERNAL_INTERVENTION_DATA_COLLECTION_FOR_BETA_CHANNEL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "config": {
            "profiles": {"n_grid": [64, 96, 128], "seeds_per_n": 14, "dmax": 24},
            "constraint_lambda_beta": 12.0,
        },
        "beta_prior_from_qw1922": prior,
        "fit_unconstrained": fit_un,
        "fit_beta_constrained": fit_con,
        "comparison": {
            "objective_unconstrained": obj_un,
            "objective_constrained": obj_con,
            "relative_objective_increase": rel_loss,
            "beta_unconstrained": beta_un,
            "beta_constrained": beta_con,
            "beta_boundary_relief": bool(beta_boundary_relief),
            "fit_preserved": bool(fit_preserved),
        },
        "verdict": verdict,
        "required_next_step": next_step,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1923: BETA-CONSTRAINED TRIAD DERIVATION",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Required next step: `{next_step}`",
        "",
        "## Beta Prior (from QW-1922)",
        f"- effect_obs_primary: {prior['effect_obs_primary']:.4f}",
        f"- effect_obs_stress: {prior['effect_obs_stress']:.4f}",
        f"- beta_target: {prior['beta_target']:.4f}",
        f"- beta_interval: [{prior['beta_interval'][0]:.4f}, {prior['beta_interval'][1]:.4f}]",
        f"- beta_scale: {prior['beta_scale']:.4f}",
        "",
        "## Fit Comparison",
        f"- unconstrained beta: {beta_un:.4f}",
        f"- constrained beta: {beta_con:.4f}",
        f"- unconstrained objective: {obj_un:.4f}",
        f"- constrained objective: {obj_con:.4f}",
        f"- relative objective increase: {rel_loss:.4f}",
        f"- beta_boundary_relief: {beta_boundary_relief}",
        f"- fit_preserved: {fit_preserved}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1923] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1923] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1923] verdict={verdict}")


if __name__ == "__main__":
    main()
