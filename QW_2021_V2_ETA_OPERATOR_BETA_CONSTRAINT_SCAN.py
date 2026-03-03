#!/usr/bin/env python3
"""
QW-2021: V2 eta-operator scan under beta constraint.

Goal:
- test minimal structural repair of triad operator:
    K(d) = cos(omega*d + phi) / (1 + beta * d^eta)
- check whether beta-pressure can be satisfied without catastrophic fit loss.

This is a strict continuation after QW-2020 (lambda-only tuning failed).
"""

from __future__ import annotations

import importlib.util
import json
import math
import hashlib
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2021_v2_eta_operator_beta_constraint_scan.json"
OUT_MD = ROOT / "RAPORT_QW2021_V2_ETA_OPERATOR_BETA_CONSTRAINT_SCAN.md"


def wrap_phi(phi: float) -> float:
    return float((phi + math.pi) % (2.0 * math.pi) - math.pi)


def load_qw1739_module():
    path = ROOT / "QW_1739_SIGNED_DYNAMIC_MICROMODEL_DERIVATION.py"
    spec = importlib.util.spec_from_file_location("qw1739_mod_2021", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["qw1739_mod_2021"] = mod
    spec.loader.exec_module(mod)
    return mod


def build_profiles(n_grid: List[int], seeds_per_n: int, dmax: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    mod = load_qw1739_module()
    rows = []
    for n in n_grid:
        for sidx in range(seeds_per_n):
            seed = 202100 + 1000 * n + sidx
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


def calibrate_beta_prior_from_qw2017(n_rep: int = 1000) -> Dict[str, float]:
    d2017 = json.loads((ROOT / "report_qw2017_v2_beta_observable_blind_external_intervention.json").read_text(encoding="utf-8"))
    e_obs_primary = float(d2017["datasets"]["primary"]["holdout"]["effect_beta"])
    e_obs_stress = float(d2017["datasets"]["stress"]["holdout"]["effect_beta"])
    e_obs = float(np.median([e_obs_primary, e_obs_stress]))

    pta = pd.read_csv(ROOT / "external_confirmatory_v2/beta_channel_true_external_v2/beta_channel_pairs.csv")
    theta = pta["theta_deg"].to_numpy(dtype=float)
    t = (theta - float(np.min(theta))) / max(float(np.max(theta) - np.min(theta)), 1e-12)

    rng = np.random.default_rng(202101)
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


def kernel_fn(d: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d + phi) / (1.0 + beta * (d**eta))


def objective(
    theta: Tuple[float, float, float, float],
    d: np.ndarray,
    y: np.ndarray,
    weights: np.ndarray,
    beta_target: float | None,
    beta_scale: float | None,
    lambda_beta: float,
) -> float:
    omega, phi, beta, eta = theta
    if omega < 0.02 or omega > 1.90 or beta <= 1e-4 or beta > 2.0 or eta < 0.40 or eta > 1.80:
        return float("inf")

    b = kernel_fn(d, omega=omega, phi=phi, beta=beta, eta=eta)
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
    start: Tuple[float, float, float, float],
    d: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    beta_target: float | None,
    beta_scale: float | None,
    lambda_beta: float,
) -> Tuple[Tuple[float, float, float, float], float]:
    omega, phi, beta, eta = start
    cur = (
        float(np.clip(omega, 0.02, 1.90)),
        wrap_phi(phi),
        float(np.clip(beta, 1e-4, 2.0)),
        float(np.clip(eta, 0.40, 1.80)),
    )
    fcur = objective(cur, d, y, w, beta_target, beta_scale, lambda_beta)

    step_plan = [
        (0.30, 1.00, 0.50, 0.35),
        (0.12, 0.45, 0.22, 0.16),
        (0.045, 0.16, 0.08, 0.06),
        (0.015, 0.07, 0.03, 0.02),
    ]

    for so, sp, sb, se in step_plan:
        improved = True
        while improved:
            improved = False

            best = (cur, fcur)
            for om in np.linspace(max(0.02, cur[0] - so), min(1.90, cur[0] + so), 9):
                th = (float(om), cur[1], cur[2], cur[3])
                f = objective(th, d, y, w, beta_target, beta_scale, lambda_beta)
                if f < best[1]:
                    best = (th, f)
            if best[1] < fcur:
                cur, fcur = best
                improved = True

            best = (cur, fcur)
            for ph in np.linspace(cur[1] - sp, cur[1] + sp, 9):
                th = (cur[0], wrap_phi(float(ph)), cur[2], cur[3])
                f = objective(th, d, y, w, beta_target, beta_scale, lambda_beta)
                if f < best[1]:
                    best = (th, f)
            if best[1] < fcur:
                cur, fcur = best
                improved = True

            best = (cur, fcur)
            for be in np.linspace(max(1e-4, cur[2] - sb), min(2.0, cur[2] + sb), 9):
                th = (cur[0], cur[1], float(be), cur[3])
                f = objective(th, d, y, w, beta_target, beta_scale, lambda_beta)
                if f < best[1]:
                    best = (th, f)
            if best[1] < fcur:
                cur, fcur = best
                improved = True

            best = (cur, fcur)
            for et in np.linspace(max(0.40, cur[3] - se), min(1.80, cur[3] + se), 9):
                th = (cur[0], cur[1], cur[2], float(et))
                f = objective(th, d, y, w, beta_target, beta_scale, lambda_beta)
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
        (0.118, 0.475, 0.300, 1.00),
        (0.088, 0.890, 0.400, 1.00),
        (0.75, 0.50, 0.02, 1.00),
        (1.10, -0.30, 0.20, 1.00),
        (0.20, -0.90, 0.80, 1.00),
        (0.20, -0.90, 0.50, 0.75),
        (0.20, -0.90, 0.50, 1.35),
    ]

    rng = np.random.default_rng(202102)
    for _ in range(14):
        starts.append(
            (
                float(rng.uniform(0.02, 1.90)),
                float(rng.uniform(-math.pi, math.pi)),
                float(10 ** rng.uniform(-3.0, math.log10(2.0))),
                float(rng.uniform(0.40, 1.80)),
            )
        )

    sols = []
    for st in starts:
        th, f = coordinate_refine(st, d, y, w, beta_target, beta_scale, lambda_beta)
        sols.append((th, f))
    sols = sorted(sols, key=lambda x: x[1])

    th, f = sols[0]
    return {
        "optimum": {
            "omega": float(th[0]),
            "phi": float(th[1]),
            "beta": float(th[2]),
            "eta": float(th[3]),
            "objective": float(f),
        },
        "top_solutions": [
            {
                "theta": [float(x) for x in a],
                "objective": float(b),
            }
            for a, b in sols[:10]
        ],
    }


def eval_external_holdout(triad: Dict[str, float], pta_path: Path) -> Dict[str, float]:
    df = pd.read_csv(pta_path)
    pair_id = df["pair_id"].astype(str).to_numpy()
    theta = df["theta_deg"].to_numpy(dtype=float)
    y = df["hxy"].to_numpy(dtype=float)

    tmin, tmax = float(np.min(theta)), float(np.max(theta))
    d_eff = 1.0 + 11.0 * (theta - tmin) / max(tmax - tmin, 1e-12)

    om = float(triad["omega"])
    ph = float(triad["phi"])
    be = float(triad["beta"])
    et = float(triad["eta"])
    k = kernel_fn(d_eff, omega=om, phi=ph, beta=be, eta=et)

    fold = np.array([int(hashlib.sha256(x.encode("utf-8")).hexdigest()[-8:], 16) % 2 for x in pair_id], dtype=int)
    disc = fold == 0
    hold = fold == 1

    X = np.column_stack([k[disc], np.ones(int(np.sum(disc)), dtype=float)])
    coef, *_ = np.linalg.lstsq(X, y[disc], rcond=None)
    a, b = float(coef[0]), float(coef[1])

    yhat = a * k + b
    yt = y[hold]
    yp = yhat[hold]

    corr = float(np.corrcoef(yt, yp)[0, 1])
    rmse = float(np.sqrt(np.mean((yt - yp) ** 2)))
    base = float(np.mean(y[disc]))
    rmse0 = float(np.sqrt(np.mean((yt - base) ** 2)))
    gain = float((rmse0 - rmse) / max(rmse0, 1e-12))

    return {
        "corr": corr,
        "rmse_gain": gain,
    }


def main() -> None:
    prior = calibrate_beta_prior_from_qw2017(n_rep=1000)
    d, Y, W = build_profiles(n_grid=[64, 96, 128], seeds_per_n=14, dmax=24)

    p_primary = ROOT / "external_confirmatory_v2/beta_channel_true_external_v2/beta_channel_pairs.csv"

    lambda_grid = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0]

    rows: List[Dict[str, object]] = []
    for lam in lambda_grid:
        fit = fit_global(d, Y, W, prior["beta_target"], prior["beta_scale"], float(lam))
        tr = eval_external_holdout(fit["optimum"], p_primary)
        rows.append({"lambda_beta": float(lam), "fit": fit["optimum"], "transfer_primary": tr})

    base = rows[0]
    obj0 = float(base["fit"]["objective"])
    corr0 = float(base["transfer_primary"]["corr"])
    gain0 = float(base["transfer_primary"]["rmse_gain"])

    for r in rows:
        obj = float(r["fit"]["objective"])
        beta = float(r["fit"]["beta"])
        corr = float(r["transfer_primary"]["corr"])
        gain = float(r["transfer_primary"]["rmse_gain"])

        rel_loss = float((obj - obj0) / max(abs(obj0), 1e-12))
        corr_ratio = float(corr / max(corr0, 1e-9))
        gain_ratio = float(gain / max(gain0, 1e-9))

        fit_score = float(max(0.0, 1.0 - rel_loss / 0.30))
        beta_score = float(max(0.0, min(1.0, (2.0 - beta) / 1.0)))
        transfer_score = float(max(0.0, min(1.0, 0.5 * corr_ratio + 0.5 * gain_ratio)))

        # small regularization preference for eta close to 1 (minimal extension)
        eta = float(r["fit"]["eta"])
        eta_score = float(max(0.0, 1.0 - abs(eta - 1.0) / 0.8))

        total = 0.35 * fit_score + 0.25 * beta_score + 0.25 * transfer_score + 0.15 * eta_score

        r["comparison"] = {
            "relative_objective_increase": rel_loss,
            "corr_ratio_vs_unconstrained": corr_ratio,
            "gain_ratio_vs_unconstrained": gain_ratio,
            "fit_score": fit_score,
            "beta_score": beta_score,
            "transfer_score": transfer_score,
            "eta_score": eta_score,
            "total_score": float(total),
        }

    best = max(rows, key=lambda r: r["comparison"]["total_score"])

    flags = {
        "beta_below_1p2": bool(float(best["fit"]["beta"]) <= 1.2),
        "relative_loss_le_0p35": bool(float(best["comparison"]["relative_objective_increase"]) <= 0.35),
        "corr_ratio_ge_0p85": bool(float(best["comparison"]["corr_ratio_vs_unconstrained"]) >= 0.85),
        "gain_ratio_ge_0p85": bool(float(best["comparison"]["gain_ratio_vs_unconstrained"]) >= 0.85),
    }

    # Check if ANY lambda solves all hard constraints.
    any_full = False
    first_full = None
    for r in rows:
        ok = (
            float(r["fit"]["beta"]) <= 1.2
            and float(r["comparison"]["relative_objective_increase"]) <= 0.35
            and float(r["comparison"]["corr_ratio_vs_unconstrained"]) >= 0.85
            and float(r["comparison"]["gain_ratio_vs_unconstrained"]) >= 0.85
        )
        if ok:
            any_full = True
            first_full = r
            break

    if any_full:
        verdict = "ETA_OPERATOR_STRUCTURAL_REPAIR_PASS"
        required_next = "RERUN_STAGE_B_GATE_WITH_ETA_OPERATOR_FROZEN"
        selected = first_full
    elif sum(1 for v in flags.values() if v) >= 3:
        verdict = "ETA_OPERATOR_STRUCTURAL_REPAIR_PARTIAL"
        required_next = "ESCALATE_TO_NEXT_ORDER_OPERATOR_BASIS"
        selected = best
    else:
        verdict = "ETA_OPERATOR_STRUCTURAL_REPAIR_FAIL"
        required_next = "ESCALATE_TO_INFORMATION_COUPLED_OPERATOR_EXTENSION"
        selected = best

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "beta_prior": prior,
        "lambda_grid_results": rows,
        "selected": selected,
        "pass_flags_selected": flags,
        "any_full_pass": bool(any_full),
        "verdict": verdict,
        "required_next_step": required_next,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2021: V2 ETA-OPERATOR BETA-CONSTRAINT SCAN",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Selected lambda_beta: {selected['lambda_beta']}",
        "",
        "## Selected Solution",
        f"- omega: {selected['fit']['omega']:.6f}",
        f"- phi: {selected['fit']['phi']:.6f}",
        f"- beta: {selected['fit']['beta']:.6f}",
        f"- eta: {selected['fit']['eta']:.6f}",
        f"- relative_objective_increase: {selected['comparison']['relative_objective_increase']:.4f}",
        f"- corr_ratio_vs_unconstrained: {selected['comparison']['corr_ratio_vs_unconstrained']:.4f}",
        f"- gain_ratio_vs_unconstrained: {selected['comparison']['gain_ratio_vs_unconstrained']:.4f}",
        f"- total_score: {selected['comparison']['total_score']:.4f}",
        "",
        "## Pass Flags (Selected)",
        f"- beta_below_1p2: {flags['beta_below_1p2']}",
        f"- relative_loss_le_0p35: {flags['relative_loss_le_0p35']}",
        f"- corr_ratio_ge_0p85: {flags['corr_ratio_ge_0p85']}",
        f"- gain_ratio_ge_0p85: {flags['gain_ratio_ge_0p85']}",
        f"- any_full_pass: {any_full}",
        "",
        "## Required Next Step",
        f"- {required_next}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2021] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2021] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2021] verdict={verdict} selected_lambda={selected['lambda_beta']}")


if __name__ == "__main__":
    main()
