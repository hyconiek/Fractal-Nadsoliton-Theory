#!/usr/bin/env python3
"""
QW-1924: Tune beta-constraint weight and retest external transfer.

Next step from QW-1923:
TUNE_CONSTRAINT_WEIGHT_AND_RETEST_TRANSFER
"""

from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1924_lambda_tuning_and_transfer_retest.json"
OUT_MD = ROOT / "RAPORT_QW1924_LAMBDA_TUNING_AND_TRANSFER_RETEST.md"


def wrap_phi(phi: float) -> float:
    return float((phi + math.pi) % (2.0 * math.pi) - math.pi)


def split_index(pair_id: str, k: int = 2) -> int:
    h = hashlib.sha256(pair_id.encode("utf-8")).hexdigest()
    return int(h[-8:], 16) % k


def load_1923_inputs() -> Dict[str, float]:
    d1923 = json.loads((ROOT / "report_qw1923_beta_constrained_triad_derivation.json").read_text(encoding="utf-8"))
    prior = d1923["beta_prior_from_qw1922"]
    return {
        "beta_target": float(prior["beta_target"]),
        "beta_scale": float(prior["beta_scale"]),
    }


def load_profiles() -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    # Reuse deterministic profile package from QW-1923 by rerunning lightweight generator copy.
    import importlib.util
    import sys

    path = ROOT / "QW_1739_SIGNED_DYNAMIC_MICROMODEL_DERIVATION.py"
    spec = importlib.util.spec_from_file_location("qw1739_mod_1924", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["qw1739_mod_1924"] = mod
    spec.loader.exec_module(mod)

    rows = []
    for n in [64, 96, 128]:
        for sidx in range(14):
            seed = 192000 + 1000 * n + sidx
            cfg = mod.MicroCfg(n_nodes=n, seed=seed)
            rng = np.random.default_rng(seed)
            theta, q = mod.build_micro_state(cfg, rng)
            w = mod.build_signed_dynamic_matrix(cfg, theta, q, rng)
            g = mod.effective_path_sum(w, lam=cfg.dyn_lambda, kmax=cfg.dyn_kmax)
            y = mod.profile_from_matrix(g, dmax=min(24, n // 2))
            rows.append(np.array(y, dtype=float))

    Y = np.array(rows, dtype=float)
    d = np.arange(1, Y.shape[1] + 1, dtype=float)
    W = np.array([1.0 / max(float(np.var(y)), 1e-5) for y in Y], dtype=float)
    return d, Y, W


def objective(
    theta: Tuple[float, float, float],
    d: np.ndarray,
    y: np.ndarray,
    weights: np.ndarray,
    beta_target: float,
    beta_scale: float,
    lambda_beta: float,
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
    sse = float(np.sum(weights[:, None] * ((y - pred) ** 2)))

    z = (beta - beta_target) / max(beta_scale, 1e-9)
    sse += float(lambda_beta) * float(z * z)
    return float(sse)


def coordinate_refine(
    start: Tuple[float, float, float],
    d: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    beta_target: float,
    beta_scale: float,
    lambda_beta: float,
) -> Tuple[Tuple[float, float, float], float]:
    omega, phi, beta = start
    cur = (float(np.clip(omega, 0.02, 1.90)), wrap_phi(phi), float(np.clip(beta, 1e-4, 2.0)))
    fcur = objective(cur, d, y, w, beta_target, beta_scale, lambda_beta)

    for so, sp, sb in [(0.30, 1.00, 0.50), (0.12, 0.45, 0.22), (0.045, 0.16, 0.08), (0.015, 0.07, 0.03)]:
        improved = True
        while improved:
            improved = False

            best = (cur, fcur)
            for om in np.linspace(max(0.02, cur[0] - so), min(1.90, cur[0] + so), 9):
                th = (float(om), cur[1], cur[2])
                f = objective(th, d, y, w, beta_target, beta_scale, lambda_beta)
                if f < best[1]:
                    best = (th, f)
            if best[1] < fcur:
                cur, fcur = best
                improved = True

            best = (cur, fcur)
            for ph in np.linspace(cur[1] - sp, cur[1] + sp, 9):
                th = (cur[0], wrap_phi(float(ph)), cur[2])
                f = objective(th, d, y, w, beta_target, beta_scale, lambda_beta)
                if f < best[1]:
                    best = (th, f)
            if best[1] < fcur:
                cur, fcur = best
                improved = True

            best = (cur, fcur)
            for be in np.linspace(max(1e-4, cur[2] - sb), min(2.0, cur[2] + sb), 9):
                th = (cur[0], cur[1], float(be))
                f = objective(th, d, y, w, beta_target, beta_scale, lambda_beta)
                if f < best[1]:
                    best = (th, f)
            if best[1] < fcur:
                cur, fcur = best
                improved = True

    return cur, fcur


def fit_for_lambda(
    d: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    beta_target: float,
    beta_scale: float,
    lambda_beta: float,
) -> Dict[str, float]:
    starts = [
        (0.118, 0.475, 0.300),
        (0.088, 0.890, 0.400),
        (0.75, 0.50, 0.02),
        (1.10, -0.30, 0.20),
        (0.20, -0.90, 0.80),
    ]
    rng = np.random.default_rng(192401)
    for _ in range(16):
        starts.append(
            (
                float(rng.uniform(0.02, 1.90)),
                float(rng.uniform(-math.pi, math.pi)),
                float(10 ** rng.uniform(-3.0, math.log10(2.0))),
            )
        )

    sols = []
    for st in starts:
        th, f = coordinate_refine(st, d, y, w, beta_target, beta_scale, lambda_beta)
        sols.append((th, f))
    sols = sorted(sols, key=lambda x: x[1])

    th, f = sols[0]
    return {
        "omega": float(th[0]),
        "phi": float(th[1]),
        "beta": float(th[2]),
        "objective": float(f),
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
    k = np.cos(om * d_eff + ph) / (1.0 + be * d_eff)

    fold = np.array([split_index(x, k=2) for x in pair_id], dtype=int)
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
    prior = load_1923_inputs()
    d, Y, W = load_profiles()

    p_primary = ROOT / "external_confirmatory_v2/confirmatory_dataset_external_source_rebuild_v2_1831cfg/pta_v2_pairs.csv"

    lambda_grid = [0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0]

    rows: List[Dict[str, object]] = []
    for lam in lambda_grid:
        fit = fit_for_lambda(d, Y, W, prior["beta_target"], prior["beta_scale"], lam)
        tr = eval_external_holdout(fit, p_primary)
        rows.append(
            {
                "lambda_beta": float(lam),
                "fit": fit,
                "transfer_primary": tr,
            }
        )

    base = rows[0]
    obj0 = float(base["fit"]["objective"])
    corr0 = float(base["transfer_primary"]["corr"])
    gain0 = float(base["transfer_primary"]["rmse_gain"])

    for r in rows:
        obj = float(r["fit"]["objective"])
        rel_loss = float((obj - obj0) / max(abs(obj0), 1e-12))
        beta = float(r["fit"]["beta"])
        corr = float(r["transfer_primary"]["corr"])
        gain = float(r["transfer_primary"]["rmse_gain"])

        fit_score = float(max(0.0, 1.0 - rel_loss / 0.30))
        beta_score = float(max(0.0, min(1.0, (2.0 - beta) / 1.0)))
        transfer_score = float(max(0.0, min(1.0, 0.5 * (corr / max(corr0, 1e-9)) + 0.5 * (gain / max(gain0, 1e-9)))))

        total = 0.40 * fit_score + 0.30 * beta_score + 0.30 * transfer_score

        r["comparison"] = {
            "relative_objective_increase": rel_loss,
            "corr_ratio_vs_unconstrained": float(corr / max(corr0, 1e-9)),
            "gain_ratio_vs_unconstrained": float(gain / max(gain0, 1e-9)),
            "fit_score": fit_score,
            "beta_score": beta_score,
            "transfer_score": transfer_score,
            "total_score": float(total),
        }

    best = max(rows, key=lambda r: r["comparison"]["total_score"])

    beta_best = float(best["fit"]["beta"])
    rel_loss_best = float(best["comparison"]["relative_objective_increase"])
    corr_ratio = float(best["comparison"]["corr_ratio_vs_unconstrained"])
    gain_ratio = float(best["comparison"]["gain_ratio_vs_unconstrained"])

    flags = {
        "beta_below_1p2": bool(beta_best <= 1.2),
        "relative_loss_le_0p35": bool(rel_loss_best <= 0.35),
        "corr_ratio_ge_0p85": bool(corr_ratio >= 0.85),
        "gain_ratio_ge_0p85": bool(gain_ratio >= 0.85),
    }

    if all(flags.values()):
        verdict = "LAMBDA_TUNING_BALANCED_PASS"
        required_next = "RERUN_STAGE_B_GATE_WITH_SELECTED_LAMBDA"
    elif sum(1 for v in flags.values() if v) >= 3:
        verdict = "LAMBDA_TUNING_PARTIAL"
        if (float(best["lambda_beta"]) <= 0.0 + 1e-12) and (not flags["beta_below_1p2"]):
            required_next = "COLLECT_TRUE_EXTERNAL_INTERVENTION_DATA_FOR_BETA_CHANNEL"
        else:
            required_next = "RERUN_STAGE_B_GATE_WITH_SELECTED_LAMBDA"
    else:
        verdict = "LAMBDA_TUNING_FAIL"
        required_next = "COLLECT_TRUE_EXTERNAL_INTERVENTION_DATA_FOR_BETA_CHANNEL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "beta_prior": prior,
        "lambda_grid_results": rows,
        "selected": best,
        "pass_flags": flags,
        "verdict": verdict,
        "required_next_step": required_next,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1924: LAMBDA TUNING AND TRANSFER RETEST",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Selected lambda_beta: {best['lambda_beta']}",
        "",
        "## Selected Solution",
        f"- omega: {best['fit']['omega']:.6f}",
        f"- phi: {best['fit']['phi']:.6f}",
        f"- beta: {best['fit']['beta']:.6f}",
        f"- relative_objective_increase: {rel_loss_best:.4f}",
        f"- corr_ratio_vs_unconstrained: {corr_ratio:.4f}",
        f"- gain_ratio_vs_unconstrained: {gain_ratio:.4f}",
        f"- total_score: {best['comparison']['total_score']:.4f}",
        "",
        "## Pass Flags",
        f"- beta_below_1p2: {flags['beta_below_1p2']}",
        f"- relative_loss_le_0p35: {flags['relative_loss_le_0p35']}",
        f"- corr_ratio_ge_0p85: {flags['corr_ratio_ge_0p85']}",
        f"- gain_ratio_ge_0p85: {flags['gain_ratio_ge_0p85']}",
        "",
        "## Required Next Step",
        f"- {out['required_next_step']}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1924] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1924] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1924] verdict={verdict} selected_lambda={best['lambda_beta']}")


if __name__ == "__main__":
    main()
