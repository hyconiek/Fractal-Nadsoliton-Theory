#!/usr/bin/env python3
"""
QW-1931: Strict triad feasibility frontier under external-beta constraints.

Goal:
- test whether one triad can satisfy simultaneously:
  1) interior beta (non-boundary),
  2) limited micromodel fit loss,
  3) preserved blind external transfer.

This script provides a hard feasibility verdict and a minimal-relaxation map.
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
OUT_JSON = ROOT / "report_qw1931_strict_triad_feasibility_frontier.json"
OUT_MD = ROOT / "RAPORT_QW1931_STRICT_TRIAD_FEASIBILITY_FRONTIER.md"


def wrap_phi(phi: float) -> float:
    return float((phi + math.pi) % (2.0 * math.pi) - math.pi)


def load_qw1739_module():
    path = ROOT / "QW_1739_SIGNED_DYNAMIC_MICROMODEL_DERIVATION.py"
    spec = importlib.util.spec_from_file_location("qw1739_mod_1931", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["qw1739_mod_1931"] = mod
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
    res = y - pred
    sse = float(np.sum(weights[:, None] * (res ** 2)))

    if lambda_beta > 0.0:
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
    rng = np.random.default_rng(193101)
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


def split_index(pair_id: str, k: int = 2) -> int:
    import hashlib

    h = hashlib.sha256(pair_id.encode("utf-8")).hexdigest()
    return int(h[-8:], 16) % k


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
    return {"corr": corr, "rmse_gain": gain}


def main() -> None:
    d1923 = json.loads((ROOT / "report_qw1923_beta_constrained_triad_derivation.json").read_text(encoding="utf-8"))
    prior = d1923["beta_prior_from_qw1922"]
    beta_target = float(prior["beta_target"])
    beta_scale = float(prior["beta_scale"])

    d, Y, W = build_profiles(n_grid=[64, 96, 128], seeds_per_n=14, dmax=24)

    p_primary = ROOT / "external_confirmatory_v2" / "beta_channel_true_external" / "beta_channel_pairs.csv"
    p_stress = ROOT / "external_confirmatory_v2" / "confirmatory_dataset_external_source_alpha6_1831cfg" / "pta_v2_pairs.csv"
    if not p_primary.exists():
        raise RuntimeError(f"Primary true external pairs not found: {p_primary}")
    if not p_stress.exists():
        raise RuntimeError(f"Stress pairs not found: {p_stress}")

    lambda_grid = [0.0, 0.25, 0.50, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0]
    rows: List[Dict[str, object]] = []
    for lam in lambda_grid:
        fit = fit_for_lambda(d, Y, W, beta_target, beta_scale, float(lam))
        tr_p = eval_external_holdout(fit, p_primary)
        tr_s = eval_external_holdout(fit, p_stress)
        rows.append(
            {
                "lambda_beta": float(lam),
                "fit": fit,
                "transfer_primary": tr_p,
                "transfer_stress": tr_s,
            }
        )

    base = rows[0]
    obj0 = float(base["fit"]["objective"])
    cp0 = float(base["transfer_primary"]["corr"])
    gp0 = float(base["transfer_primary"]["rmse_gain"])
    cs0 = float(base["transfer_stress"]["corr"])
    gs0 = float(base["transfer_stress"]["rmse_gain"])

    strict_candidates = []
    for r in rows:
        obj = float(r["fit"]["objective"])
        beta = float(r["fit"]["beta"])
        rel_loss = float((obj - obj0) / max(abs(obj0), 1e-12))
        cp = float(r["transfer_primary"]["corr"])
        gp = float(r["transfer_primary"]["rmse_gain"])
        cs = float(r["transfer_stress"]["corr"])
        gs = float(r["transfer_stress"]["rmse_gain"])

        flags = {
            "beta_le_1p20": bool(beta <= 1.20),
            "rel_loss_le_0p35": bool(rel_loss <= 0.35),
            "primary_corr_ratio_ge_0p90": bool(cp / max(cp0, 1e-9) >= 0.90),
            "primary_gain_ratio_ge_0p90": bool(gp / max(gp0, 1e-9) >= 0.90),
            "stress_corr_ratio_ge_0p80": bool(cs / max(cs0, 1e-9) >= 0.80),
            "stress_gain_ratio_ge_0p80": bool(gs / max(gs0, 1e-9) >= 0.80),
            "omega_interior": bool(0.05 <= float(r["fit"]["omega"]) <= 1.50),
        }
        r["comparison"] = {
            "rel_loss_vs_unconstrained": rel_loss,
            "primary_corr_ratio": float(cp / max(cp0, 1e-9)),
            "primary_gain_ratio": float(gp / max(gp0, 1e-9)),
            "stress_corr_ratio": float(cs / max(cs0, 1e-9)),
            "stress_gain_ratio": float(gs / max(gs0, 1e-9)),
            "strict_flags": flags,
            "strict_all_pass": bool(all(flags.values())),
        }
        if r["comparison"]["strict_all_pass"]:
            strict_candidates.append(r)

    if strict_candidates:
        best = min(strict_candidates, key=lambda r: r["comparison"]["rel_loss_vs_unconstrained"])
        verdict = "STRICT_TRIAD_FEASIBILITY_PASS"
        required_next = "RUN_STAGE_B_DERIVATIONAL_GATE_WITH_STRICT_TRIAD"
    else:
        feasible_partial = [
            r
            for r in rows
            if r["comparison"]["strict_flags"]["beta_le_1p20"]
            and r["comparison"]["strict_flags"]["primary_corr_ratio_ge_0p90"]
            and r["comparison"]["strict_flags"]["primary_gain_ratio_ge_0p90"]
        ]
        min_rel_partial = None
        if feasible_partial:
            min_rel_partial = float(min(r["comparison"]["rel_loss_vs_unconstrained"] for r in feasible_partial))

        best = max(
            rows,
            key=lambda r: (
                int(r["comparison"]["strict_flags"]["beta_le_1p20"]),
                int(r["comparison"]["strict_flags"]["primary_corr_ratio_ge_0p90"]),
                int(r["comparison"]["strict_flags"]["primary_gain_ratio_ge_0p90"]),
                -r["comparison"]["rel_loss_vs_unconstrained"],
            ),
        )
        verdict = "STRICT_TRIAD_FEASIBILITY_FAIL"
        required_next = "REFORMULATE_MICROMODEL_OR_PARAMETERIZATION_TO_REMOVE_STRUCTURAL_TRADEOFF"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "beta_prior": {
            "beta_target": beta_target,
            "beta_scale": beta_scale,
        },
        "lambda_grid_results": rows,
        "strict_pass_count": int(len(strict_candidates)),
        "selected": best,
        "minimal_rel_loss_if_beta_and_primary_transfer_hold": (
            None
            if "min_rel_partial" not in locals() or min_rel_partial is None
            else float(min_rel_partial)
        ),
        "verdict": verdict,
        "required_next_step": required_next,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    sel = out["selected"]
    comp = sel["comparison"]
    lines = [
        "# RAPORT QW-1931: STRICT TRIAD FEASIBILITY FRONTIER",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- strict_pass_count: {out['strict_pass_count']}",
        "",
        "## Selected Candidate",
        f"- lambda_beta: {sel['lambda_beta']}",
        f"- omega: {sel['fit']['omega']:.6f}",
        f"- phi: {sel['fit']['phi']:.6f}",
        f"- beta: {sel['fit']['beta']:.6f}",
        f"- rel_loss_vs_unconstrained: {comp['rel_loss_vs_unconstrained']:.4f}",
        f"- primary corr/gain ratios: {comp['primary_corr_ratio']:.4f} / {comp['primary_gain_ratio']:.4f}",
        f"- stress corr/gain ratios: {comp['stress_corr_ratio']:.4f} / {comp['stress_gain_ratio']:.4f}",
        "",
        "## Strict Flags (selected)",
    ]
    for k, v in comp["strict_flags"].items():
        lines.append(f"- {k}: {v}")
    lines.extend(
        [
            "",
            "## Required Next Step",
            f"- {required_next}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1931] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1931] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1931] verdict={verdict} strict_pass_count={out['strict_pass_count']}")


if __name__ == "__main__":
    main()

