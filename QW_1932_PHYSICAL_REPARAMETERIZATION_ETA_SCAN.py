#!/usr/bin/env python3
"""
QW-1932: Physical reparameterization scan for triad closure.

Reparameterized kernel family:
    K_eta(d) = cos(omega*d + phi) / (1 + beta * d**eta)

with eta=1 recovering the canonical legacy form.

Objective:
- test whether scale-dependent torsion exponent eta can remove the structural
  beta-interior vs fit/transfer trade-off observed in QW-1923/1924/1931.
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
OUT_JSON = ROOT / "report_qw1932_physical_reparameterization_eta_scan.json"
OUT_MD = ROOT / "RAPORT_QW1932_PHYSICAL_REPARAMETERIZATION_ETA_SCAN.md"


def wrap_phi(phi: float) -> float:
    return float((phi + math.pi) % (2.0 * math.pi) - math.pi)


def load_qw1739_module():
    path = ROOT / "QW_1739_SIGNED_DYNAMIC_MICROMODEL_DERIVATION.py"
    spec = importlib.util.spec_from_file_location("qw1739_mod_1932", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["qw1739_mod_1932"] = mod
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


def objective_eta(
    theta: Tuple[float, float, float],
    eta: float,
    d: np.ndarray,
    y: np.ndarray,
    weights: np.ndarray,
) -> float:
    omega, phi, beta = theta
    if omega < 0.02 or omega > 1.90 or beta <= 1e-4 or beta > 2.0:
        return float("inf")
    if eta < 1.0 or eta > 3.0:
        return float("inf")

    b = np.cos(omega * d + phi) / (1.0 + beta * (d ** eta))
    bb = float(np.dot(b, b))
    if bb <= 1e-12:
        return float("inf")

    a = (y @ b) / bb
    pred = a[:, None] * b[None, :]
    res = y - pred
    sse = float(np.sum(weights[:, None] * (res ** 2)))
    return float(sse)


def coordinate_refine_eta(
    start: Tuple[float, float, float],
    eta: float,
    d: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
) -> Tuple[Tuple[float, float, float], float]:
    omega, phi, beta = start
    cur = (float(np.clip(omega, 0.02, 1.90)), wrap_phi(phi), float(np.clip(beta, 1e-4, 2.0)))
    fcur = objective_eta(cur, eta, d, y, w)

    for so, sp, sb in [(0.30, 1.00, 0.50), (0.12, 0.45, 0.22), (0.045, 0.16, 0.08), (0.015, 0.07, 0.03)]:
        improved = True
        while improved:
            improved = False

            best = (cur, fcur)
            for om in np.linspace(max(0.02, cur[0] - so), min(1.90, cur[0] + so), 7):
                th = (float(om), cur[1], cur[2])
                f = objective_eta(th, eta, d, y, w)
                if f < best[1]:
                    best = (th, f)
            if best[1] < fcur:
                cur, fcur = best
                improved = True

            best = (cur, fcur)
            for ph in np.linspace(cur[1] - sp, cur[1] + sp, 7):
                th = (cur[0], wrap_phi(float(ph)), cur[2])
                f = objective_eta(th, eta, d, y, w)
                if f < best[1]:
                    best = (th, f)
            if best[1] < fcur:
                cur, fcur = best
                improved = True

            best = (cur, fcur)
            for be in np.linspace(max(1e-4, cur[2] - sb), min(2.0, cur[2] + sb), 7):
                th = (cur[0], cur[1], float(be))
                f = objective_eta(th, eta, d, y, w)
                if f < best[1]:
                    best = (th, f)
            if best[1] < fcur:
                cur, fcur = best
                improved = True

    return cur, fcur


def fit_eta(eta: float, d: np.ndarray, y: np.ndarray, w: np.ndarray) -> Dict[str, float]:
    starts = [
        (0.05625, 0.9425, 2.0),
        (0.3000, 2.4000, 0.5000),
        (0.1200, 0.5000, 0.5000),
        (0.2000, -0.8000, 0.3000),
    ]
    rng = np.random.default_rng(193201 + int(100 * eta))
    for _ in range(24):
        starts.append(
            (
                float(rng.uniform(0.02, 1.90)),
                float(rng.uniform(-math.pi, math.pi)),
                float(10 ** rng.uniform(-3.0, math.log10(2.0))),
            )
        )

    sols = []
    for st in starts:
        th, f = coordinate_refine_eta(st, eta, d, y, w)
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


def eval_external_holdout_eta(triad: Dict[str, float], eta: float, pta_path: Path) -> Dict[str, float]:
    df = pd.read_csv(pta_path)
    pair_id = df["pair_id"].astype(str).to_numpy()
    theta = df["theta_deg"].to_numpy(dtype=float)
    y = df["hxy"].to_numpy(dtype=float)

    tmin, tmax = float(np.min(theta)), float(np.max(theta))
    d_eff = 1.0 + 11.0 * (theta - tmin) / max(tmax - tmin, 1e-12)

    om = float(triad["omega"])
    ph = float(triad["phi"])
    be = float(triad["beta"])
    k = np.cos(om * d_eff + ph) / (1.0 + be * (d_eff ** eta))

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
    d, Y, W = build_profiles(n_grid=[64, 96, 128], seeds_per_n=14, dmax=24)
    n_obs = int(Y.size)

    p_primary = ROOT / "external_confirmatory_v2" / "beta_channel_true_external" / "beta_channel_pairs.csv"
    p_stress = ROOT / "external_confirmatory_v2" / "confirmatory_dataset_external_source_alpha6_1831cfg" / "pta_v2_pairs.csv"
    if not p_primary.exists() or not p_stress.exists():
        raise RuntimeError("Required external datasets not found.")

    eta_grid = [round(x, 1) for x in np.arange(1.0, 3.01, 0.2)]
    rows: List[Dict[str, object]] = []
    for eta in eta_grid:
        fit = fit_eta(float(eta), d, Y, W)
        tr_p = eval_external_holdout_eta(fit, float(eta), p_primary)
        tr_s = eval_external_holdout_eta(fit, float(eta), p_stress)
        rows.append(
            {
                "eta": float(eta),
                "fit": fit,
                "transfer_primary": tr_p,
                "transfer_stress": tr_s,
            }
        )

    base = next(r for r in rows if abs(float(r["eta"]) - 1.0) < 1e-12)
    obj0 = float(base["fit"]["objective"])
    cp0 = float(base["transfer_primary"]["corr"])
    gp0 = float(base["transfer_primary"]["rmse_gain"])
    cs0 = float(base["transfer_stress"]["corr"])
    gs0 = float(base["transfer_stress"]["rmse_gain"])

    for r in rows:
        obj = float(r["fit"]["objective"])
        eta = float(r["eta"])
        beta = float(r["fit"]["beta"])
        rel_loss = float((obj - obj0) / max(abs(obj0), 1e-12))
        cp = float(r["transfer_primary"]["corr"])
        gp = float(r["transfer_primary"]["rmse_gain"])
        cs = float(r["transfer_stress"]["corr"])
        gs = float(r["transfer_stress"]["rmse_gain"])

        # Approximate BIC with +1 complexity for eta-selection branch.
        bic = float(n_obs * np.log(max(obj / max(n_obs, 1), 1e-18)) + 4.0 * np.log(max(n_obs, 2)))
        bic0 = float(n_obs * np.log(max(obj0 / max(n_obs, 1), 1e-18)) + 3.0 * np.log(max(n_obs, 2)))
        delta_bic_vs_eta1 = float(bic - bic0)

        flags = {
            "beta_le_1p20": bool(beta <= 1.20),
            "rel_loss_le_0p35": bool(rel_loss <= 0.35),
            "primary_corr_ratio_ge_0p95": bool(cp / max(cp0, 1e-9) >= 0.95),
            "primary_gain_ratio_ge_1p00": bool(gp / max(gp0, 1e-9) >= 1.00),
            "stress_corr_ratio_ge_0p95": bool(cs / max(cs0, 1e-9) >= 0.95),
            "stress_gain_ratio_ge_1p00": bool(gs / max(gs0, 1e-9) >= 1.00),
            "omega_interior": bool(0.05 <= float(r["fit"]["omega"]) <= 1.50),
            "eta_not_grid_edge": bool(eta > min(eta_grid) and eta < max(eta_grid)),
            "delta_bic_le_10": bool(delta_bic_vs_eta1 <= 10.0),
        }

        r["comparison"] = {
            "rel_loss_vs_eta1": rel_loss,
            "primary_corr_ratio": float(cp / max(cp0, 1e-9)),
            "primary_gain_ratio": float(gp / max(gp0, 1e-9)),
            "stress_corr_ratio": float(cs / max(cs0, 1e-9)),
            "stress_gain_ratio": float(gs / max(gs0, 1e-9)),
            "delta_bic_vs_eta1": delta_bic_vs_eta1,
            "strict_flags": flags,
            "strict_all_pass": bool(all(flags.values())),
        }

    strict = [r for r in rows if r["comparison"]["strict_all_pass"]]
    if strict:
        selected = min(strict, key=lambda r: r["comparison"]["rel_loss_vs_eta1"])
        verdict = "PHYSICAL_REPARAMETERIZATION_STRICT_PASS"
        required_next = "INTEGRATE_REPARAM_BRANCH_IN_STAGE_B_GATE_AND_RETEST_IDENTIFIABILITY"
    else:
        selected = max(
            rows,
            key=lambda r: (
                int(r["comparison"]["strict_flags"]["beta_le_1p20"]),
                int(r["comparison"]["strict_flags"]["primary_gain_ratio_ge_1p00"]),
                int(r["comparison"]["strict_flags"]["stress_gain_ratio_ge_1p00"]),
                -r["comparison"]["rel_loss_vs_eta1"],
            ),
        )
        verdict = "PHYSICAL_REPARAMETERIZATION_PARTIAL_OR_FAIL"
        required_next = "EXTEND_ETA_GRID_OR_REFORMULATE_TORSION_TERM_WITHIN_MICRODYNAMICS"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "eta_grid": eta_grid,
        "eta1_baseline": base,
        "results": rows,
        "strict_pass_count": int(len(strict)),
        "selected": selected,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    sel = selected
    comp = sel["comparison"]
    lines = [
        "# RAPORT QW-1932: PHYSICAL REPARAMETERIZATION ETA SCAN",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- strict_pass_count: {out['strict_pass_count']}",
        "",
        "## Selected",
        f"- eta: {sel['eta']}",
        f"- omega: {sel['fit']['omega']:.6f}",
        f"- phi: {sel['fit']['phi']:.6f}",
        f"- beta: {sel['fit']['beta']:.6f}",
        f"- rel_loss_vs_eta1: {comp['rel_loss_vs_eta1']:.4f}",
        f"- primary corr/gain ratios: {comp['primary_corr_ratio']:.4f}/{comp['primary_gain_ratio']:.4f}",
        f"- stress corr/gain ratios: {comp['stress_corr_ratio']:.4f}/{comp['stress_gain_ratio']:.4f}",
        f"- delta_bic_vs_eta1: {comp['delta_bic_vs_eta1']:.4f}",
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

    print(f"[QW-1932] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1932] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1932] verdict={verdict} strict_pass_count={out['strict_pass_count']}")


if __name__ == "__main__":
    main()

