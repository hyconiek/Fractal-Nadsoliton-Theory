#!/usr/bin/env python3
"""
QW-1753: Orthogonal-triad joint inference.

Objective:
- Combine beta evidence from QW-1749 with omega/phi evidence from QW-1752.
- Quantify joint identifiability under weighted empirical posteriors.
- Keep strict separation between observables to avoid ansatz-driven coupling.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1753_orthogonal_triad_joint_inference.json"
OUT_MD = ROOT / "RAPORT_QW1753_ORTHOGONAL_TRIAD_JOINT_INFERENCE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def normalize_weights(w: np.ndarray) -> np.ndarray:
    z = np.array(w, dtype=float)
    z = np.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0)
    z = np.clip(z, 0.0, None)
    s = float(np.sum(z))
    if s <= 1e-15:
        return np.ones_like(z, dtype=float) / max(len(z), 1)
    return z / s


def effective_sample_size(w: np.ndarray) -> float:
    ww = normalize_weights(w)
    denom = float(np.sum(ww ** 2))
    if denom <= 1e-15:
        return 0.0
    return float(1.0 / denom)


def weighted_quantile(values: np.ndarray, weights: np.ndarray, probs: np.ndarray) -> np.ndarray:
    x = np.array(values, dtype=float)
    w = normalize_weights(np.array(weights, dtype=float))
    p = np.array(probs, dtype=float)

    idx = np.argsort(x)
    xs = x[idx]
    ws = w[idx]
    cw = np.cumsum(ws)
    out = np.interp(np.clip(p, 0.0, 1.0), cw, xs)
    return out


def wrap_phi(x: float) -> float:
    return float((x + math.pi) % (2.0 * math.pi) - math.pi)


def weighted_circular_mean(phi: np.ndarray, w: np.ndarray) -> float:
    ww = normalize_weights(w)
    s = float(np.sum(ww * np.sin(phi)))
    c = float(np.sum(ww * np.cos(phi)))
    return float(math.atan2(s, c))


def weighted_circular_std(phi: np.ndarray, w: np.ndarray) -> float:
    ww = normalize_weights(w)
    s = float(np.sum(ww * np.sin(phi)))
    c = float(np.sum(ww * np.cos(phi)))
    r = float(np.sqrt(s * s + c * c))
    r = max(r, 1e-12)
    return float(np.sqrt(-2.0 * np.log(r)))


def rankdata(x: np.ndarray) -> np.ndarray:
    order = np.argsort(x)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(len(x), dtype=float)
    return ranks


def spearman_corr(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3 or len(y) < 3:
        return float("nan")
    rx = rankdata(x)
    ry = rankdata(y)
    c = np.corrcoef(rx, ry)
    return float(c[0, 1])


def weighted_polyfit_1d(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> Tuple[float, float]:
    xx = np.array(x, dtype=float)
    yy = np.array(y, dtype=float)
    ww = normalize_weights(np.array(w, dtype=float))
    X = np.column_stack([xx, np.ones_like(xx)])
    W = np.diag(ww)
    XtWX = X.T @ W @ X
    XtWy = X.T @ W @ yy
    beta = np.linalg.pinv(XtWX) @ XtWy
    return float(beta[0]), float(beta[1])


def main() -> None:
    d1749 = load("report_qw1749_beta_orthogonal_observable.json")
    d1752 = load("report_qw1752_omega_orthogonal_observable.json")

    r49 = d1749.get("rows", [])
    if len(r49) < 10:
        raise RuntimeError("QW-1749 rows too short for weighted posterior construction.")

    beta = np.array([float(r["beta_hat"]) for r in r49], dtype=float)
    rmse_b = np.array([float(r["rmse_env"]) for r in r49], dtype=float)
    dsc = np.array([float(r["delta_beta_scramble"]) for r in r49], dtype=float)

    # Prefer low envelope RMSE and low scramble sensitivity.
    w_beta = (1.0 / np.clip(rmse_b ** 2, 1e-6, None)) * np.exp(-dsc / 0.02)
    w_beta = normalize_weights(w_beta)
    beta_q = weighted_quantile(beta, w_beta, np.array([0.025, 0.50, 0.975], dtype=float))
    beta_neff = effective_sample_size(w_beta)
    beta_boundary_rate = float(np.sum(w_beta * ((beta <= 0.0005) | (beta >= 0.23)).astype(float)))

    rows52 = d1752.get("rows", [])
    p = d1752.get("protocol", {}).get("good_filter", {})
    smin = float(p.get("snr_like_min", 4.0))
    rcmax = float(p.get("phase_rmse_complex_max", 0.45))
    rhmax = float(p.get("phase_rmse_hilbert_max", 0.55))

    good52 = [
        r
        for r in rows52
        if np.isfinite(float(r.get("omega_hat", float("nan"))))
        and float(r.get("snr_like", 0.0)) >= smin
        and float(r.get("phase_rmse_complex", 1e9)) <= rcmax
        and float(r.get("phase_rmse_hilbert", 1e9)) <= rhmax
    ]
    if len(good52) < 8:
        raise RuntimeError("QW-1752 good subset too short for orthogonal-triad inference.")

    omega = np.array([float(r["omega_hat"]) for r in good52], dtype=float)
    phi = np.array([float(r["phi_hat"]) for r in good52], dtype=float)
    snr = np.array([float(r["snr_like"]) for r in good52], dtype=float)
    rmse_c = np.array([float(r["phase_rmse_complex"]) for r in good52], dtype=float)
    rmse_h = np.array([float(r["phase_rmse_hilbert"]) for r in good52], dtype=float)
    dbe = np.array([float(r["delta_omega_beta_perturb"]) for r in good52], dtype=float)

    # Prefer high SNR, low phase fit errors, and low beta-perturb sensitivity.
    w_omega = (snr / np.clip((0.35 + rmse_c + rmse_h) ** 2, 1e-6, None)) * np.exp(-dbe / 0.08)
    w_omega = normalize_weights(w_omega)

    omega_q = weighted_quantile(omega, w_omega, np.array([0.025, 0.50, 0.975], dtype=float))
    omega_neff = effective_sample_size(w_omega)
    omega_boundary_rate = float(np.sum(w_omega * ((omega <= 0.10) | (omega >= 1.6)).astype(float)))
    phi_mean = weighted_circular_mean(phi, w_omega)
    phi_std = weighted_circular_std(phi, w_omega)

    rho_omega_dbe = spearman_corr(omega, dbe)
    slope_dbe_omega, intercept_dbe_omega = weighted_polyfit_1d(omega, dbe, w_omega)

    # Weighted empirical posterior sampling (separated observables).
    rng = np.random.default_rng(1753)
    n_mc = 20000
    idx_b = rng.choice(len(beta), size=n_mc, replace=True, p=w_beta)
    idx_o = rng.choice(len(omega), size=n_mc, replace=True, p=w_omega)
    b_mc = beta[idx_b]
    o_mc = omega[idx_o]
    p_mc = np.array([wrap_phi(v) for v in phi[idx_o]], dtype=float)

    joint_ci = {
        "beta": [float(np.quantile(b_mc, 0.025)), float(np.quantile(b_mc, 0.975))],
        "omega": [float(np.quantile(o_mc, 0.025)), float(np.quantile(o_mc, 0.975))],
        "phi": [float(np.quantile(p_mc, 0.025)), float(np.quantile(p_mc, 0.975))],
    }
    joint_nonboundary_rate = float(
        np.mean((b_mc > 0.0005) & (b_mc < 0.23) & (o_mc > 0.10) & (o_mc < 1.6))
    )

    p49 = d1749.get("pass_flags", {})
    p52 = d1752.get("pass_flags", {})

    pass_beta_obs = bool(p49.get("fit_quality")) and bool(p49.get("spread_control")) and bool(p49.get("phase_orthogonality"))
    pass_omega_obs = bool(p52.get("fit_quality")) and bool(p52.get("spread_control"))
    pass_info = (beta_neff >= 10.0) and (omega_neff >= 6.0) and (len(good52) >= 10)
    pass_orth = (
        np.isfinite(rho_omega_dbe)
        and abs(float(rho_omega_dbe)) <= 0.35
        and float(np.quantile(dbe, 0.9)) <= 0.10
        and abs(slope_dbe_omega) <= 0.20
    )
    pass_phi = phi_std <= 1.10
    pass_nonboundary = (
        beta_q[0] > 0.0005
        and omega_q[0] > 0.10
        and omega_q[2] < 1.6
        and joint_nonboundary_rate >= 0.75
    )

    if pass_beta_obs and pass_omega_obs and pass_info and pass_orth and pass_phi and pass_nonboundary:
        verdict = "ORTHOGONAL_TRIAD_IDENTIFIABILITY_SUPPORTED"
    elif pass_beta_obs and pass_omega_obs and pass_info and pass_orth:
        verdict = "ORTHOGONAL_TRIAD_IDENTIFIABILITY_PARTIAL"
    else:
        verdict = "ORTHOGONAL_TRIAD_IDENTIFIABILITY_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "beta_source": "QW-1749",
            "omega_phi_source": "QW-1752",
            "n_beta_rows": len(r49),
            "n_omega_rows_total": len(rows52),
            "n_omega_rows_good": len(good52),
        },
        "weighted_posteriors": {
            "beta": {
                "neff": beta_neff,
                "ci95": [float(beta_q[0]), float(beta_q[2])],
                "median": float(beta_q[1]),
                "boundary_rate": beta_boundary_rate,
            },
            "omega": {
                "neff": omega_neff,
                "ci95": [float(omega_q[0]), float(omega_q[2])],
                "median": float(omega_q[1]),
                "boundary_rate": omega_boundary_rate,
            },
            "phi": {
                "circular_mean": phi_mean,
                "circular_std": phi_std,
            },
        },
        "orthogonality_tests": {
            "spearman_omega_vs_delta_beta_perturb": rho_omega_dbe,
            "weighted_slope_delta_beta_perturb_vs_omega": slope_dbe_omega,
            "weighted_intercept_delta_beta_perturb_vs_omega": intercept_dbe_omega,
            "delta_beta_perturb_q90": float(np.quantile(dbe, 0.9)),
        },
        "joint_empirical_mc": {
            "n_samples": n_mc,
            "ci95": joint_ci,
            "joint_nonboundary_rate": joint_nonboundary_rate,
        },
        "pass_flags": {
            "beta_observable_quality": bool(pass_beta_obs),
            "omega_observable_quality": bool(pass_omega_obs),
            "information_content": bool(pass_info),
            "orthogonality": bool(pass_orth),
            "phi_stability": bool(pass_phi),
            "nonboundary_joint": bool(pass_nonboundary),
        },
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines: List[str] = [
        "# RAPORT QW-1753: ORTHOGONAL TRIAD JOINT INFERENCE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Input rows: beta={len(r49)}, omega_total={len(rows52)}, omega_good={len(good52)}",
        f"- beta median/CI95: {beta_q[1]:.6f} / [{beta_q[0]:.6f}, {beta_q[2]:.6f}]",
        f"- omega median/CI95: {omega_q[1]:.6f} / [{omega_q[0]:.6f}, {omega_q[2]:.6f}]",
        f"- phi circular std: {phi_std:.6f}",
        f"- rho_s(omega,delta_beta_perturb): {rho_omega_dbe:.6f}",
        f"- joint nonboundary rate: {joint_nonboundary_rate:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- beta_observable_quality: {pass_beta_obs}",
        f"- omega_observable_quality: {pass_omega_obs}",
        f"- information_content: {pass_info}",
        f"- orthogonality: {pass_orth}",
        f"- phi_stability: {pass_phi}",
        f"- nonboundary_joint: {pass_nonboundary}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1753] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1753] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
