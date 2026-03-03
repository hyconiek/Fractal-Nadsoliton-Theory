#!/usr/bin/env python3
"""
QW-1750: Separated inference pipeline
Step 1: beta from QW-1749 (orthogonal observable),
Step 2: omega,phi from dynamic phase observables with beta fixed.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1750_separated_beta_then_omega_phi.json"
OUT_MD = ROOT / "RAPORT_QW1750_SEPARATED_BETA_THEN_OMEGA_PHI.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def wrap_phi(x: float) -> float:
    return float((x + math.pi) % (2 * math.pi) - math.pi)


def circ_mean(x: np.ndarray) -> float:
    return float(np.arctan2(np.mean(np.sin(x)), np.mean(np.cos(x))))


def circ_std(x: np.ndarray) -> float:
    r = np.sqrt(np.mean(np.sin(x)) ** 2 + np.mean(np.cos(x)) ** 2)
    if r <= 1e-12:
        return float("inf")
    return float(np.sqrt(-2 * np.log(r)))


def objective_2p(theta: Tuple[float, float], omega_phase: np.ndarray, omega_zero: np.ndarray, phi_obs: np.ndarray, w_phase: np.ndarray, w_zero: np.ndarray, w_phi: np.ndarray) -> float:
    omega, phi = theta
    if not (0.10 <= omega <= 1.6 and -math.pi <= phi <= math.pi):
        return float("inf")
    res = 0.0
    res += float(np.sum(w_phase * (omega_phase - omega) ** 2))
    res += float(np.sum(w_zero * (omega_zero - omega) ** 2))
    dphi = np.array([wrap_phi(float(x - phi)) for x in phi_obs], dtype=float)
    res += float(np.sum(w_phi * (dphi ** 2)))
    return res


def refine_2p(start: Tuple[float, float], om_p: np.ndarray, om_z: np.ndarray, ph: np.ndarray, wp: np.ndarray, wz: np.ndarray, wh: np.ndarray) -> Tuple[Tuple[float, float], float]:
    cur = (float(start[0]), wrap_phi(float(start[1])))
    fcur = objective_2p(cur, om_p, om_z, ph, wp, wz, wh)
    steps = [(0.20, 0.80), (0.08, 0.30), (0.03, 0.12), (0.012, 0.05)]
    for so, sp in steps:
        improved = True
        while improved:
            improved = False
            # omega
            best = (cur, fcur)
            for om in np.linspace(max(0.10, cur[0] - so), min(1.6, cur[0] + so), 9):
                th = (float(om), cur[1])
                f = objective_2p(th, om_p, om_z, ph, wp, wz, wh)
                if f < best[1]:
                    best = (th, f)
            if best[1] < fcur:
                cur, fcur = best
                improved = True

            # phi
            best = (cur, fcur)
            for phv in np.linspace(cur[1] - sp, cur[1] + sp, 9):
                th = (cur[0], wrap_phi(float(phv)))
                f = objective_2p(th, om_p, om_z, ph, wp, wz, wh)
                if f < best[1]:
                    best = (th, f)
            if best[1] < fcur:
                cur, fcur = best
                improved = True
    return cur, fcur


def hessian_2p(theta: Tuple[float, float], om_p: np.ndarray, om_z: np.ndarray, ph: np.ndarray, wp: np.ndarray, wz: np.ndarray, wh: np.ndarray) -> np.ndarray:
    x = np.array(theta, dtype=float)
    h = np.array([0.006, 0.02], dtype=float)
    f0 = objective_2p((x[0], x[1]), om_p, om_z, ph, wp, wz, wh)
    H = np.zeros((2, 2), dtype=float)

    for i in range(2):
        xp = x.copy(); xm = x.copy()
        xp[i] += h[i]; xm[i] -= h[i]
        if i == 1:
            xp[i] = wrap_phi(float(xp[i])); xm[i] = wrap_phi(float(xm[i]))
        fp = objective_2p((xp[0], xp[1]), om_p, om_z, ph, wp, wz, wh)
        fm = objective_2p((xm[0], xm[1]), om_p, om_z, ph, wp, wz, wh)
        H[i, i] = (fp - 2 * f0 + fm) / (h[i] ** 2)

    xpp = x.copy(); xpm = x.copy(); xmp = x.copy(); xmm = x.copy()
    xpp[0] += h[0]; xpp[1] += h[1]
    xpm[0] += h[0]; xpm[1] -= h[1]
    xmp[0] -= h[0]; xmp[1] += h[1]
    xmm[0] -= h[0]; xmm[1] -= h[1]
    xpp[1] = wrap_phi(float(xpp[1])); xpm[1] = wrap_phi(float(xpm[1]))
    xmp[1] = wrap_phi(float(xmp[1])); xmm[1] = wrap_phi(float(xmm[1]))
    fpp = objective_2p((xpp[0], xpp[1]), om_p, om_z, ph, wp, wz, wh)
    fpm = objective_2p((xpm[0], xpm[1]), om_p, om_z, ph, wp, wz, wh)
    fmp = objective_2p((xmp[0], xmp[1]), om_p, om_z, ph, wp, wz, wh)
    fmm = objective_2p((xmm[0], xmm[1]), om_p, om_z, ph, wp, wz, wh)
    v = (fpp - fpm - fmp + fmm) / (4 * h[0] * h[1])
    H[0, 1] = v; H[1, 0] = v
    return H


def bootstrap_2p(om_p: np.ndarray, om_z: np.ndarray, ph: np.ndarray, wp: np.ndarray, wz: np.ndarray, wh: np.ndarray, th0: Tuple[float, float], n_boot: int = 200, seed: int = 1750) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n = len(om_p)
    out = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        th, _ = refine_2p(th0, om_p[idx], om_z[idx], ph[idx], wp[idx], wz[idx], wh[idx])
        out.append(th)
    return np.array(out, dtype=float)


def fit_omega_phi_fixed_beta_response(y: np.ndarray, beta_fixed: float) -> Tuple[float, float, float, float]:
    d = np.arange(1, len(y) + 1, dtype=float)
    om_grid_coarse = np.linspace(0.10, 1.6, 101)
    ph_grid_coarse = np.linspace(-math.pi, math.pi, 121)
    best = {"omega": 0.5, "phi": 0.0, "amp": 0.0, "rmse": float("inf")}

    def scan(om_grid: np.ndarray, ph_grid: np.ndarray, cur_best: Dict[str, float]) -> Dict[str, float]:
        out = dict(cur_best)
        for om in om_grid:
            wd = om * d
            for ph in ph_grid:
                b = np.cos(wd + ph) / (1.0 + beta_fixed * d)
                bb = float(np.dot(b, b))
                if bb <= 1e-12:
                    continue
                a = float(np.dot(y, b) / bb)
                pred = a * b
                rmse = float(np.sqrt(np.mean((y - pred) ** 2)))
                if rmse < out["rmse"]:
                    out = {"omega": float(om), "phi": float(ph), "amp": a, "rmse": rmse}
        return out

    best = scan(om_grid_coarse, ph_grid_coarse, best)
    om_c, ph_c = best["omega"], best["phi"]
    om_grid_ref = np.linspace(max(0.10, om_c - 0.12), min(1.6, om_c + 0.12), 61)
    ph_grid_ref = np.linspace(ph_c - 0.45, ph_c + 0.45, 61)
    best = scan(om_grid_ref, ph_grid_ref, best)
    best["phi"] = wrap_phi(float(best["phi"]))
    return best["omega"], best["phi"], best["amp"], best["rmse"]


def main() -> None:
    d1749 = load("report_qw1749_beta_orthogonal_observable.json")
    d1746 = load("report_qw1746_dynamic_observables_derivation.json")

    beta_fix = float(d1749["summary"]["beta_median"])
    beta_ci = d1749["summary"]["beta_ci95_empirical"]
    beta_verdict = d1749["verdict"]

    rows = d1746.get("rows", [])
    # use relaxed informative subset from 1746
    good = [r for r in rows if bool(r.get("is_good_relaxed", False))]
    if len(good) < 8:
        raise RuntimeError("Too few relaxed rows from 1746 for separated inference.")

    om_p = np.array([r["omega_phase_hat"] for r in good], dtype=float)
    om_z = np.array([r["omega_zero_hat"] for r in good], dtype=float)
    ph = np.array([r["phi_phase_hat"] for r in good], dtype=float)
    wp = 1.0 / np.clip(np.array([r["phase_fit_rmse"] for r in good], dtype=float) ** 2, 1e-5, None)
    wz = 1.0 / np.clip((0.6 + np.array([r["phase_fit_rmse"] for r in good], dtype=float)) ** 2, 1e-5, None)
    wh = wp.copy()

    starts = [
        (float(np.median(om_p)), float(circ_mean(ph))),
        (float(np.median(om_z)), float(circ_mean(ph) + 0.4)),
        (0.65, -0.2),
        (0.35, 0.8),
    ]
    sols = sorted([refine_2p(st, om_p, om_z, ph, wp, wz, wh) for st in starts], key=lambda x: x[1])
    th, fbest = sols[0]

    H = hessian_2p(th, om_p, om_z, ph, wp, wz, wh)
    H = 0.5 * (H + H.T)
    H = np.nan_to_num(H, nan=1e6, posinf=1e6, neginf=-1e6)
    eig = np.linalg.eigvalsh(H)
    min_e = float(np.min(eig))
    reg = 0.0
    if min_e <= 1e-8:
        reg = (1e-8 - min_e) + 1e-6
    Hr = H + reg * np.eye(2)
    cond = float(np.linalg.cond(Hr))
    cov = np.linalg.inv(Hr)
    std = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    ci_h = [float(1.96 * std[0]), float(1.96 * std[1])]
    corr = cov / np.sqrt(np.outer(np.diag(cov), np.diag(cov)))
    corr = np.clip(corr, -1.0, 1.0)
    max_corr = float(np.max(np.abs(corr - np.eye(2))))

    boot = bootstrap_2p(om_p, om_z, ph, wp, wz, wh, th, n_boot=200, seed=1750)
    bo, bp = boot[:, 0], boot[:, 1]
    ci_boot = {
        "omega": [float(np.quantile(bo, 0.025)), float(np.quantile(bo, 0.975))],
        "phi": [float(np.quantile(bp, 0.025)), float(np.quantile(bp, 0.975))],
    }

    width_omega = ci_boot["omega"][1] - ci_boot["omega"][0]
    phi_std = float(circ_std(bp))

    # Sensitivity over beta CI from 1749 using direct response fitting:
    # y(d) ~= A*cos(omega*d+phi)/(1+beta*d), beta fixed.
    beta_lo, beta_hi = float(beta_ci[0]), float(beta_ci[1])
    if beta_hi <= beta_lo:
        beta_hi = beta_lo + 1e-6
    beta_grid = np.linspace(beta_lo, beta_hi, 9)
    resp_rows = []
    for r in sorted(good, key=lambda x: float(x.get("snr_like", 0.0)), reverse=True):
        pr = r.get("response_real", {})
        if pr:
            vec = np.array([float(pr[str(i)]) for i in range(1, 13)], dtype=float)
            if np.max(np.abs(vec)) > 1e-12:
                vec = vec / np.max(np.abs(vec))
            resp_rows.append(vec)
        if len(resp_rows) >= 8:
            break

    sweep = []
    for bfix in beta_grid:
        if resp_rows:
            omega_fit = []
            phi_fit = []
            rmse_fit = []
            for yrow in resp_rows:
                omr, phr, _a, rmr = fit_omega_phi_fixed_beta_response(yrow, float(bfix))
                omega_fit.append(omr)
                phi_fit.append(phr)
                rmse_fit.append(rmr)
            omega_med_b = float(np.median(np.array(omega_fit, dtype=float)))
            phi_cmean_b = float(circ_mean(np.array(phi_fit, dtype=float)))
            rmse_med_b = float(np.median(np.array(rmse_fit, dtype=float)))
        else:
            omega_med_b = float("nan")
            phi_cmean_b = float("nan")
            rmse_med_b = float("nan")

        on_boundary = (not np.isfinite(omega_med_b)) or (omega_med_b <= 0.105 or omega_med_b >= 1.595)
        sweep.append(
            {
                "beta_fixed": float(bfix),
                "omega_hat_median_response": omega_med_b,
                "phi_hat_cmean_response": phi_cmean_b,
                "rmse_median_response": rmse_med_b,
                "omega_on_boundary": bool(on_boundary),
            }
        )
    boundary_rate = float(np.mean([1.0 if x["omega_on_boundary"] else 0.0 for x in sweep])) if sweep else 1.0

    pass_beta_quality = beta_verdict in ("BETA_ORTHOGONAL_OBSERVABLE_SUPPORTED", "BETA_ORTHOGONAL_OBSERVABLE_PARTIAL")
    pass_omega_width = width_omega <= 0.30
    pass_phi_width = phi_std <= 1.0
    pass_cond = cond <= 1e5
    pass_corr = max_corr <= 0.90
    pass_nonboundary = 0.10 < th[0] < 1.6
    pass_sweep = boundary_rate <= 0.40

    if pass_beta_quality and pass_omega_width and pass_phi_width and pass_nonboundary and (pass_cond or pass_corr) and pass_sweep:
        verdict = "SEPARATED_PIPELINE_IDENTIFIABILITY_MODERATE_OR_BETTER"
    elif pass_beta_quality and pass_nonboundary and (pass_omega_width or pass_phi_width):
        verdict = "SEPARATED_PIPELINE_PARTIAL"
    else:
        verdict = "SEPARATED_PIPELINE_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "beta_from_1749": {
            "beta_fixed": beta_fix,
            "beta_ci95_empirical": beta_ci,
            "beta_source_verdict": beta_verdict,
        },
        "n_good_rows_from_1746_relaxed": len(good),
        "omega_phi_optimum": {"omega": th[0], "phi": th[1], "objective": fbest},
        "hessian": {
            "condition_number": cond,
            "eigenvalues": [float(v) for v in eig],
            "regularization": reg,
            "ci95_halfwidth_hessian": {"omega": ci_h[0], "phi": ci_h[1]},
            "max_abs_offdiag_corr": max_corr,
            "corr": corr.tolist(),
        },
        "bootstrap_ci95": ci_boot,
        "bootstrap_phi_circular_std": phi_std,
        "beta_sensitivity_sweep": {
            "beta_grid_min": beta_lo,
            "beta_grid_max": beta_hi,
            "n_grid": len(beta_grid),
            "omega_boundary_rate": boundary_rate,
            "rows": sweep,
        },
        "pass_flags": {
            "beta_quality": bool(pass_beta_quality),
            "omega_width": bool(pass_omega_width),
            "phi_width": bool(pass_phi_width),
            "conditioning_or_corr": bool(pass_cond or pass_corr),
            "nonboundary_omega": bool(pass_nonboundary),
            "beta_sweep_nonboundary": bool(pass_sweep),
        },
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1750: SEPARATED BETA THEN OMEGA/PHI",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- beta_fixed (from 1749): {beta_fix:.6f}",
        f"- Good rows (1746 relaxed): {len(good)}",
        f"- Optimum: omega={th[0]:.6f}, phi={th[1]:.6f}",
        f"- Hessian cond: {cond:.3e}, max|corr_offdiag|={max_corr:.4f}",
        f"- Bootstrap omega width: {width_omega:.4f}, phi circular std: {phi_std:.4f}",
        f"- Beta-sweep omega boundary rate: {boundary_rate:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- beta_quality: {pass_beta_quality}",
        f"- omega_width: {pass_omega_width}",
        f"- phi_width: {pass_phi_width}",
        f"- conditioning_or_corr: {pass_cond or pass_corr}",
        f"- nonboundary_omega: {pass_nonboundary}",
        f"- beta_sweep_nonboundary: {pass_sweep}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1750] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1750] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
