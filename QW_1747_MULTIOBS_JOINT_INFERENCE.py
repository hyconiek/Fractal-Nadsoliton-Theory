#!/usr/bin/env python3
"""
QW-1747: Joint inference from dynamic multi-observables (from QW-1746).

Combines independent estimators:
- omega_phase_hat
- omega_zero_hat
- beta_env_hat
- phi_phase_hat
into a global posterior-like estimate with bootstrap uncertainty.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1747_multiobs_joint_inference.json"
OUT_MD = ROOT / "RAPORT_QW1747_MULTIOBS_JOINT_INFERENCE.md"


def load_1746() -> Dict:
    p = ROOT / "report_qw1746_dynamic_observables_derivation.json"
    return json.loads(p.read_text(encoding="utf-8"))


def wrap_phi(x: float) -> float:
    return float((x + math.pi) % (2 * math.pi) - math.pi)


def circ_mean(x: np.ndarray) -> float:
    return float(np.arctan2(np.mean(np.sin(x)), np.mean(np.cos(x))))


def circ_std(x: np.ndarray) -> float:
    r = np.sqrt(np.mean(np.sin(x)) ** 2 + np.mean(np.cos(x)) ** 2)
    if r <= 1e-12:
        return float("inf")
    return float(np.sqrt(-2.0 * np.log(r)))


def objective(theta: Tuple[float, float, float], data: Dict[str, np.ndarray]) -> float:
    omega, phi, beta = theta
    if not (0.10 <= omega <= 1.6 and -math.pi <= phi <= math.pi and 0.001 <= beta <= 0.25):
        return float("inf")

    om_p = data["omega_phase"]
    om_z = data["omega_zero"]
    be = data["beta"]
    ph = data["phi"]
    wp = data["w_phase"]
    wz = data["w_zero"]
    wb = data["w_beta"]
    wh = data["w_phi"]

    res = 0.0
    res += float(np.sum(wp * (om_p - omega) ** 2))
    res += float(np.sum(wz * (om_z - omega) ** 2))
    res += float(np.sum(wb * (be - beta) ** 2))
    dphi = np.array([wrap_phi(float(x - phi)) for x in ph], dtype=float)
    res += float(np.sum(wh * (dphi ** 2)))
    return res


def refine(start: Tuple[float, float, float], data: Dict[str, np.ndarray]) -> Tuple[Tuple[float, float, float], float]:
    cur = (float(start[0]), wrap_phi(float(start[1])), float(start[2]))
    fcur = objective(cur, data)
    steps = [(0.20, 0.70, 0.04), (0.08, 0.25, 0.015), (0.03, 0.10, 0.006), (0.012, 0.04, 0.0025)]
    for so, sp, sb in steps:
        improved = True
        while improved:
            improved = False
            # omega
            best = (cur, fcur)
            for om in np.linspace(max(0.10, cur[0] - so), min(1.6, cur[0] + so), 9):
                th = (float(om), cur[1], cur[2])
                f = objective(th, data)
                if f < best[1]:
                    best = (th, f)
            if best[1] < fcur:
                cur, fcur = best
                improved = True

            # phi
            best = (cur, fcur)
            for ph in np.linspace(cur[1] - sp, cur[1] + sp, 9):
                th = (cur[0], wrap_phi(float(ph)), cur[2])
                f = objective(th, data)
                if f < best[1]:
                    best = (th, f)
            if best[1] < fcur:
                cur, fcur = best
                improved = True

            # beta
            best = (cur, fcur)
            for be in np.linspace(max(0.001, cur[2] - sb), min(0.25, cur[2] + sb), 9):
                th = (cur[0], cur[1], float(be))
                f = objective(th, data)
                if f < best[1]:
                    best = (th, f)
            if best[1] < fcur:
                cur, fcur = best
                improved = True
    return cur, fcur


def hessian(theta: Tuple[float, float, float], data: Dict[str, np.ndarray]) -> np.ndarray:
    x = np.array(theta, dtype=float)
    h = np.array([0.006, 0.02, 0.003], dtype=float)
    f0 = objective(tuple(x), data)
    H = np.zeros((3, 3), dtype=float)
    for i in range(3):
        xp = x.copy(); xm = x.copy()
        xp[i] += h[i]; xm[i] -= h[i]
        if i == 1:
            xp[i] = wrap_phi(float(xp[i])); xm[i] = wrap_phi(float(xm[i]))
        fp = objective(tuple(xp), data); fm = objective(tuple(xm), data)
        H[i, i] = (fp - 2 * f0 + fm) / (h[i] ** 2)
    for i in range(3):
        for j in range(i + 1, 3):
            xpp = x.copy(); xpm = x.copy(); xmp = x.copy(); xmm = x.copy()
            xpp[i] += h[i]; xpp[j] += h[j]
            xpm[i] += h[i]; xpm[j] -= h[j]
            xmp[i] -= h[i]; xmp[j] += h[j]
            xmm[i] -= h[i]; xmm[j] -= h[j]
            if i == 1:
                xpp[i] = wrap_phi(float(xpp[i])); xpm[i] = wrap_phi(float(xpm[i]))
                xmp[i] = wrap_phi(float(xmp[i])); xmm[i] = wrap_phi(float(xmm[i]))
            if j == 1:
                xpp[j] = wrap_phi(float(xpp[j])); xpm[j] = wrap_phi(float(xpm[j]))
                xmp[j] = wrap_phi(float(xmp[j])); xmm[j] = wrap_phi(float(xmm[j]))
            fpp = objective(tuple(xpp), data); fpm = objective(tuple(xpm), data)
            fmp = objective(tuple(xmp), data); fmm = objective(tuple(xmm), data)
            v = (fpp - fpm - fmp + fmm) / (4 * h[i] * h[j])
            H[i, j] = v; H[j, i] = v
    return H


def bootstrap(data: Dict[str, np.ndarray], th0: Tuple[float, float, float], n_boot: int = 200, seed: int = 1747) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n = len(data["omega_phase"])
    out = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        db = {k: v[idx] for k, v in data.items()}
        th, _ = refine(th0, db)
        out.append(th)
    return np.array(out, dtype=float)


def main() -> None:
    d1746 = load_1746()
    rows = d1746.get("rows", [])
    chosen_label = d1746.get("subset_sizes", {}).get("chosen_label", "strict")
    if chosen_label == "strict":
        key = "is_good_strict"
    else:
        key = "is_good_relaxed"

    good = [r for r in rows if bool(r.get(key, False))]

    if len(good) < 8:
        raise RuntimeError("Too few good rows from QW-1746 for joint inference.")

    om_p = np.array([r["omega_phase_hat"] for r in good], dtype=float)
    om_z = np.array([r["omega_zero_hat"] for r in good], dtype=float)
    be = np.array([r["beta_env_hat"] for r in good], dtype=float)
    ph = np.array([r["phi_phase_hat"] for r in good], dtype=float)

    # Weights from per-run fit quality and residuals.
    wp = 1.0 / np.clip(np.array([r["phase_fit_rmse"] for r in good], dtype=float) ** 2, 1e-5, None)
    wz = 1.0 / np.clip((0.6 + np.array([r["phase_fit_rmse"] for r in good], dtype=float)) ** 2, 1e-5, None)
    wb = 1.0 / np.clip(np.array([r["beta_fit_rmse"] for r in good], dtype=float) ** 2, 1e-5, None)
    wh = wp.copy()

    data = {"omega_phase": om_p, "omega_zero": om_z, "beta": be, "phi": ph, "w_phase": wp, "w_zero": wz, "w_beta": wb, "w_phi": wh}

    starts = [
        (float(np.median(om_p)), float(circ_mean(ph)), float(np.median(be))),
        (float(np.median(om_z)), float(circ_mean(ph) + 0.4), float(np.median(be) * 0.8)),
        (0.75, -0.2, 0.03),
        (0.45, 0.9, 0.10),
    ]
    sols = sorted([refine(st, data) for st in starts], key=lambda x: x[1])
    th, fbest = sols[0]

    # Hessian/Fisher local diagnostics
    H = hessian(th, data)
    H = 0.5 * (H + H.T)
    H = np.nan_to_num(H, nan=1e6, posinf=1e6, neginf=-1e6)
    eig = np.linalg.eigvalsh(H)
    min_e = float(np.min(eig))
    reg = 0.0
    if min_e <= 1e-8:
        reg = (1e-8 - min_e) + 1e-6
    Hr = H + reg * np.eye(3)
    cond = float(np.linalg.cond(Hr))
    cov = np.linalg.inv(Hr)
    std = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    ci_h = [float(1.96 * std[0]), float(1.96 * std[1]), float(1.96 * std[2])]
    corr = cov / np.sqrt(np.outer(np.diag(cov), np.diag(cov)))
    corr = np.clip(corr, -1.0, 1.0)
    max_corr = float(np.max(np.abs(corr - np.eye(3))))

    # Bootstrap
    boot = bootstrap(data, th, n_boot=200, seed=1747)
    bo, bp, bb = boot[:, 0], boot[:, 1], boot[:, 2]
    ci_boot = {
        "omega": [float(np.quantile(bo, 0.025)), float(np.quantile(bo, 0.975))],
        "phi": [float(np.quantile(bp, 0.025)), float(np.quantile(bp, 0.975))],
        "beta": [float(np.quantile(bb, 0.025)), float(np.quantile(bb, 0.975))],
    }

    pass_cond = cond <= 1e5
    pass_corr = max_corr <= 0.90
    pass_ci = (ci_boot["omega"][1] - ci_boot["omega"][0]) <= 0.24 and (ci_boot["beta"][1] - ci_boot["beta"][0]) <= 0.09
    pass_nonboundary = 0.15 < th[0] < 1.5 and 0.001 < th[2] < 0.25

    if pass_cond and pass_corr and pass_ci and pass_nonboundary:
        verdict = "MULTIOBS_JOINT_IDENTIFIABILITY_STRONG"
    elif pass_ci and (pass_cond or pass_corr):
        verdict = "MULTIOBS_JOINT_IDENTIFIABILITY_MODERATE"
    else:
        verdict = "MULTIOBS_JOINT_IDENTIFIABILITY_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_good_rows": len(good),
        "subset_label_from_1746": chosen_label,
        "optimum": {"omega": th[0], "phi": th[1], "beta": th[2], "objective": fbest},
        "hessian": {
            "condition_number": cond,
            "eigenvalues": [float(v) for v in eig],
            "regularization": reg,
            "ci95_halfwidth_hessian": {"omega": ci_h[0], "phi": ci_h[1], "beta": ci_h[2]},
            "max_abs_offdiag_corr": max_corr,
            "corr": corr.tolist(),
        },
        "bootstrap_ci95": ci_boot,
        "pass_flags": {
            "conditioning": pass_cond,
            "correlation": pass_corr,
            "bootstrap_width": pass_ci,
            "nonboundary": pass_nonboundary,
        },
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1747: MULTI-OBS JOINT INFERENCE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Good rows: {len(good)}",
        f"- Optimum: omega={th[0]:.6f}, phi={th[1]:.6f}, beta={th[2]:.6f}",
        f"- Hessian cond: {cond:.3e}, max|corr_offdiag|={max_corr:.4f}",
        f"- Bootstrap widths: domega={ci_boot['omega'][1]-ci_boot['omega'][0]:.4f}, dbeta={ci_boot['beta'][1]-ci_boot['beta'][0]:.4f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- conditioning: {pass_cond}",
        f"- correlation: {pass_corr}",
        f"- bootstrap_width: {pass_ci}",
        f"- nonboundary: {pass_nonboundary}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1747] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1747] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
