#!/usr/bin/env python3
"""
QW-1744: Identifiability audit on oscillatory cohort (from QW-1743).
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1744_oscillatory_cohort_identifiability.json"
OUT_MD = ROOT / "RAPORT_QW1744_OSCILLATORY_COHORT_IDENTIFIABILITY.md"


def load_1743() -> Dict:
    p = ROOT / "report_qw1743_oscillatory_cohort_derivation.json"
    return json.loads(p.read_text(encoding="utf-8"))


def build_dataset(d1743: Dict) -> Tuple[np.ndarray, np.ndarray]:
    rows = d1743.get("selected_rows_head_120", [])
    y = []
    for r in rows:
        prof = r["profile"]
        y.append([float(prof[str(i)]) for i in range(1, 13)])
    y = np.array(y, dtype=float)
    w = 1.0 / np.clip(np.var(y, axis=1), 1e-5, None)
    return y, w


def robust_loss(r: np.ndarray, delta: float = 0.08) -> np.ndarray:
    a = np.abs(r)
    return np.where(a <= delta, 0.5 * a * a, delta * (a - 0.5 * delta))


def basis(d: np.ndarray, om: float, ph: float, be: float) -> np.ndarray:
    return np.cos(om * d + ph) / (1.0 + be * d)


def obj(theta: Tuple[float, float, float], y: np.ndarray, w: np.ndarray) -> float:
    om, ph, be = theta
    if not (0.20 <= om <= 1.5 and -math.pi <= ph <= math.pi and 0.001 <= be <= 0.20):
        return float("inf")
    d = np.arange(1, y.shape[1] + 1, dtype=float)
    b = basis(d, om, ph, be)
    bb = float(np.dot(b, b))
    if bb <= 1e-12:
        return float("inf")
    a = (y @ b) / bb
    p = a[:, None] * b[None, :]
    r = y - p
    return float(np.sum(w[:, None] * robust_loss(r)))


def wrap_phi(x: float) -> float:
    return float((x + math.pi) % (2 * math.pi) - math.pi)


def opt_given_fixed(name: str, val: float, start: Tuple[float, float, float], y: np.ndarray, w: np.ndarray) -> Tuple[Tuple[float, float, float], float]:
    idx = {"omega": 0, "phi": 1, "beta": 2}[name]
    cur = list(start)
    cur[idx] = val
    cur = (float(cur[0]), float(cur[1]), float(cur[2]))
    fcur = obj(cur, y, w)

    steps = [(0.10, 0.30, 0.015), (0.04, 0.12, 0.006), (0.015, 0.05, 0.0025)]
    for so, sp, sb in steps:
        improved = True
        while improved:
            improved = False
            for j in [0, 1, 2]:
                if j == idx:
                    continue
                best = (cur, fcur)
                if j == 0:
                    for x in np.linspace(max(0.20, cur[0] - so), min(1.5, cur[0] + so), 9):
                        th = [cur[0], cur[1], cur[2]]
                        th[0] = float(x)
                        th[idx] = val
                        f = obj((th[0], th[1], th[2]), y, w)
                        if f < best[1]:
                            best = ((th[0], th[1], th[2]), f)
                elif j == 1:
                    for x in np.linspace(cur[1] - sp, cur[1] + sp, 9):
                        th = [cur[0], cur[1], cur[2]]
                        th[1] = wrap_phi(float(x))
                        th[idx] = val
                        f = obj((th[0], th[1], th[2]), y, w)
                        if f < best[1]:
                            best = ((th[0], th[1], th[2]), f)
                else:
                    for x in np.linspace(max(0.001, cur[2] - sb), min(0.20, cur[2] + sb), 9):
                        th = [cur[0], cur[1], cur[2]]
                        th[2] = float(x)
                        th[idx] = val
                        f = obj((th[0], th[1], th[2]), y, w)
                        if f < best[1]:
                            best = ((th[0], th[1], th[2]), f)
                if best[1] < fcur:
                    cur, fcur = best
                    improved = True
    return cur, fcur


def interval(x: np.ndarray, dlt: np.ndarray, thr: float) -> List[float]:
    m = dlt <= thr
    if not np.any(m):
        return [float("nan"), float("nan")]
    return [float(np.min(x[m])), float(np.max(x[m]))]


def hessian(theta: Tuple[float, float, float], y: np.ndarray, w: np.ndarray) -> np.ndarray:
    x = np.array(theta, dtype=float)
    h = np.array([0.006, 0.02, 0.003], dtype=float)
    f0 = obj(tuple(x), y, w)
    H = np.zeros((3, 3), dtype=float)
    for i in range(3):
        xp = x.copy(); xm = x.copy()
        xp[i] += h[i]; xm[i] -= h[i]
        if i == 1:
            xp[i] = wrap_phi(float(xp[i])); xm[i] = wrap_phi(float(xm[i]))
        fp = obj(tuple(xp), y, w); fm = obj(tuple(xm), y, w)
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
            fpp = obj(tuple(xpp), y, w); fpm = obj(tuple(xpm), y, w)
            fmp = obj(tuple(xmp), y, w); fmm = obj(tuple(xmm), y, w)
            v = (fpp - fpm - fmp + fmm) / (4 * h[i] * h[j])
            H[i, j] = v; H[j, i] = v
    return H


def main() -> None:
    d1743 = load_1743()
    y, w = build_dataset(d1743)
    th0 = (
        float(d1743["optimum"]["omega"]),
        float(d1743["optimum"]["phi"]),
        float(d1743["optimum"]["beta"]),
    )
    f0 = obj(th0, y, w)

    og = np.linspace(max(0.20, th0[0] - 0.25), min(1.5, th0[0] + 0.25), 121)
    pg = np.linspace(th0[1] - 1.0, th0[1] + 1.0, 121)
    bg = np.linspace(max(0.001, th0[2] - 0.07), min(0.20, th0[2] + 0.07), 121)

    pwo, pph, pbe = [], [], []
    for v in og:
        _, f = opt_given_fixed("omega", float(v), th0, y, w)
        pwo.append((float(v), float(f)))
    for v in pg:
        vv = wrap_phi(float(v))
        _, f = opt_given_fixed("phi", vv, th0, y, w)
        pph.append((vv, float(f)))
    for v in bg:
        _, f = opt_given_fixed("beta", float(v), th0, y, w)
        pbe.append((float(v), float(f)))

    d_om = np.array([x[1] - f0 for x in pwo], dtype=float)
    d_ph = np.array([x[1] - f0 for x in pph], dtype=float)
    d_be = np.array([x[1] - f0 for x in pbe], dtype=float)

    ci_om = interval(np.array([x[0] for x in pwo]), d_om, 3.84)
    ci_ph = interval(np.array([x[0] for x in pph]), d_ph, 3.84)
    ci_be = interval(np.array([x[0] for x in pbe]), d_be, 3.84)

    w_om = ci_om[1] - ci_om[0] if np.isfinite(ci_om[0]) and np.isfinite(ci_om[1]) else float("inf")
    w_ph = ci_ph[1] - ci_ph[0] if np.isfinite(ci_ph[0]) and np.isfinite(ci_ph[1]) else float("inf")
    w_be = ci_be[1] - ci_be[0] if np.isfinite(ci_be[0]) and np.isfinite(ci_be[1]) else float("inf")

    H = hessian(th0, y, w)
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
    corr = cov / np.sqrt(np.outer(np.diag(cov), np.diag(cov)))
    corr = np.clip(corr, -1.0, 1.0)
    mc = float(np.max(np.abs(corr - np.eye(3))))

    pass_prof = w_om <= 0.22 and w_ph <= 0.95 and w_be <= 0.08
    pass_cond = cond <= 1e5
    pass_corr = mc <= 0.90

    if pass_prof and pass_cond and pass_corr:
        verdict = "OSCILLATORY_IDENTIFIABILITY_STRONG"
    elif pass_prof and (pass_cond or pass_corr):
        verdict = "OSCILLATORY_IDENTIFIABILITY_MODERATE"
    else:
        verdict = "OSCILLATORY_IDENTIFIABILITY_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_optimum": {"omega": th0[0], "phi": th0[1], "beta": th0[2], "objective": f0},
        "profile_ci95_like": {"omega": ci_om, "phi": ci_ph, "beta": ci_be},
        "profile_width": {"omega": w_om, "phi": w_ph, "beta": w_be},
        "hessian": {
            "condition_number": cond,
            "eigenvalues": [float(v) for v in eig],
            "regularization": reg,
            "max_abs_offdiag_corr": mc,
            "corr": corr.tolist(),
        },
        "pass_flags": {"profile_width": pass_prof, "conditioning": pass_cond, "correlation": pass_corr},
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1744: OSCILLATORY COHORT IDENTIFIABILITY",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Optimum: omega={th0[0]:.6f}, phi={th0[1]:.6f}, beta={th0[2]:.6f}",
        f"- Widths: domega={w_om:.4f}, dphi={w_ph:.4f}, dbeta={w_be:.4f}",
        f"- Hessian cond: {cond:.3e}, max|corr_offdiag|={mc:.4f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- profile_width: {pass_prof}",
        f"- conditioning: {pass_cond}",
        f"- correlation: {pass_corr}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1744] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1744] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
