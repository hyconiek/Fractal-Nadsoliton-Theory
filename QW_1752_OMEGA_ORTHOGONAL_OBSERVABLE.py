#!/usr/bin/env python3
"""
QW-1752: Omega-orthogonal observable from spatial phase structure.

Core idea:
- Estimate omega from phase-vs-distance (complex lock-in response + analytic signal),
  which is weakly sensitive to amplitude envelope (beta-like effects).
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1752_omega_orthogonal_observable.json"
OUT_MD = ROOT / "RAPORT_QW1752_OMEGA_ORTHOGONAL_OBSERVABLE.md"


@dataclass
class Cfg:
    n: int
    seed: int
    xi: float
    p_short: float
    eta: float
    tau: float
    rho: float
    c2: float
    damp: float
    omega_drive: float
    drive_amp: float
    n_steps: int = 1600
    burn: int = 700


def ring_dist(i: int, j: int, n: int) -> int:
    d = abs(i - j)
    return min(d, n - d)


def smooth_periodic(x: np.ndarray, iters: int = 7) -> np.ndarray:
    y = x.copy()
    for _ in range(iters):
        y = 0.25 * np.roll(y, -1) + 0.5 * y + 0.25 * np.roll(y, 1)
    return y


def build_w(cfg: Cfg, rng: np.random.Generator) -> np.ndarray:
    n = cfg.n
    th = math.pi * np.tanh(smooth_periodic(rng.normal(size=n), iters=8) / 1.25)
    q = rng.integers(-2, 3, size=n).astype(float)

    w = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            d = ring_dist(i, j, n)
            amp = math.exp(-d / cfg.xi)
            if rng.random() < cfg.p_short / (d ** cfg.eta):
                amp += 0.9 * (d ** -0.95) * (1.0 + 0.15 * rng.normal())

            dt = abs(th[i] - th[j]) / math.pi
            dq = abs(q[i] - q[j])
            tors = dq + 0.35 * dt

            sgn = 1.0 if (math.cos(th[i] - th[j]) + 0.22 * (q[i] - q[j])) >= 0 else -1.0
            sym = sgn * amp * math.exp(-cfg.tau * tors)
            anti = cfg.rho * math.sin(th[i] - th[j]) / (d ** 1.1)
            w[i, j] = sym + anti
            w[j, i] = sym - anti

    smax = float(np.linalg.norm(w, 2))
    if smax > 1e-12:
        w /= smax
    return w


def lockin_response(cfg: Cfg, w: np.ndarray, src: int) -> np.ndarray:
    n = cfg.n
    lap = w - np.eye(n)
    x_prev = np.zeros(n, dtype=float)
    x_curr = np.zeros(n, dtype=float)
    x_curr[src] = 1.0

    re_acc = np.zeros(n, dtype=float)
    im_acc = np.zeros(n, dtype=float)
    count = 0

    for t in range(cfg.n_steps):
        drive = cfg.drive_amp * math.sin(cfg.omega_drive * t)
        forcing = np.zeros(n, dtype=float)
        forcing[src] = drive
        x_next = 2.0 * x_curr - x_prev + cfg.c2 * (lap @ x_curr) - cfg.damp * (x_curr - x_prev) + forcing
        x_next = np.tanh(1.3 * x_next)
        x_prev, x_curr = x_curr, x_next

        if t >= cfg.burn:
            c = math.cos(cfg.omega_drive * t)
            s = math.sin(cfg.omega_drive * t)
            re_acc += x_curr * c
            im_acc -= x_curr * s
            count += 1

    z = (re_acc + 1j * im_acc) / max(count, 1)
    return z


def distance_response(z: np.ndarray, src: int, dmax: int = 12) -> np.ndarray:
    n = len(z)
    out = []
    for d in range(1, min(dmax, n // 2) + 1):
        j1 = (src + d) % n
        j2 = (src - d) % n
        out.append(0.5 * (z[j1] + z[j2]))
    return np.array(out, dtype=complex)


def fit_phase_slope_complex(resp: np.ndarray) -> Tuple[float, float, float]:
    d = np.arange(1, len(resp) + 1, dtype=float)
    amp = np.abs(resp)
    ph = np.unwrap(np.angle(resp))
    m = amp >= np.quantile(amp, 0.35)
    if np.sum(m) < 4:
        m = amp >= np.quantile(amp, 0.2)
    x = d[m]
    y = ph[m]
    A = np.column_stack([x, np.ones_like(x)])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    slope, intercept = float(coef[0]), float(coef[1])
    rmse = float(np.sqrt(np.mean((y - (slope * x + intercept)) ** 2)))
    return slope, intercept, rmse


def analytic_signal_real(y: np.ndarray) -> np.ndarray:
    n = len(y)
    Y = np.fft.fft(y)
    H = np.zeros(n)
    if n % 2 == 0:
        H[0] = 1
        H[n // 2] = 1
        H[1 : n // 2] = 2
    else:
        H[0] = 1
        H[1 : (n + 1) // 2] = 2
    z = np.fft.ifft(Y * H)
    return z


def fit_phase_slope_hilbert(y_real: np.ndarray) -> Tuple[float, float, float]:
    d = np.arange(1, len(y_real) + 1, dtype=float)
    z = analytic_signal_real(y_real)
    amp = np.abs(z)
    ph = np.unwrap(np.angle(z))
    m = amp >= np.quantile(amp, 0.35)
    if np.sum(m) < 4:
        m = amp >= np.quantile(amp, 0.2)
    x = d[m]
    y = ph[m]
    A = np.column_stack([x, np.ones_like(x)])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    slope, intercept = float(coef[0]), float(coef[1])
    rmse = float(np.sqrt(np.mean((y - (slope * x + intercept)) ** 2)))
    return slope, intercept, rmse


def omega_from_phase_slopes(s1: float, s2: float) -> float:
    vals = [abs(s1), abs(s2)]
    vals = [v for v in vals if np.isfinite(v)]
    if not vals:
        return float("nan")
    return float(np.median(np.array(vals, dtype=float)))


def beta_envelope_perturb(resp: np.ndarray, beta_val: float) -> np.ndarray:
    d = np.arange(1, len(resp) + 1, dtype=float)
    env = 1.0 / (1.0 + beta_val * d)
    return resp * env


def main() -> None:
    rows: List[Dict[str, object]] = []
    for n in [96, 120]:
        for k in range(18):
            seed = 175200 + 100 * n + k
            cfg = Cfg(
                n=n,
                seed=seed,
                xi=1.55 + 0.22 * ((k % 3) - 1),
                p_short=0.10 + 0.05 * (k % 2),
                eta=1.45 + 0.18 * ((k + 1) % 3),
                tau=0.45 + 0.16 * (k % 3),
                rho=0.18 + 0.10 * (k % 2),
                c2=0.30 + 0.08 * (k % 3),
                damp=0.03 + 0.02 * (k % 2),
                omega_drive=0.22 + 0.06 * (k % 4),
                drive_amp=0.9 + 0.3 * (k % 2),
            )
            rng = np.random.default_rng(seed)
            w = build_w(cfg, rng)
            src = seed % n
            z = lockin_response(cfg, w, src=src)
            resp = distance_response(z, src=src, dmax=12)

            slope_c, phi_c, rmse_c = fit_phase_slope_complex(resp)
            y_real = np.real(resp)
            slope_h, phi_h, rmse_h = fit_phase_slope_hilbert(y_real)
            omega_hat = omega_from_phase_slopes(slope_c, slope_h)
            phi_hat = float(0.5 * (phi_c + phi_h))

            # beta perturbation stability test: omega should remain stable under envelope scaling
            beta_grid = np.linspace(0.0, 0.20, 9)
            omega_pert = []
            for b in beta_grid:
                rp = beta_envelope_perturb(resp, float(b))
                sc2, _pc2, _rc2 = fit_phase_slope_complex(rp)
                sh2, _ph2, _rh2 = fit_phase_slope_hilbert(np.real(rp))
                omega_pert.append(omega_from_phase_slopes(sc2, sh2))
            omega_pert = np.array(omega_pert, dtype=float)
            delta_beta = float(np.nanmax(np.abs(omega_pert - omega_hat)))

            amp = np.abs(resp)
            snr_like = float(np.max(amp) / max(np.median(np.abs(amp - np.median(amp))), 1e-6))

            rows.append(
                {
                    "n": n,
                    "seed": seed,
                    "omega_drive": cfg.omega_drive,
                    "omega_hat": float(omega_hat),
                    "phi_hat": float(phi_hat),
                    "slope_complex": float(slope_c),
                    "slope_hilbert": float(slope_h),
                    "phase_rmse_complex": float(rmse_c),
                    "phase_rmse_hilbert": float(rmse_h),
                    "delta_omega_beta_perturb": delta_beta,
                    "snr_like": snr_like,
                    "response_real": {str(i + 1): float(np.real(resp[i])) for i in range(len(resp))},
                    "response_abs": {str(i + 1): float(np.abs(resp[i])) for i in range(len(resp))},
                }
            )

    good = [
        r for r in rows
        if np.isfinite(r["omega_hat"])
        and r["snr_like"] >= 4.0
        and r["phase_rmse_complex"] <= 0.45
        and r["phase_rmse_hilbert"] <= 0.55
    ]

    om = np.array([r["omega_hat"] for r in good], dtype=float) if good else np.array([], dtype=float)
    ph = np.array([r["phi_hat"] for r in good], dtype=float) if good else np.array([], dtype=float)
    dbe = np.array([r["delta_omega_beta_perturb"] for r in good], dtype=float) if good else np.array([], dtype=float)
    rmse_c = np.array([r["phase_rmse_complex"] for r in good], dtype=float) if good else np.array([], dtype=float)
    rmse_h = np.array([r["phase_rmse_hilbert"] for r in good], dtype=float) if good else np.array([], dtype=float)

    if len(good) >= 8:
        summary = {
            "n_good": len(good),
            "omega_median": float(np.median(om)),
            "omega_iqr": float(np.quantile(om, 0.75) - np.quantile(om, 0.25)),
            "omega_ci95_empirical": [float(np.quantile(om, 0.025)), float(np.quantile(om, 0.975))],
            "phi_circular_mean": float(math.atan2(np.mean(np.sin(ph)), np.mean(np.cos(ph)))),
            "phi_circular_std": float(np.sqrt(-2.0 * np.log(max(np.sqrt(np.mean(np.sin(ph)) ** 2 + np.mean(np.cos(ph)) ** 2), 1e-12)))),
            "median_phase_rmse_complex": float(np.median(rmse_c)),
            "median_phase_rmse_hilbert": float(np.median(rmse_h)),
            "delta_omega_beta_perturb_q90": float(np.quantile(dbe, 0.9)),
        }
    else:
        summary = {
            "n_good": len(good),
            "omega_median": None,
            "omega_iqr": None,
            "omega_ci95_empirical": [None, None],
            "phi_circular_mean": None,
            "phi_circular_std": None,
            "median_phase_rmse_complex": None,
            "median_phase_rmse_hilbert": None,
            "delta_omega_beta_perturb_q90": None,
        }

    pass_good = len(good) >= 12
    pass_spread = summary["omega_iqr"] is not None and summary["omega_iqr"] <= 0.35
    pass_beta_orth = summary["delta_omega_beta_perturb_q90"] is not None and summary["delta_omega_beta_perturb_q90"] <= 0.08
    pass_fit = summary["median_phase_rmse_complex"] is not None and summary["median_phase_rmse_complex"] <= 0.35

    if pass_good and pass_spread and pass_beta_orth and pass_fit:
        verdict = "OMEGA_ORTHOGONAL_OBSERVABLE_SUPPORTED"
    elif pass_fit and (pass_spread or pass_beta_orth):
        verdict = "OMEGA_ORTHOGONAL_OBSERVABLE_PARTIAL"
    else:
        verdict = "OMEGA_ORTHOGONAL_OBSERVABLE_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "n_total_runs": len(rows),
            "good_filter": {
                "snr_like_min": 4.0,
                "phase_rmse_complex_max": 0.45,
                "phase_rmse_hilbert_max": 0.55,
            },
            "beta_perturb_grid": [0.0, 0.2, 9],
        },
        "summary": summary,
        "pass_flags": {
            "n_good_enough": bool(pass_good),
            "spread_control": bool(pass_spread),
            "beta_orthogonality": bool(pass_beta_orth),
            "fit_quality": bool(pass_fit),
        },
        "verdict": verdict,
        "rows": rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1752: OMEGA ORTHOGONAL OBSERVABLE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Total runs: {len(rows)}",
        f"- Good runs: {summary['n_good']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- n_good_enough: {pass_good}",
        f"- spread_control: {pass_spread}",
        f"- beta_orthogonality: {pass_beta_orth}",
        f"- fit_quality: {pass_fit}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1752] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1752] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
