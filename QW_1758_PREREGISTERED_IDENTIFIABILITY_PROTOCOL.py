#!/usr/bin/env python3
"""
QW-1758: Preregistered identifiability protocol (beta, omega, phi).

Design principles:
1) Fixed thresholds declared before execution.
2) Multi-source measurements per topology sample.
3) Extended distance window for both envelope and phase observables.
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
OUT_JSON = ROOT / "report_qw1758_preregistered_identifiability_protocol.json"
OUT_MD = ROOT / "RAPORT_QW1758_PREREGISTERED_IDENTIFIABILITY_PROTOCOL.md"


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
    n_steps_imp: int = 900
    tail_start: int = 320
    tail_end: int = 860
    n_steps_lock: int = 1400
    burn_lock: int = 620
    dmax: int = 16


PREREG = {
    "good_row_filter": {
        "snr_like_min": 5.0,
        "beta_rmse_max": 0.10,
        "phase_rmse_complex_max": 0.40,
        "phase_rmse_hilbert_max": 0.50,
    },
    "pass_thresholds": {
        "n_good_min": 18,
        "beta_evidence_median_dbic_min": 0.5,
        "beta_evidence_positive_rate_min": 0.55,
        "beta_nonboundary_rate_min": 0.65,
        "omega_nonboundary_rate_min": 0.75,
        "beta_iqr_max": 0.015,
        "omega_iqr_max": 0.30,
        "rho_beta_omega_abs_max": 0.35,
        "delta_omega_beta_pert_q90_max": 0.06,
        "median_phase_rmse_complex_max": 0.35,
        "median_phase_rmse_hilbert_max": 0.45,
    },
}


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
                amp += 0.8 * (d ** -0.95) * (1.0 + 0.14 * rng.normal())

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


def impulse_response_timeseries(cfg: Cfg, w: np.ndarray, src: int) -> np.ndarray:
    n = cfg.n
    lap = w - np.eye(n)
    x_prev = np.zeros(n, dtype=float)
    x_curr = np.zeros(n, dtype=float)
    x_curr[src] = 1.0

    hist = []
    for _ in range(cfg.n_steps_imp):
        x_next = 2.0 * x_curr - x_prev + cfg.c2 * (lap @ x_curr) - cfg.damp * (x_curr - x_prev)
        x_next = np.tanh(1.35 * x_next)
        x_prev, x_curr = x_curr, x_next
        hist.append(x_curr.copy())
    return np.array(hist, dtype=float)


def distance_tail_observable(ts: np.ndarray, src: int, tail_start: int, tail_end: int, dmax: int) -> np.ndarray:
    t1 = max(0, tail_start)
    t2 = min(ts.shape[0], tail_end)
    tail = np.abs(ts[t1:t2])
    n = ts.shape[1]

    vals = []
    for d in range(1, min(dmax, n // 2) + 1):
        j1 = (src + d) % n
        j2 = (src - d) % n
        arr = 0.5 * (tail[:, j1] + tail[:, j2])
        vals.append(float(np.median(arr)))
    y = np.array(vals, dtype=float)
    if y[0] > 1e-12:
        y /= y[0]
    return y


def fit_beta_envelope(y: np.ndarray) -> Tuple[float, float, float]:
    d = np.arange(1, len(y) + 1, dtype=float)
    # Preregistered log-grid includes near-zero but avoids hard 0.
    beta_grid = np.logspace(np.log10(1e-5), np.log10(0.25), 3200)

    best_beta = 1e-3
    best_rmse = float("inf")
    best_sse = float("inf")
    for b in beta_grid:
        env = 1.0 / (1.0 + b * d)
        den = float(np.dot(env, env))
        if den <= 1e-12:
            continue
        a = float(np.dot(y, env) / den)
        pred = a * env
        sse = float(np.sum((y - pred) ** 2))
        rmse = float(np.sqrt(np.mean((y - pred) ** 2)))
        if rmse < best_rmse:
            best_rmse = rmse
            best_beta = float(b)
            best_sse = sse

    # Null model (constant envelope).
    a0 = float(np.mean(y))
    sse0 = float(np.sum((y - a0) ** 2))
    n = len(y)
    sse0 = max(sse0, 1e-15)
    sse1 = max(best_sse, 1e-15)
    bic0 = float(n * np.log(sse0 / n) + 1.0 * np.log(n))
    bic1 = float(n * np.log(sse1 / n) + 2.0 * np.log(n))
    delta_bic_null_minus_beta = bic0 - bic1

    return best_beta, best_rmse, float(delta_bic_null_minus_beta)


def lockin_response(cfg: Cfg, w: np.ndarray, src: int) -> np.ndarray:
    n = cfg.n
    lap = w - np.eye(n)
    x_prev = np.zeros(n, dtype=float)
    x_curr = np.zeros(n, dtype=float)
    x_curr[src] = 1.0

    re_acc = np.zeros(n, dtype=float)
    im_acc = np.zeros(n, dtype=float)
    cnt = 0

    for t in range(cfg.n_steps_lock):
        drive = cfg.drive_amp * math.sin(cfg.omega_drive * t)
        forcing = np.zeros(n, dtype=float)
        forcing[src] = drive
        x_next = 2.0 * x_curr - x_prev + cfg.c2 * (lap @ x_curr) - cfg.damp * (x_curr - x_prev) + forcing
        x_next = np.tanh(1.3 * x_next)
        x_prev, x_curr = x_curr, x_next
        if t >= cfg.burn_lock:
            c = math.cos(cfg.omega_drive * t)
            s = math.sin(cfg.omega_drive * t)
            re_acc += x_curr * c
            im_acc -= x_curr * s
            cnt += 1

    return (re_acc + 1j * im_acc) / max(cnt, 1)


def distance_response(z: np.ndarray, src: int, dmax: int) -> np.ndarray:
    n = len(z)
    out = []
    for d in range(1, min(dmax, n // 2) + 1):
        j1 = (src + d) % n
        j2 = (src - d) % n
        out.append(0.5 * (z[j1] + z[j2]))
    return np.array(out, dtype=complex)


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
    return np.fft.ifft(Y * H)


def fit_phase_slope_complex(resp: np.ndarray) -> Tuple[float, float, float]:
    d = np.arange(1, len(resp) + 1, dtype=float)
    amp = np.abs(resp)
    ph = np.unwrap(np.angle(resp))
    m = amp >= np.quantile(amp, 0.35)
    if np.sum(m) < 5:
        m = amp >= np.quantile(amp, 0.20)
    x = d[m]
    y = ph[m]
    A = np.column_stack([x, np.ones_like(x)])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    slope, intercept = float(coef[0]), float(coef[1])
    rmse = float(np.sqrt(np.mean((y - (slope * x + intercept)) ** 2)))
    return slope, intercept, rmse


def fit_phase_slope_hilbert(y_real: np.ndarray) -> Tuple[float, float, float]:
    d = np.arange(1, len(y_real) + 1, dtype=float)
    z = analytic_signal_real(y_real)
    amp = np.abs(z)
    ph = np.unwrap(np.angle(z))
    m = amp >= np.quantile(amp, 0.35)
    if np.sum(m) < 5:
        m = amp >= np.quantile(amp, 0.20)
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


def rankdata(x: np.ndarray) -> np.ndarray:
    idx = np.argsort(x)
    out = np.empty_like(idx, dtype=float)
    out[idx] = np.arange(len(x), dtype=float)
    return out


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3:
        return float("nan")
    rx = rankdata(x)
    ry = rankdata(y)
    c = np.corrcoef(rx, ry)
    return float(c[0, 1])


def main() -> None:
    rows: List[Dict[str, object]] = []
    for n in [96, 120]:
        for k in range(8):
            seed = 175800 + 100 * n + k
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
                omega_drive=0.20 + 0.05 * (k % 5),
                drive_amp=0.9 + 0.3 * (k % 2),
            )
            rng = np.random.default_rng(seed)
            w = build_w(cfg, rng)

            srcs = [seed % n, (seed * 3 + 7) % n]
            for sidx, src in enumerate(srcs):
                # beta observable
                ts = impulse_response_timeseries(cfg, w, src=src)
                yb = distance_tail_observable(ts, src=src, tail_start=cfg.tail_start, tail_end=cfg.tail_end, dmax=cfg.dmax)
                beta_hat, rmse_beta, delta_bic = fit_beta_envelope(yb)

                # omega/phi observable
                z = lockin_response(cfg, w, src=src)
                resp = distance_response(z, src=src, dmax=cfg.dmax)
                slope_c, phi_c, rmse_c = fit_phase_slope_complex(resp)
                slope_h, phi_h, rmse_h = fit_phase_slope_hilbert(np.real(resp))
                omega_hat = omega_from_phase_slopes(slope_c, slope_h)
                phi_hat = float(0.5 * (phi_c + phi_h))

                beta_grid = np.linspace(0.0, 0.20, 9)
                omega_pert = []
                for b in beta_grid:
                    rp = beta_envelope_perturb(resp, float(b))
                    sc2, _pc2, _rc2 = fit_phase_slope_complex(rp)
                    sh2, _ph2, _rh2 = fit_phase_slope_hilbert(np.real(rp))
                    omega_pert.append(omega_from_phase_slopes(sc2, sh2))
                omega_pert = np.array(omega_pert, dtype=float)
                delta_omega_beta_pert = float(np.nanmax(np.abs(omega_pert - omega_hat)))

                amp = np.abs(resp)
                snr_like = float(np.max(amp) / max(np.median(np.abs(amp - np.median(amp))), 1e-6))

                rows.append(
                    {
                        "n": n,
                        "seed": seed,
                        "source_index": sidx,
                        "source_node": int(src),
                        "omega_drive": cfg.omega_drive,
                        "beta_hat": float(beta_hat),
                        "beta_rmse": float(rmse_beta),
                        "beta_delta_bic_null_minus_beta": float(delta_bic),
                        "omega_hat": float(omega_hat),
                        "phi_hat": float(phi_hat),
                        "phase_rmse_complex": float(rmse_c),
                        "phase_rmse_hilbert": float(rmse_h),
                        "delta_omega_beta_perturb": float(delta_omega_beta_pert),
                        "snr_like": snr_like,
                    }
                )

    gf = PREREG["good_row_filter"]
    good = [
        r
        for r in rows
        if np.isfinite(r["omega_hat"])
        and r["snr_like"] >= gf["snr_like_min"]
        and r["beta_rmse"] <= gf["beta_rmse_max"]
        and r["phase_rmse_complex"] <= gf["phase_rmse_complex_max"]
        and r["phase_rmse_hilbert"] <= gf["phase_rmse_hilbert_max"]
    ]

    if len(good) >= 8:
        beta = np.array([r["beta_hat"] for r in good], dtype=float)
        omega = np.array([r["omega_hat"] for r in good], dtype=float)
        phi = np.array([r["phi_hat"] for r in good], dtype=float)
        rmse_b = np.array([r["beta_rmse"] for r in good], dtype=float)
        rmse_c = np.array([r["phase_rmse_complex"] for r in good], dtype=float)
        rmse_h = np.array([r["phase_rmse_hilbert"] for r in good], dtype=float)
        dbic = np.array([r["beta_delta_bic_null_minus_beta"] for r in good], dtype=float)
        dpert = np.array([r["delta_omega_beta_perturb"] for r in good], dtype=float)

        beta_iqr = float(np.quantile(beta, 0.75) - np.quantile(beta, 0.25))
        omega_iqr = float(np.quantile(omega, 0.75) - np.quantile(omega, 0.25))
        beta_ci = [float(np.quantile(beta, 0.025)), float(np.quantile(beta, 0.975))]
        omega_ci = [float(np.quantile(omega, 0.025)), float(np.quantile(omega, 0.975))]
        phi_cmean = float(math.atan2(np.mean(np.sin(phi)), np.mean(np.cos(phi))))
        phi_cstd = float(np.sqrt(-2.0 * np.log(max(np.sqrt(np.mean(np.sin(phi)) ** 2 + np.mean(np.cos(phi)) ** 2), 1e-12))))

        rho_bo = spearman(beta, omega)
        beta_nonboundary_rate = float(np.mean((beta > 5e-4) & (beta < 0.23)))
        omega_nonboundary_rate = float(np.mean((omega > 0.10) & (omega < 1.6)))
        beta_evidence_positive_rate = float(np.mean(dbic > 0.0))
        med_dbic = float(np.median(dbic))
        q90_dpert = float(np.quantile(dpert, 0.9))

        summary = {
            "n_good": len(good),
            "beta_median": float(np.median(beta)),
            "beta_iqr": beta_iqr,
            "beta_ci95_empirical": beta_ci,
            "omega_median": float(np.median(omega)),
            "omega_iqr": omega_iqr,
            "omega_ci95_empirical": omega_ci,
            "phi_circular_mean": phi_cmean,
            "phi_circular_std": phi_cstd,
            "median_beta_rmse": float(np.median(rmse_b)),
            "median_phase_rmse_complex": float(np.median(rmse_c)),
            "median_phase_rmse_hilbert": float(np.median(rmse_h)),
            "median_beta_delta_bic_null_minus_beta": med_dbic,
            "beta_evidence_positive_rate": beta_evidence_positive_rate,
            "beta_nonboundary_rate": beta_nonboundary_rate,
            "omega_nonboundary_rate": omega_nonboundary_rate,
            "spearman_beta_vs_omega": rho_bo,
            "delta_omega_beta_perturb_q90": q90_dpert,
        }
    else:
        summary = {
            "n_good": len(good),
            "beta_median": None,
            "beta_iqr": None,
            "beta_ci95_empirical": [None, None],
            "omega_median": None,
            "omega_iqr": None,
            "omega_ci95_empirical": [None, None],
            "phi_circular_mean": None,
            "phi_circular_std": None,
            "median_beta_rmse": None,
            "median_phase_rmse_complex": None,
            "median_phase_rmse_hilbert": None,
            "median_beta_delta_bic_null_minus_beta": None,
            "beta_evidence_positive_rate": None,
            "beta_nonboundary_rate": None,
            "omega_nonboundary_rate": None,
            "spearman_beta_vs_omega": None,
            "delta_omega_beta_perturb_q90": None,
        }

    th = PREREG["pass_thresholds"]
    pass_n = summary["n_good"] >= th["n_good_min"]
    pass_beta_evidence = (
        summary["median_beta_delta_bic_null_minus_beta"] is not None
        and summary["median_beta_delta_bic_null_minus_beta"] >= th["beta_evidence_median_dbic_min"]
        and summary["beta_evidence_positive_rate"] >= th["beta_evidence_positive_rate_min"]
    )
    pass_beta_nonboundary = summary["beta_nonboundary_rate"] is not None and summary["beta_nonboundary_rate"] >= th["beta_nonboundary_rate_min"]
    pass_omega_nonboundary = summary["omega_nonboundary_rate"] is not None and summary["omega_nonboundary_rate"] >= th["omega_nonboundary_rate_min"]
    pass_beta_spread = summary["beta_iqr"] is not None and summary["beta_iqr"] <= th["beta_iqr_max"]
    pass_omega_spread = summary["omega_iqr"] is not None and summary["omega_iqr"] <= th["omega_iqr_max"]
    pass_coupling = summary["spearman_beta_vs_omega"] is not None and abs(summary["spearman_beta_vs_omega"]) <= th["rho_beta_omega_abs_max"]
    pass_pert = summary["delta_omega_beta_perturb_q90"] is not None and summary["delta_omega_beta_perturb_q90"] <= th["delta_omega_beta_pert_q90_max"]
    pass_fit = (
        summary["median_phase_rmse_complex"] is not None
        and summary["median_phase_rmse_complex"] <= th["median_phase_rmse_complex_max"]
        and summary["median_phase_rmse_hilbert"] <= th["median_phase_rmse_hilbert_max"]
    )

    if all([pass_n, pass_beta_evidence, pass_beta_nonboundary, pass_omega_nonboundary, pass_beta_spread, pass_omega_spread, pass_coupling, pass_pert, pass_fit]):
        verdict = "PREREGISTERED_IDENTIFIABILITY_SUPPORTED"
    elif pass_n and pass_fit and pass_coupling and (pass_beta_nonboundary or pass_omega_nonboundary):
        verdict = "PREREGISTERED_IDENTIFIABILITY_PARTIAL"
    else:
        verdict = "PREREGISTERED_IDENTIFIABILITY_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "preregistered_protocol": PREREG,
        "n_total_rows": len(rows),
        "summary": summary,
        "pass_flags": {
            "n_good_enough": bool(pass_n),
            "beta_evidence": bool(pass_beta_evidence),
            "beta_nonboundary": bool(pass_beta_nonboundary),
            "omega_nonboundary": bool(pass_omega_nonboundary),
            "beta_spread": bool(pass_beta_spread),
            "omega_spread": bool(pass_omega_spread),
            "beta_omega_coupling_low": bool(pass_coupling),
            "beta_perturb_stability": bool(pass_pert),
            "fit_quality": bool(pass_fit),
        },
        "verdict": verdict,
        "rows": rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1758: PREREGISTERED IDENTIFIABILITY PROTOCOL",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Total rows: {len(rows)}",
        f"- Good rows: {summary['n_good']}",
        f"- beta median/CI95: {summary['beta_median']} / {summary['beta_ci95_empirical']}",
        f"- omega median/CI95: {summary['omega_median']} / {summary['omega_ci95_empirical']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- n_good_enough: {pass_n}",
        f"- beta_evidence: {pass_beta_evidence}",
        f"- beta_nonboundary: {pass_beta_nonboundary}",
        f"- omega_nonboundary: {pass_omega_nonboundary}",
        f"- beta_spread: {pass_beta_spread}",
        f"- omega_spread: {pass_omega_spread}",
        f"- beta_omega_coupling_low: {pass_coupling}",
        f"- beta_perturb_stability: {pass_pert}",
        f"- fit_quality: {pass_fit}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1758] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1758] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
