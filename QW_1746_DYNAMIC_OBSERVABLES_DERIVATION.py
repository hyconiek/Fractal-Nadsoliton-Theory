#!/usr/bin/env python3
"""
QW-1746: Dynamic observables derivation to break omega-beta degeneracy.

Observables per run:
1) omega_phase: slope of distance-phase lag at driven frequency.
2) omega_zero: inferred from zero-cross spacing of real transfer response.
3) beta_env: envelope decay from |transfer(d)| ~ A/(1+beta*d).
4) phi_phase: intercept from phase regression.
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
OUT_JSON = ROOT / "report_qw1746_dynamic_observables_derivation.json"
OUT_MD = ROOT / "RAPORT_QW1746_DYNAMIC_OBSERVABLES_DERIVATION.md"


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


def smooth_periodic(x: np.ndarray, it: int = 7) -> np.ndarray:
    y = x.copy()
    for _ in range(it):
        y = 0.25 * np.roll(y, -1) + 0.5 * y + 0.25 * np.roll(y, 1)
    return y


def build_w(cfg: Cfg, rng: np.random.Generator) -> np.ndarray:
    n = cfg.n
    th = math.pi * np.tanh(smooth_periodic(rng.normal(size=n)) / 1.2)
    q = rng.integers(-2, 3, size=n).astype(float)

    w = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            d = ring_dist(i, j, n)
            amp = math.exp(-d / cfg.xi)
            if rng.random() < cfg.p_short / (d ** cfg.eta):
                amp += 0.95 * (d ** -0.95) * (1.0 + 0.15 * rng.normal())

            dt = abs(th[i] - th[j]) / math.pi
            dq = abs(q[i] - q[j])
            tors = dq + 0.38 * dt

            sgn = 1.0 if (math.cos(th[i] - th[j]) + 0.22 * (q[i] - q[j])) >= 0 else -1.0
            sym = sgn * amp * math.exp(-cfg.tau * tors)
            anti = cfg.rho * math.sin(th[i] - th[j]) / (d ** 1.1)
            w[i, j] = sym + anti
            w[j, i] = sym - anti

    smax = float(np.linalg.norm(w, 2))
    if smax > 1e-12:
        w /= smax
    return w


def lockin_response(cfg: Cfg, w: np.ndarray) -> np.ndarray:
    n = cfg.n
    lap = w - np.eye(n)
    x_prev = np.zeros(n, dtype=float)
    x_curr = np.zeros(n, dtype=float)
    src = 0
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


def distance_complex_response(z: np.ndarray, dmax: int = 12) -> np.ndarray:
    n = len(z)
    src = 0
    out = []
    for d in range(1, min(dmax, n // 2) + 1):
        j1 = (src + d) % n
        j2 = (src - d) % n
        out.append(0.5 * (z[j1] + z[j2]))
    return np.array(out, dtype=complex)


def fit_phase_slope(resp: np.ndarray) -> Tuple[float, float, float]:
    d = np.arange(1, len(resp) + 1, dtype=float)
    amp = np.abs(resp)
    ph = np.unwrap(np.angle(resp))
    # Keep informative points only.
    m = amp >= np.quantile(amp, 0.35)
    if np.sum(m) < 4:
        m = amp >= np.quantile(amp, 0.2)
    x = d[m]
    y = ph[m]
    A = np.column_stack([x, np.ones_like(x)])
    coeff, *_ = np.linalg.lstsq(A, y, rcond=None)
    slope, intercept = float(coeff[0]), float(coeff[1])
    pred = slope * x + intercept
    rmse = float(np.sqrt(np.mean((y - pred) ** 2)))
    return slope, intercept, rmse


def zero_cross_spacing(re_resp: np.ndarray) -> float:
    d = np.arange(1, len(re_resp) + 1, dtype=float)
    zc = []
    for i in range(len(re_resp) - 1):
        y1, y2 = re_resp[i], re_resp[i + 1]
        if y1 == 0:
            zc.append(float(d[i]))
        elif y1 * y2 < 0:
            # linear interpolation
            frac = abs(y1) / (abs(y1) + abs(y2))
            zc.append(float(d[i] + frac))
    if len(zc) < 2:
        return float("nan")
    spacings = np.diff(np.array(zc, dtype=float))
    return float(np.median(spacings))


def omega_from_spacing(spacing: float) -> float:
    # cos node spacing ~ pi/omega
    if not np.isfinite(spacing) or spacing <= 1e-6:
        return float("nan")
    return float(math.pi / spacing)


def fit_beta_envelope(resp: np.ndarray) -> Tuple[float, float]:
    d = np.arange(1, len(resp) + 1, dtype=float)
    a = np.abs(resp)
    if a[0] > 1e-12:
        a = a / a[0]
    beta_grid = np.linspace(0.001, 0.25, 2500)
    best_beta = 0.02
    best_rmse = float("inf")
    for b in beta_grid:
        e = 1.0 / (1.0 + b * d)
        s = float(np.dot(a, e) / max(np.dot(e, e), 1e-12))
        p = s * e
        rmse = float(np.sqrt(np.mean((a - p) ** 2)))
        if rmse < best_rmse:
            best_rmse = rmse
            best_beta = float(b)
    return best_beta, best_rmse


def circular_mean(x: np.ndarray) -> float:
    return float(np.arctan2(np.mean(np.sin(x)), np.mean(np.cos(x))))


def circular_std(x: np.ndarray) -> float:
    r = np.sqrt(np.mean(np.sin(x)) ** 2 + np.mean(np.cos(x)) ** 2)
    if r <= 1e-12:
        return float("inf")
    return float(np.sqrt(-2.0 * np.log(r)))


def main() -> None:
    runs: List[Dict[str, object]] = []
    for n in [96, 120]:
        for k in range(18):
            seed = 174600 + 100 * n + k
            cfg = Cfg(
                n=n,
                seed=seed,
                xi=1.6 + 0.25 * ((k % 3) - 1),
                p_short=0.10 + 0.05 * (k % 2),
                eta=1.45 + 0.2 * ((k + 1) % 3),
                tau=0.45 + 0.18 * (k % 3),
                rho=0.20 + 0.10 * (k % 2),
                c2=0.32 + 0.08 * (k % 3),
                damp=0.03 + 0.02 * (k % 2),
                omega_drive=0.18 + 0.05 * (k % 4),
                drive_amp=0.9 + 0.3 * (k % 2),
            )
            rng = np.random.default_rng(seed)
            w = build_w(cfg, rng)
            z = lockin_response(cfg, w)
            resp = distance_complex_response(z, dmax=12)

            phase_slope, phase_intercept, phase_rmse = fit_phase_slope(resp)
            spacing = zero_cross_spacing(np.real(resp))
            omega_zero = omega_from_spacing(spacing)
            beta_hat, beta_rmse = fit_beta_envelope(resp)

            amp = np.abs(resp)
            snr_like = float(np.max(amp) / max(np.median(np.abs(amp - np.median(amp))), 1e-6))

            runs.append(
                {
                    "n": n,
                    "seed": seed,
                    "omega_drive": cfg.omega_drive,
                    "omega_phase_hat": phase_slope,
                    "phi_phase_hat": phase_intercept,
                    "phase_fit_rmse": phase_rmse,
                    "zero_spacing": spacing,
                    "omega_zero_hat": omega_zero,
                    "beta_env_hat": beta_hat,
                    "beta_fit_rmse": beta_rmse,
                    "snr_like": snr_like,
                    "response_real": {str(i + 1): float(np.real(resp[i])) for i in range(len(resp))},
                    "response_abs": {str(i + 1): float(np.abs(resp[i])) for i in range(len(resp))},
                }
            )

    def is_good_strict(r: Dict[str, object]) -> bool:
        return (
            np.isfinite(r["omega_zero_hat"])
            and r["snr_like"] >= 4.0
            and r["phase_fit_rmse"] <= 0.45
            and r["beta_fit_rmse"] <= 0.25
        )

    def is_good_relaxed(r: Dict[str, object]) -> bool:
        return (
            np.isfinite(r["omega_zero_hat"])
            and r["snr_like"] >= 4.0
            and r["phase_fit_rmse"] <= 0.60
            and r["beta_fit_rmse"] <= 0.25
        )

    for r in runs:
        r["is_good_strict"] = bool(is_good_strict(r))
        r["is_good_relaxed"] = bool(is_good_relaxed(r))

    strict = [r for r in runs if r["is_good_strict"]]
    relaxed = [r for r in runs if r["is_good_relaxed"]]

    if len(strict) >= 12:
        chosen_label = "strict"
        good = strict
    elif len(relaxed) >= 10:
        chosen_label = "relaxed"
        good = relaxed
    else:
        chosen_label = "fallback_relaxed_small"
        good = relaxed

    om_p = np.array([r["omega_phase_hat"] for r in good], dtype=float) if good else np.array([], dtype=float)
    om_z = np.array([r["omega_zero_hat"] for r in good], dtype=float) if good else np.array([], dtype=float)
    be = np.array([r["beta_env_hat"] for r in good], dtype=float) if good else np.array([], dtype=float)
    ph = np.array([r["phi_phase_hat"] for r in good], dtype=float) if good else np.array([], dtype=float)

    if len(good) >= 8:
        summary = {
            "omega_phase_median": float(np.median(om_p)),
            "omega_phase_iqr": float(np.quantile(om_p, 0.75) - np.quantile(om_p, 0.25)),
            "omega_zero_median": float(np.median(om_z)),
            "omega_zero_iqr": float(np.quantile(om_z, 0.75) - np.quantile(om_z, 0.25)),
            "beta_median": float(np.median(be)),
            "beta_iqr": float(np.quantile(be, 0.75) - np.quantile(be, 0.25)),
            "phi_circular_mean": float(circular_mean(ph)),
            "phi_circular_std": float(circular_std(ph)),
            "n_good": len(good),
        }
        omega_consistency = abs(summary["omega_phase_median"] - summary["omega_zero_median"])
    else:
        summary = {
            "omega_phase_median": None,
            "omega_phase_iqr": None,
            "omega_zero_median": None,
            "omega_zero_iqr": None,
            "beta_median": None,
            "beta_iqr": None,
            "phi_circular_mean": None,
            "phi_circular_std": None,
            "n_good": len(good),
        }
        omega_consistency = float("inf")

    pass_good = len(good) >= 12
    pass_consistency = np.isfinite(omega_consistency) and omega_consistency <= 0.20
    pass_spread = (
        summary["omega_phase_iqr"] is not None
        and summary["omega_zero_iqr"] is not None
        and summary["beta_iqr"] is not None
        and summary["omega_phase_iqr"] <= 0.30
        and summary["omega_zero_iqr"] <= 0.30
        and summary["beta_iqr"] <= 0.10
    )

    if pass_good and pass_consistency and pass_spread:
        verdict = "DYNAMIC_OBSERVABLES_DERIVATION_INFORMATIVE"
    elif pass_good:
        verdict = "DYNAMIC_OBSERVABLES_DERIVATION_PARTIAL"
    else:
        verdict = "DYNAMIC_OBSERVABLES_DERIVATION_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "n_total_runs": len(runs),
            "quality_filter": {
                "snr_like_min": 4.0,
                "phase_fit_rmse_max_strict": 0.45,
                "phase_fit_rmse_max_relaxed": 0.60,
                "beta_fit_rmse_max": 0.25,
                "finite_omega_zero_required": True,
            },
        },
        "subset_sizes": {
            "strict": len(strict),
            "relaxed": len(relaxed),
            "chosen_label": chosen_label,
        },
        "summary_good_subset": summary,
        "omega_consistency_absdiff": omega_consistency,
        "pass_flags": {
            "n_good_enough": bool(pass_good),
            "omega_consistency": bool(pass_consistency),
            "spread_control": bool(pass_spread),
        },
        "verdict": verdict,
        "rows": runs,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1746: DYNAMIC OBSERVABLES DERIVATION",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Total runs: {len(runs)}",
        f"- Good runs: {summary['n_good']}",
        f"- omega consistency |phase-zero|: {omega_consistency}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- n_good_enough: {pass_good}",
        f"- omega_consistency: {pass_consistency}",
        f"- spread_control: {pass_spread}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1746] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1746] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
