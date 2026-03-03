#!/usr/bin/env python3
"""
QW-1770: Kernel bridge from non-envelope branch.

Objective:
- Map leakage-controlled non-envelope observables to kernel-side parameters:
  beta_bridge, omega, phi.
- Use grouped CV by topology to avoid leakage.
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
OUT_JSON = ROOT / "report_qw1770_kernel_bridge_from_nonenvelope.json"
OUT_MD = ROOT / "RAPORT_QW1770_KERNEL_BRIDGE_FROM_NONENVELOPE.md"


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
    mem_mu: float
    mem_rate: float
    omega_drive: float
    n_steps: int = 700


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
    th = math.pi * np.tanh(smooth_periodic(rng.normal(size=n), iters=8) / 1.22)
    q = rng.integers(-2, 3, size=n).astype(float)
    w = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            d = ring_dist(i, j, n)
            amp = math.exp(-d / cfg.xi)
            if rng.random() < cfg.p_short / (d ** cfg.eta):
                amp += 0.8 * (d ** -0.95) * (1.0 + 0.13 * rng.normal())
            dt = abs(th[i] - th[j]) / math.pi
            dq = abs(q[i] - q[j])
            tors = dq + 0.35 * dt
            sgn = 1.0 if (math.cos(th[i] - th[j]) + 0.2 * (q[i] - q[j])) >= 0 else -1.0
            sym = sgn * amp * math.exp(-cfg.tau * tors)
            anti = cfg.rho * math.sin(th[i] - th[j]) / (d ** 1.1)
            w[i, j] = sym + anti
            w[j, i] = sym - anti
    smax = float(np.linalg.norm(w, 2))
    if smax > 1e-12:
        w /= smax
    return w


def run_driven(cfg: Cfg, w: np.ndarray, src: int, drive_amp: float) -> Tuple[np.ndarray, np.ndarray]:
    n = cfg.n
    lap = w - np.eye(n)
    x_prev = np.zeros(n, dtype=float)
    x_curr = np.zeros(n, dtype=float)
    mem = np.zeros(n, dtype=float)
    xs = np.zeros(cfg.n_steps, dtype=float)
    fs = np.zeros(cfg.n_steps, dtype=float)
    for t in range(cfg.n_steps):
        f = drive_amp * math.sin(cfg.omega_drive * t)
        forcing = np.zeros(n, dtype=float)
        forcing[src] = f
        mem = (1.0 - cfg.mem_rate) * mem + cfg.mem_rate * np.abs(x_curr)
        damping = cfg.mem_mu * mem * x_curr
        x_next = 2.0 * x_curr - x_prev + cfg.c2 * (lap @ x_curr) - cfg.damp * (x_curr - x_prev) - damping + forcing
        x_next = np.tanh(1.30 * x_next)
        x_prev, x_curr = x_curr, x_next
        xs[t] = x_curr[src]
        fs[t] = f
    return xs, fs


def cycle_peaks(x: np.ndarray, omega: float) -> np.ndarray:
    cl = max(int(round(2.0 * math.pi / max(omega, 1e-6))), 4)
    ncyc = len(x) // cl
    if ncyc < 4:
        return np.array([], dtype=float)
    vals = []
    for i in range(ncyc):
        seg = x[i * cl : (i + 1) * cl]
        vals.append(float(np.max(np.abs(seg))))
    return np.array(vals, dtype=float)


def adaptation_slope(x: np.ndarray, omega: float) -> float:
    p = cycle_peaks(x, omega)
    if len(p) < 6:
        return 0.0
    idx = np.arange(len(p), dtype=float)
    q1 = max(1, len(p) // 4)
    q3 = max(q1 + 2, (3 * len(p)) // 4)
    m = np.zeros(len(p), dtype=bool)
    m[q1:q3] = True
    xx = idx[m]
    yy = p[m]
    A = np.column_stack([xx, np.ones_like(xx)])
    coef, *_ = np.linalg.lstsq(A, yy, rcond=None)
    den = max(float(np.mean(yy)), 1e-8)
    return float(coef[0] / den)


def lockin_gain_phase(x: np.ndarray, omega: float, burn_frac: float = 0.40) -> Tuple[float, float]:
    t = np.arange(len(x), dtype=float)
    b = int(len(x) * burn_frac)
    xx = x[b:]
    tt = t[b:]
    c = np.cos(omega * tt)
    s = np.sin(omega * tt)
    re = float(np.mean(xx * c))
    im = float(-np.mean(xx * s))
    amp = float(np.sqrt(re * re + im * im))
    phase = float(math.atan2(im, re))
    return amp, phase


def hysteresis_area(x: np.ndarray, f: np.ndarray, omega: float, burn_frac: float = 0.40) -> float:
    b = int(len(x) * burn_frac)
    xx = x[b:]
    ff = f[b:]
    if len(xx) < 4:
        return 0.0
    area = 0.5 * abs(float(np.sum(ff[:-1] * xx[1:] - ff[1:] * xx[:-1])))
    ncyc = max(int((len(xx) * omega) / (2.0 * math.pi)), 1)
    return float(area / ncyc)


def wrap_phase(x: float) -> float:
    return float((x + math.pi) % (2.0 * math.pi) - math.pi)


def make_features(cfg: Cfg, w: np.ndarray, src: int) -> Dict[str, float]:
    x_lo, f_lo = run_driven(cfg, w, src=src, drive_amp=0.8)
    x_hi, f_hi = run_driven(cfg, w, src=src, drive_amp=1.4)
    g_lo, ph_lo = lockin_gain_phase(x_lo, cfg.omega_drive)
    g_hi, ph_hi = lockin_gain_phase(x_hi, cfg.omega_drive)
    h_lo = hysteresis_area(x_lo, f_lo, cfg.omega_drive)
    h_hi = hysteresis_area(x_hi, f_hi, cfg.omega_drive)
    a_lo = adaptation_slope(x_lo, cfg.omega_drive)
    a_hi = adaptation_slope(x_hi, cfg.omega_drive)
    return {
        "gain_low": float(g_lo),
        "gain_high": float(g_hi),
        "gain_compression": float(g_hi / max(g_lo, 1e-8)),
        "phase_low": float(ph_lo),
        "phase_shift": float(wrap_phase(ph_hi - ph_lo)),
        "hysteresis_ratio": float(h_hi / max(h_lo, 1e-8)),
        "adapt_low": float(a_lo),
        "adapt_diff": float(a_hi - a_lo),
    }


def ridge_fit_predict(Xtr: np.ndarray, ytr: np.ndarray, Xte: np.ndarray, lam: float = 0.8) -> np.ndarray:
    n, p = Xtr.shape
    Xa = np.column_stack([Xtr, np.ones(n)])
    I = np.eye(p + 1)
    I[-1, -1] = 0.0
    w = np.linalg.pinv(Xa.T @ Xa + lam * I) @ (Xa.T @ ytr)
    Xea = np.column_stack([Xte, np.ones(Xte.shape[0])])
    return Xea @ w


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


def r2_score(y: np.ndarray, yhat: np.ndarray) -> float:
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    if ss_tot <= 1e-15:
        return float("nan")
    return float(1.0 - ss_res / ss_tot)


def circular_mae(phi_true: np.ndarray, phi_pred: np.ndarray) -> float:
    d = np.array([wrap_phase(float(a - b)) for a, b in zip(phi_true, phi_pred)], dtype=float)
    return float(np.mean(np.abs(d)))


def grouped_kfold_indices(group_ids: np.ndarray, k: int = 5, seed: int = 1770) -> List[np.ndarray]:
    rng = np.random.default_rng(seed)
    groups = np.array(sorted(set(group_ids.tolist())), dtype=int)
    rng.shuffle(groups)
    chunks = [groups[i::k] for i in range(k)]
    folds = []
    for i in range(k):
        gset = set(chunks[i].tolist())
        idx = np.array([j for j, g in enumerate(group_ids.tolist()) if g in gset], dtype=int)
        folds.append(np.sort(idx))
    return folds


def main() -> None:
    rows: List[Dict[str, object]] = []
    mu_schedule = [0.0, 0.0, 0.0, 0.08, 0.14, 0.20]

    for n in [96, 120]:
        for k in range(18):
            seed = 177000 + 100 * n + k
            mu_target = mu_schedule[k % len(mu_schedule)]
            base_cfg = Cfg(
                n=n,
                seed=seed,
                xi=1.55 + 0.22 * ((k % 3) - 1),
                p_short=0.10 + 0.05 * (k % 2),
                eta=1.45 + 0.18 * ((k + 1) % 3),
                tau=0.45 + 0.16 * (k % 3),
                rho=0.18 + 0.10 * (k % 2),
                c2=0.30 + 0.08 * (k % 3),
                damp=0.03 + 0.02 * (k % 2),
                mem_mu=0.0,
                mem_rate=[0.02, 0.04, 0.06][k % 3],
                omega_drive=0.21 + 0.06 * (k % 4),
                n_steps=700,
            )
            rng = np.random.default_rng(seed)
            w = build_w(base_cfg, rng)
            srcs = [seed % n, (seed * 3 + 11) % n]
            for sidx, src in enumerate(srcs):
                cfg0 = Cfg(**{**base_cfg.__dict__, "mem_mu": 0.0})
                cfg1 = Cfg(**{**base_cfg.__dict__, "mem_mu": float(mu_target)})
                f0 = make_features(cfg0, w, src=src)
                f1 = make_features(cfg1, w, src=src)

                row = {
                    "n": n,
                    "seed": seed,
                    "topology_group": int(seed),
                    "source_index": sidx,
                    "source_node": int(src),
                    "beta_bridge_true": float(mu_target),
                    "omega_true": float(base_cfg.omega_drive),
                    "phi_true": float(f0["phase_low"]),
                    "d_gain_low": float(f1["gain_low"] - f0["gain_low"]),
                    "d_gain_high": float(f1["gain_high"] - f0["gain_high"]),
                    "d_gain_compression": float(f1["gain_compression"] - f0["gain_compression"]),
                    "d_phase_low": float(wrap_phase(f1["phase_low"] - f0["phase_low"])),
                    "d_phase_shift": float(wrap_phase(f1["phase_shift"] - f0["phase_shift"])),
                    "d_hysteresis_ratio": float(f1["hysteresis_ratio"] - f0["hysteresis_ratio"]),
                    "d_adapt_low": float(f1["adapt_low"] - f0["adapt_low"]),
                    "d_adapt_diff": float(f1["adapt_diff"] - f0["adapt_diff"]),
                    "b_gain_low": float(f0["gain_low"]),
                    "b_gain_high": float(f0["gain_high"]),
                    "b_gain_compression": float(f0["gain_compression"]),
                    "b_phase_low": float(f0["phase_low"]),
                    "b_phase_shift": float(f0["phase_shift"]),
                    "b_hysteresis_ratio": float(f0["hysteresis_ratio"]),
                    "b_adapt_low": float(f0["adapt_low"]),
                    "b_adapt_diff": float(f0["adapt_diff"]),
                }
                rows.append(row)

    feats = [
        "d_gain_low",
        "d_gain_high",
        "d_gain_compression",
        "d_phase_low",
        "d_phase_shift",
        "d_hysteresis_ratio",
        "d_adapt_low",
        "d_adapt_diff",
        "b_gain_low",
        "b_gain_high",
        "b_gain_compression",
        "b_phase_low",
        "b_phase_shift",
        "b_hysteresis_ratio",
        "b_adapt_low",
        "b_adapt_diff",
    ]
    X = np.array([[float(r[f]) for f in feats] for r in rows], dtype=float)
    yb = np.array([float(r["beta_bridge_true"]) for r in rows], dtype=float)
    yo = np.array([float(r["omega_true"]) for r in rows], dtype=float)
    yp = np.array([float(r["phi_true"]) for r in rows], dtype=float)
    group_ids = np.array([int(r["topology_group"]) for r in rows], dtype=int)

    folds = grouped_kfold_indices(group_ids, k=5, seed=1770)
    pb = np.zeros(len(rows), dtype=float)
    po = np.zeros(len(rows), dtype=float)
    pp = np.zeros(len(rows), dtype=float)

    fold_r2_beta = []
    fold_r2_omega = []
    fold_phi_mae = []

    for test_idx in folds:
        mask = np.ones(len(rows), dtype=bool)
        mask[test_idx] = False
        tr = np.where(mask)[0]
        te = test_idx

        Xtr = X[tr]
        Xte = X[te]
        m = np.mean(Xtr, axis=0)
        s = np.std(Xtr, axis=0)
        s = np.where(s < 1e-8, 1.0, s)
        Xtrn = (Xtr - m) / s
        Xten = (Xte - m) / s

        # beta and omega direct regression
        pb_te = ridge_fit_predict(Xtrn, yb[tr], Xten, lam=0.8)
        po_te = ridge_fit_predict(Xtrn, yo[tr], Xten, lam=0.8)
        pb[te] = np.clip(pb_te, 0.0, 0.25)
        po[te] = np.clip(po_te, 0.10, 1.6)

        # phi via sin/cos regression
        ys = np.sin(yp[tr])
        yc = np.cos(yp[tr])
        ps = ridge_fit_predict(Xtrn, ys, Xten, lam=0.8)
        pc = ridge_fit_predict(Xtrn, yc, Xten, lam=0.8)
        pp_te = np.array([wrap_phase(float(math.atan2(sv, cv))) for sv, cv in zip(ps, pc)], dtype=float)
        pp[te] = pp_te

        fold_r2_beta.append(r2_score(yb[te], pb[te]))
        fold_r2_omega.append(r2_score(yo[te], po[te]))
        fold_phi_mae.append(circular_mae(yp[te], pp[te]))

    r2_beta = r2_score(yb, pb)
    r2_omega = r2_score(yo, po)
    rho_beta = spearman(yb, pb)
    rho_omega = spearman(yo, po)
    phi_mae = circular_mae(yp, pp)
    phi_cerr_q90 = float(np.quantile(np.abs(np.array([wrap_phase(float(a - b)) for a, b in zip(yp, pp)], dtype=float)), 0.9))

    beta_nonboundary = float(np.mean((pb > 0.01) & (pb < 0.23)))
    omega_nonboundary = float(np.mean((po > 0.10) & (po < 1.6)))

    fb = np.array([x for x in fold_r2_beta if np.isfinite(x)], dtype=float)
    fo = np.array([x for x in fold_r2_omega if np.isfinite(x)], dtype=float)
    fp = np.array([x for x in fold_phi_mae if np.isfinite(x)], dtype=float)

    st_beta = float(np.std(fb)) if len(fb) else float("nan")
    st_omega = float(np.std(fo)) if len(fo) else float("nan")
    st_phi = float(np.std(fp)) if len(fp) else float("nan")

    pass_beta = np.isfinite(r2_beta) and r2_beta >= 0.60 and np.isfinite(rho_beta) and rho_beta >= 0.75
    pass_omega = np.isfinite(r2_omega) and r2_omega >= 0.70 and np.isfinite(rho_omega) and rho_omega >= 0.75
    pass_phi = phi_mae <= 0.55 and phi_cerr_q90 <= 0.90
    pass_nonboundary = beta_nonboundary >= 0.70 and omega_nonboundary >= 0.95
    pass_stability = (
        np.isfinite(st_beta)
        and np.isfinite(st_omega)
        and np.isfinite(st_phi)
        and st_beta <= 0.20
        and st_omega <= 0.18
        and st_phi <= 0.18
    )

    if pass_beta and pass_omega and pass_phi and pass_nonboundary and pass_stability:
        verdict = "KERNEL_BRIDGE_FROM_NONENVELOPE_SUPPORTED"
    elif pass_beta and pass_omega and (pass_phi or pass_nonboundary):
        verdict = "KERNEL_BRIDGE_FROM_NONENVELOPE_PARTIAL"
    else:
        verdict = "KERNEL_BRIDGE_FROM_NONENVELOPE_WEAK"

    summary = {
        "n_rows": len(rows),
        "n_groups": int(len(set(group_ids.tolist()))),
        "r2_beta_bridge": float(r2_beta),
        "rho_beta_bridge": float(rho_beta),
        "r2_omega": float(r2_omega),
        "rho_omega": float(rho_omega),
        "phi_circular_mae": float(phi_mae),
        "phi_circular_error_q90": float(phi_cerr_q90),
        "beta_pred_nonboundary_rate": beta_nonboundary,
        "omega_pred_nonboundary_rate": omega_nonboundary,
        "fold_r2_beta_std": st_beta,
        "fold_r2_omega_std": st_omega,
        "fold_phi_mae_std": st_phi,
        "beta_pred_ci95_empirical": [float(np.quantile(pb, 0.025)), float(np.quantile(pb, 0.975))],
        "omega_pred_ci95_empirical": [float(np.quantile(po, 0.025)), float(np.quantile(po, 0.975))],
    }

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "grouped_cv_by_topology": True,
            "targets": ["beta_bridge", "omega", "phi"],
            "feature_set": feats,
            "cv_folds": 5,
        },
        "summary": summary,
        "pass_flags": {
            "beta_bridge_regression": bool(pass_beta),
            "omega_regression": bool(pass_omega),
            "phi_regression": bool(pass_phi),
            "nonboundary_predictions": bool(pass_nonboundary),
            "cv_stability": bool(pass_stability),
        },
        "verdict": verdict,
        "rows": rows,
        "predictions": {
            "beta_bridge_pred": [float(v) for v in pb],
            "omega_pred": [float(v) for v in po],
            "phi_pred": [float(v) for v in pp],
        },
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1770: KERNEL BRIDGE FROM NONENVELOPE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Rows/groups: {summary['n_rows']} / {summary['n_groups']}",
        f"- R2 beta/omega: {summary['r2_beta_bridge']:.3f} / {summary['r2_omega']:.3f}",
        f"- Phi MAE / q90: {summary['phi_circular_mae']:.3f} / {summary['phi_circular_error_q90']:.3f}",
        f"- Nonboundary beta/omega: {summary['beta_pred_nonboundary_rate']:.3f} / {summary['omega_pred_nonboundary_rate']:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- beta_bridge_regression: {pass_beta}",
        f"- omega_regression: {pass_omega}",
        f"- phi_regression: {pass_phi}",
        f"- nonboundary_predictions: {pass_nonboundary}",
        f"- cv_stability: {pass_stability}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1770] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1770] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
