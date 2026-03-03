#!/usr/bin/env python3
"""
QW-1768: Leakage-controlled paired intervention.

Key rigor upgrade over QW-1766:
- Grouped CV by topology seed to prevent train/test leakage
  across multiple sources from the same network.
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
OUT_JSON = ROOT / "report_qw1768_leakage_controlled_paired_intervention.json"
OUT_MD = ROOT / "RAPORT_QW1768_LEAKAGE_CONTROLLED_PAIRED_INTERVENTION.md"


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


def sigmoid(z: np.ndarray) -> np.ndarray:
    x = np.clip(z, -35.0, 35.0)
    return 1.0 / (1.0 + np.exp(-x))


def train_logistic(X: np.ndarray, y: np.ndarray, lr: float = 0.08, epochs: int = 2600, l2: float = 1e-2) -> np.ndarray:
    n, p = X.shape
    Xa = np.column_stack([X, np.ones(n)])
    w = np.zeros(p + 1, dtype=float)
    for _ in range(epochs):
        pr = sigmoid(Xa @ w)
        grad = (Xa.T @ (pr - y)) / n
        reg = np.concatenate([l2 * w[:-1], np.array([0.0])])
        w -= lr * (grad + reg)
    return w


def predict_logistic(X: np.ndarray, w: np.ndarray) -> np.ndarray:
    Xa = np.column_stack([X, np.ones(X.shape[0])])
    return sigmoid(Xa @ w)


def ridge_fit_predict(Xtr: np.ndarray, ytr: np.ndarray, Xte: np.ndarray, lam: float = 0.6) -> np.ndarray:
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


def auc_score(y_true: np.ndarray, y_prob: np.ndarray) -> float:
    y = y_true.astype(int)
    p = y_prob.astype(float)
    n_pos = int(np.sum(y == 1))
    n_neg = int(np.sum(y == 0))
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    r = rankdata(p) + 1.0
    r_pos = float(np.sum(r[y == 1]))
    u = r_pos - n_pos * (n_pos + 1) / 2.0
    return float(u / (n_pos * n_neg))


def balanced_accuracy(y_true: np.ndarray, y_prob: np.ndarray, thr: float = 0.5) -> float:
    y = y_true.astype(int)
    yp = (y_prob >= thr).astype(int)
    pos = y == 1
    neg = y == 0
    tpr = float(np.mean(yp[pos] == 1)) if np.any(pos) else 0.0
    tnr = float(np.mean(yp[neg] == 0)) if np.any(neg) else 0.0
    return 0.5 * (tpr + tnr)


def r2_score(y: np.ndarray, yhat: np.ndarray) -> float:
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    if ss_tot <= 1e-15:
        return float("nan")
    return float(1.0 - ss_res / ss_tot)


def grouped_stratified_folds(group_ids: np.ndarray, group_label: Dict[int, int], k: int = 5, seed: int = 1768) -> List[np.ndarray]:
    rng = np.random.default_rng(seed)
    groups = np.array(sorted(group_label.keys()), dtype=int)
    pos = np.array([g for g in groups if group_label[g] == 1], dtype=int)
    neg = np.array([g for g in groups if group_label[g] == 0], dtype=int)
    rng.shuffle(pos)
    rng.shuffle(neg)
    pos_chunks = [pos[i::k] for i in range(k)]
    neg_chunks = [neg[i::k] for i in range(k)]
    folds = []
    for i in range(k):
        gtest = set(np.concatenate([pos_chunks[i], neg_chunks[i]]).tolist())
        idx = np.array([j for j, g in enumerate(group_ids.tolist()) if g in gtest], dtype=int)
        folds.append(np.sort(idx))
    return folds


def main() -> None:
    rows: List[Dict[str, object]] = []
    mu_schedule = [0.0, 0.0, 0.0, 0.08, 0.14, 0.20]

    for n in [96, 120]:
        for k in range(18):
            seed = 176800 + 100 * n + k
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
                    "mem_mu_target": float(mu_target),
                    "mem_on_target": 1 if mu_target > 0.02 else 0,
                    "d_gain_low": float(f1["gain_low"] - f0["gain_low"]),
                    "d_gain_high": float(f1["gain_high"] - f0["gain_high"]),
                    "d_gain_compression": float(f1["gain_compression"] - f0["gain_compression"]),
                    "d_phase_low": float(wrap_phase(f1["phase_low"] - f0["phase_low"])),
                    "d_phase_shift": float(wrap_phase(f1["phase_shift"] - f0["phase_shift"])),
                    "d_hysteresis_ratio": float(f1["hysteresis_ratio"] - f0["hysteresis_ratio"]),
                    "d_adapt_low": float(f1["adapt_low"] - f0["adapt_low"]),
                    "d_adapt_diff": float(f1["adapt_diff"] - f0["adapt_diff"]),
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
    ]
    X = np.array([[float(r[f]) for f in feats] for r in rows], dtype=float)
    y_bin = np.array([int(r["mem_on_target"]) for r in rows], dtype=int)
    y_mu = np.array([float(r["mem_mu_target"]) for r in rows], dtype=float)
    group_ids = np.array([int(r["topology_group"]) for r in rows], dtype=int)

    # group label taken from first row in each topology group
    group_label: Dict[int, int] = {}
    for r in rows:
        g = int(r["topology_group"])
        if g not in group_label:
            group_label[g] = int(r["mem_on_target"])

    folds = grouped_stratified_folds(group_ids, group_label, k=5, seed=1768)
    p_bin = np.zeros(len(rows), dtype=float)
    p_mu = np.zeros(len(rows), dtype=float)
    fold_auc = []
    fold_r2 = []

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

        wl = train_logistic(Xtrn, y_bin[tr], lr=0.08, epochs=2600, l2=1e-2)
        pb = predict_logistic(Xten, wl)
        p_bin[te] = pb

        pm = ridge_fit_predict(Xtrn, y_mu[tr], Xten, lam=0.6)
        p_mu[te] = np.clip(pm, 0.0, 0.25)

        fold_auc.append(auc_score(y_bin[te], pb))
        fold_r2.append(r2_score(y_mu[te], p_mu[te]))

    auc = auc_score(y_bin, p_bin)
    bacc = balanced_accuracy(y_bin, p_bin, thr=0.5)
    r2 = r2_score(y_mu, p_mu)
    rho = spearman(y_mu, p_mu)

    pos = y_bin == 1
    neg = y_bin == 0
    pos_nonboundary = float(np.mean((p_mu[pos] > 0.01) & (p_mu[pos] < 0.23))) if np.any(pos) else 0.0
    neg_low = float(np.mean(p_mu[neg] < 0.03)) if np.any(neg) else 0.0

    fold_auc_arr = np.array([x for x in fold_auc if np.isfinite(x)], dtype=float)
    fold_r2_arr = np.array([x for x in fold_r2 if np.isfinite(x)], dtype=float)
    auc_std = float(np.std(fold_auc_arr)) if len(fold_auc_arr) else float("nan")
    r2_std = float(np.std(fold_r2_arr)) if len(fold_r2_arr) else float("nan")

    pass_disc = np.isfinite(auc) and auc >= 0.78 and bacc >= 0.72
    pass_reg = np.isfinite(r2) and r2 >= 0.45 and np.isfinite(rho) and rho >= 0.70
    pass_cal = pos_nonboundary >= 0.70 and neg_low >= 0.70
    pass_stability = np.isfinite(auc_std) and np.isfinite(r2_std) and auc_std <= 0.12 and r2_std <= 0.22

    if pass_disc and pass_reg and pass_cal and pass_stability:
        verdict = "LEAKAGE_CONTROLLED_NONENVELOPE_SUPPORTED"
    elif pass_disc and (pass_reg or pass_cal):
        verdict = "LEAKAGE_CONTROLLED_NONENVELOPE_PARTIAL"
    else:
        verdict = "LEAKAGE_CONTROLLED_NONENVELOPE_WEAK"

    summary = {
        "n_rows": len(rows),
        "n_groups": len(group_label),
        "n_mem_on": int(np.sum(y_bin == 1)),
        "n_mem_off": int(np.sum(y_bin == 0)),
        "cv_auc_grouped": float(auc),
        "cv_balanced_accuracy_grouped": float(bacc),
        "cv_r2_mem_mu_grouped": float(r2),
        "cv_spearman_mem_mu_grouped": float(rho),
        "positive_nonboundary_rate": pos_nonboundary,
        "negative_low_mu_rate": neg_low,
        "fold_auc_std": auc_std,
        "fold_r2_std": r2_std,
        "pred_mem_mu_median": float(np.median(p_mu)),
        "pred_mem_mu_ci95_empirical": [float(np.quantile(p_mu, 0.025)), float(np.quantile(p_mu, 0.975))],
    }

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "grouped_cv_by_topology": True,
            "feature_set": feats,
            "cv_folds": 5,
            "target": "mem_mu",
            "no_envelope_fit_in_inference": True,
        },
        "summary": summary,
        "pass_flags": {
            "discrimination": bool(pass_disc),
            "continuous_regression": bool(pass_reg),
            "calibration_nonboundary": bool(pass_cal),
            "cv_stability": bool(pass_stability),
        },
        "verdict": verdict,
        "rows": rows,
        "predictions": {
            "mem_on_prob": [float(v) for v in p_bin],
            "mem_mu_pred": [float(v) for v in p_mu],
        },
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1768: LEAKAGE CONTROLLED PAIRED INTERVENTION",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Rows/groups: {summary['n_rows']} / {summary['n_groups']}",
        f"- CV-grouped AUC / BACC: {summary['cv_auc_grouped']:.3f} / {summary['cv_balanced_accuracy_grouped']:.3f}",
        f"- CV-grouped R2 / Spearman: {summary['cv_r2_mem_mu_grouped']:.3f} / {summary['cv_spearman_mem_mu_grouped']:.3f}",
        f"- Pos nonboundary / Neg low-mu rate: {summary['positive_nonboundary_rate']:.3f} / {summary['negative_low_mu_rate']:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- discrimination: {pass_disc}",
        f"- continuous_regression: {pass_reg}",
        f"- calibration_nonboundary: {pass_cal}",
        f"- cv_stability: {pass_stability}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1768] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1768] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
