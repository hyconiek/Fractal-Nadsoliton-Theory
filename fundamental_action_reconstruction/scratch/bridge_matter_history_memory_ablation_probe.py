#!/usr/bin/env python3
"""Scratch probe: memory-structure ablation on out-of-window prediction.

Compares memory families under identical holdout protocol:
- exp_1lag
- exp_2lag
- powerlaw_tail

Strictly evidence-level; no theorem claim.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
import scipy.optimize as so
import scipy.stats as ss
import sympy as sp

HERE = Path(__file__).resolve().parent
OUT = HERE / "bridge_matter_history_memory_ablation_report.json"


@dataclass
class Case:
    name: str
    train_min: float
    train_max: float
    hold_min: float
    hold_max: float
    n_train: int
    n_hold: int
    phase_shift: float
    eta_shift: float


def k_legacy(d: np.ndarray, alpha_geo: float, omega: float, phi: float, beta_tors: float) -> np.ndarray:
    return alpha_geo * np.cos(omega * d + phi) / (1.0 + beta_tors * d)


def k_strict(d: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d + phi) / (1.0 + beta * (d**eta))


def mem_exp_1lag(x: np.ndarray, rho: float) -> np.ndarray:
    y = np.empty_like(x)
    y[0] = x[0]
    a = 1.0 - rho
    for i in range(1, x.size):
        y[i] = a * x[i] + rho * y[i - 1]
    return y


def mem_exp_2lag(x: np.ndarray, rho1: float, rho2: float) -> np.ndarray:
    y = np.empty_like(x)
    y[:2] = x[:2]
    a = max(0.0, 1.0 - rho1 - rho2)
    for i in range(2, x.size):
        y[i] = a * x[i] + rho1 * y[i - 1] + rho2 * y[i - 2]
    return y


def mem_powerlaw(x: np.ndarray, nu: float, kmax: int = 32) -> np.ndarray:
    n = x.size
    y = np.zeros_like(x)
    w = 1.0 / (np.arange(1, kmax + 1, dtype=float) ** nu)
    w /= np.sum(w)
    for i in range(n):
        kk = min(kmax, i + 1)
        y[i] = float(np.dot(w[:kk], x[i - kk + 1 : i + 1][::-1]))
    return y


def fit_family(train_d, y_train, alpha_geo, omega, phi, family, seed):
    rng = np.random.default_rng(seed)

    def loss_exp1(v):
        dphi, bt, rho, gamma = v
        raw = k_legacy(train_d, alpha_geo, omega, phi + dphi, bt)
        mem = mem_exp_1lag(raw, rho)
        feat = np.vstack([raw, gamma * mem]).T
        c, *_ = np.linalg.lstsq(feat, y_train, rcond=None)
        pred = feat @ c
        return float(np.mean((pred - y_train) ** 2))

    def loss_exp2(v):
        dphi, bt, r1, r2, gamma = v
        raw = k_legacy(train_d, alpha_geo, omega, phi + dphi, bt)
        mem = mem_exp_2lag(raw, r1, r2)
        feat = np.vstack([raw, gamma * mem]).T
        c, *_ = np.linalg.lstsq(feat, y_train, rcond=None)
        pred = feat @ c
        return float(np.mean((pred - y_train) ** 2))

    def loss_pow(v):
        dphi, bt, nu, gamma = v
        raw = k_legacy(train_d, alpha_geo, omega, phi + dphi, bt)
        mem = mem_powerlaw(raw, nu)
        feat = np.vstack([raw, gamma * mem]).T
        c, *_ = np.linalg.lstsq(feat, y_train, rcond=None)
        pred = feat @ c
        return float(np.mean((pred - y_train) ** 2))

    if family == "exp_1lag":
        best = None
        for _ in range(10):
            x0 = np.array([rng.uniform(-np.pi, np.pi), rng.uniform(1e-4, 40.0), rng.uniform(1e-5, 0.95), rng.uniform(-8, 8)])
            res = so.minimize(loss_exp1, x0=x0, bounds=[(-np.pi, np.pi), (1e-5, 60.0), (1e-5, 0.999), (-10.0, 10.0)], method="L-BFGS-B")
            if best is None or res.fun < best.fun:
                best = res
        dphi, bt, rho, gamma = [float(v) for v in best.x]
        return {"family": family, "delta_phi": dphi, "beta_tors": bt, "rho": rho, "gamma": gamma}

    if family == "exp_2lag":
        best = None
        for _ in range(10):
            x0 = np.array([rng.uniform(-np.pi, np.pi), rng.uniform(1e-4, 40.0), rng.uniform(1e-5, 0.7), rng.uniform(1e-5, 0.25), rng.uniform(-8, 8)])
            res = so.minimize(loss_exp2, x0=x0, bounds=[(-np.pi, np.pi), (1e-5, 60.0), (1e-5, 0.95), (1e-5, 0.35), (-10.0, 10.0)], method="L-BFGS-B")
            if best is None or res.fun < best.fun:
                best = res
        dphi, bt, r1, r2, gamma = [float(v) for v in best.x]
        return {"family": family, "delta_phi": dphi, "beta_tors": bt, "rho1": r1, "rho2": r2, "gamma": gamma}

    if family == "powerlaw_tail":
        best = None
        for _ in range(10):
            x0 = np.array([rng.uniform(-np.pi, np.pi), rng.uniform(1e-4, 40.0), rng.uniform(0.2, 2.5), rng.uniform(-8, 8)])
            res = so.minimize(loss_pow, x0=x0, bounds=[(-np.pi, np.pi), (1e-5, 60.0), (0.1, 3.0), (-10.0, 10.0)], method="L-BFGS-B")
            if best is None or res.fun < best.fun:
                best = res
        dphi, bt, nu, gamma = [float(v) for v in best.x]
        return {"family": family, "delta_phi": dphi, "beta_tors": bt, "nu": nu, "gamma": gamma}

    raise ValueError(family)


def holdout_mse(d, y, params, alpha_geo, omega, phi):
    raw = k_legacy(d, alpha_geo, omega, phi + params["delta_phi"], params["beta_tors"])
    fam = params["family"]
    if fam == "exp_1lag":
        mem = mem_exp_1lag(raw, params["rho"])
    elif fam == "exp_2lag":
        mem = mem_exp_2lag(raw, params["rho1"], params["rho2"])
    else:
        mem = mem_powerlaw(raw, params["nu"])
    feat = np.vstack([raw, params["gamma"] * mem]).T
    coef, *_ = np.linalg.lstsq(feat, y, rcond=None)
    pred = feat @ coef
    return float(np.mean((pred - y) ** 2))


def main() -> None:
    alpha_geo = float(4.0 * np.log(2.0))
    omega0, phi0, beta0, eta0 = 0.18575, 0.16250, 1.0, 1.8
    cases = [
        Case("expand_right", 1.0, 8.0, 8.0, 14.0, 450, 350, 0.0, 0.0),
        Case("expand_left", 4.0, 11.0, 1.0, 4.0, 450, 280, 0.0, 0.0),
        Case("phase_eta_plus", 1.0, 8.0, 8.0, 14.0, 450, 350, 0.04, 0.10),
        Case("phase_eta_minus", 4.0, 11.0, 1.0, 4.0, 450, 280, -0.04, -0.10),
    ]
    families = ["exp_1lag", "exp_2lag", "powerlaw_tail"]
    seeds = [20260519 + i for i in range(6)]

    rows = []
    all_rank1 = []
    for c in cases:
        dtr = np.linspace(c.train_min, c.train_max, c.n_train)
        dho = np.linspace(c.hold_min, c.hold_max, c.n_hold)
        ytr = k_strict(dtr, omega0, phi0 + c.phase_shift, beta0, eta0 + c.eta_shift)
        yho = k_strict(dho, omega0, phi0 + c.phase_shift, beta0, eta0 + c.eta_shift)
        case_seed_rows = []
        for s in seeds:
            mses = {}
            for f in families:
                p = fit_family(dtr, ytr, alpha_geo, omega0, phi0 + c.phase_shift, f, s)
                mses[f] = holdout_mse(dho, yho, p, alpha_geo, omega0, phi0 + c.phase_shift)
            ranking = sorted(mses.items(), key=lambda kv: kv[1])
            all_rank1.append(ranking[0][0])
            case_seed_rows.append({"seed": s, "holdout_mse_by_family": mses, "rank": [k for k, _ in ranking]})

        # aggregate medians
        med = {f: float(np.median([r["holdout_mse_by_family"][f] for r in case_seed_rows])) for f in families}
        rows.append({"case": asdict(c), "seed_rows": case_seed_rows, "median_holdout_mse_by_family": med})

    winner_freq = {f: float(np.mean(np.array(all_rank1) == f)) for f in families}
    # pairwise sign tests against exp_1lag baseline
    deltas_21 = []
    deltas_pw = []
    for cr in rows:
        for sr in cr["seed_rows"]:
            mse = sr["holdout_mse_by_family"]
            deltas_21.append(mse["exp_1lag"] - mse["exp_2lag"])
            deltas_pw.append(mse["exp_1lag"] - mse["powerlaw_tail"])
    deltas_21 = np.array(deltas_21, dtype=float)
    deltas_pw = np.array(deltas_pw, dtype=float)
    p_21 = float(ss.wilcoxon(deltas_21, alternative="greater", zero_method="zsplit", method="approx").pvalue) if np.any(np.abs(deltas_21)>1e-15) else 1.0
    p_pw = float(ss.wilcoxon(deltas_pw, alternative="greater", zero_method="zsplit", method="approx").pvalue) if np.any(np.abs(deltas_pw)>1e-15) else 1.0

    z, r1, r2 = sp.symbols("z r1 r2", real=True)
    two_lag_tf = sp.simplify(1 / (1 - r1 / z - r2 / z**2))

    report = {
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "result_kind": "SCRATCH_MEMORY_ABLATION_OUT_OF_WINDOW__EVIDENCE_NOT_THEOREM",
        "cases": rows,
        "global": {
            "winner_frequency": winner_freq,
            "delta_exp1_minus_exp2_q05_q50_q95": [float(v) for v in np.quantile(deltas_21, [0.05,0.5,0.95])],
            "delta_exp1_minus_powerlaw_q05_q50_q95": [float(v) for v in np.quantile(deltas_pw, [0.05,0.5,0.95])],
            "wilcoxon_p_exp2_better_than_exp1": p_21,
            "wilcoxon_p_powerlaw_better_than_exp1": p_pw,
            "robust_new_family_gain_flag": bool((winner_freq["exp_2lag"] > 0.55 or winner_freq["powerlaw_tail"] > 0.55) and (p_21 < 0.01 or p_pw < 0.01)),
        },
        "symbolic_trace": {
            "two_lag_transfer_form": sp.sstr(two_lag_tf),
        },
        "hard_limits": [
            "Ablation ranking is model-selection evidence only.",
            "No theorem-level bridge is claimed.",
            "No legacy physical-role transfer to strict kernel.",
        ],
        "next_honest_step": "If no family shows robust holdout gain, classify matter-history memory route as non-bridge fit family and pivot effort to strict-lane global obstructions (renorm/unitarity/PO3/PO2/QW-2191).",
    }

    OUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
