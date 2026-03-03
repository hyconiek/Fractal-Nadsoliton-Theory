#!/usr/bin/env python3
"""
QW-1895: Joint Pareto gate under derivational constraints with anti-overfit audit.

Key audit:
- Intra-profile odd/even holdout prediction test (train on odd d, test on even d).
- Complexity check via parameter count vs observation count.

Compares:
- Model A: constrained baseline (QW-1891 style, linear channel size=2)
- Model B: amplitude-reformulated (QW-1894 style, linear channel size=5)
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1895_joint_pareto_gate_under_constraints.json"
OUT_MD = ROOT / "RAPORT_QW1895_JOINT_PARETO_GATE_UNDER_CONSTRAINTS.md"


def read_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def node_indicator(d: np.ndarray, nodes: List[int]) -> np.ndarray:
    out = np.zeros_like(d, dtype=float)
    node_set = {int(x) for x in nodes}
    for i, x in enumerate(d):
        out[i] = 1.0 if int(x) in node_set else 0.0
    return out


def signed_coupling_drive(d: np.ndarray, nodes: List[int], tau: float = 2.5) -> np.ndarray:
    out = np.zeros_like(d, dtype=float)
    for i, di in enumerate(d):
        s = 0.0
        for n in nodes:
            delta = float(di - n)
            sign = 1.0 if delta > 0 else (-1.0 if delta < 0 else 0.0)
            s += sign * math.exp(-abs(delta) / tau)
        out[i] = s
    m = float(np.max(np.abs(out)))
    if m > 1e-12:
        out = out / m
    return out


def simulate_state(
    d: np.ndarray,
    omega: float,
    phi: float,
    rho: float,
    xi: float,
    zeta: float,
    eta: float,
    nodes: List[int],
    tau: float,
) -> np.ndarray:
    ind = node_indicator(d, nodes)
    coupl = signed_coupling_drive(d, nodes, tau=tau)
    st = np.zeros_like(d, dtype=float)
    prev = 0.0
    for i, di in enumerate(d):
        drive = xi * math.sin(omega * float(di) + phi)
        node_term = -zeta * ind[i]
        cur = rho * prev + drive + node_term + eta * coupl[i]
        st[i] = cur
        prev = cur
    return st


def canon_score(omega: float, phi: float, beta: float, t: Dict[str, float]) -> float:
    z_o = abs(omega - t["omega"]) / 0.20
    dphi = (phi - t["phi"] + math.pi) % (2.0 * math.pi) - math.pi
    z_p = abs(dphi) / 0.30
    z_b = abs(beta - t["beta"]) / 0.030
    z = min(20.0, math.sqrt((z_o * z_o + z_p * z_p + z_b * z_b) / 3.0))
    return float(math.exp(-0.5 * z * z))


def design_matrix(
    d: np.ndarray,
    omega: float,
    phi: float,
    beta: float,
    rho: float,
    xi: float,
    zeta: float,
    eta: float,
    nodes: List[int],
    tau: float,
    model_kind: str,
) -> np.ndarray:
    env = np.cos(omega * d + phi) / (1.0 + beta * d)
    st = simulate_state(d, omega, phi, rho, xi, zeta, eta, nodes, tau)

    if model_kind == "base":
        return np.column_stack([env, st])
    if model_kind == "amp":
        return np.column_stack([np.ones_like(d), d, 1.0 / d, env, st])

    raise ValueError(model_kind)


def fit_linear(X: np.ndarray, y: np.ndarray, model_kind: str, lambda_amp: float) -> np.ndarray:
    if model_kind == "base":
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        return coef

    # ridge for first 3 amplitude baseline coefficients
    R = np.diag([lambda_amp, lambda_amp, lambda_amp, 1e-9, 1e-9])
    return np.linalg.solve(X.T @ X + R, X.T @ y)


def fit_nonlinear_on_train(
    y_train: np.ndarray,
    d_train: np.ndarray,
    t: Dict[str, float],
    nodes: List[int],
    tau: float,
    lambda_c: float,
    lambda_b: float,
    lambda_amp: float,
    model_kind: str,
    rng: np.random.Generator,
    n_samples: int,
) -> Dict:
    def penalty(omega: float, phi: float, beta: float, rho: float, xi: float, zeta: float, eta: float) -> float:
        p = 0.0
        if not (math.pi / 4 - 0.12 <= omega <= math.pi / 4 + 0.12):
            p += 1.0
        if min(abs(phi - math.pi / 6.0), abs(phi + math.pi / 6.0)) > 0.22:
            p += 1.0
        if not (0.005 <= beta <= 0.05):
            p += 1.0
        if not (0.35 < rho < 0.93):
            p += 0.5
        if not (0.03 < xi < 0.58):
            p += 0.5
        if not (0.03 < zeta < 0.55):
            p += 0.5
        if not (-0.45 < eta < 0.45):
            p += 0.5
        return p / 5.0

    def eval_one(omega: float, phi: float, beta: float, rho: float, xi: float, zeta: float, eta: float) -> Dict:
        X = design_matrix(d_train, omega, phi, beta, rho, xi, zeta, eta, nodes, tau, model_kind)
        coef = fit_linear(X, y_train, model_kind, lambda_amp)
        yhat = X @ coef
        rmse = float(np.sqrt(np.mean((y_train - yhat) ** 2)))
        cs = canon_score(omega, phi, beta, t)
        obj = rmse + lambda_c * (1.0 - cs) + lambda_b * penalty(omega, phi, beta, rho, xi, zeta, eta)
        return {
            "omega": omega,
            "phi": phi,
            "beta": beta,
            "rho": rho,
            "xi": xi,
            "zeta": zeta,
            "eta": eta,
            "coef": coef,
            "train_rmse": rmse,
            "canon_score": cs,
            "objective": obj,
        }

    anchors = [
        (math.pi / 4.0, math.pi / 6.0, 0.01, 0.60, 0.20, 0.20, 0.0),
        (math.pi / 4.0, -math.pi / 6.0, 0.01, 0.60, 0.20, 0.20, 0.0),
        (math.pi / 4.0 + 0.05, math.pi / 6.0, 0.02, 0.75, 0.20, 0.20, 0.2),
    ]

    best = None
    for a in anchors:
        c = eval_one(*a)
        if best is None or c["objective"] < best["objective"]:
            best = c

    for _ in range(n_samples):
        omega = float(np.clip(rng.normal(best["omega"], 0.03), math.pi / 4 - 0.14, math.pi / 4 + 0.14))
        phi_center = math.pi / 6.0 if rng.random() < 0.5 else -math.pi / 6.0
        phi = float(np.clip(rng.normal(phi_center, 0.11), -1.2, 1.2))
        beta = float(np.clip(rng.normal(best["beta"], 0.010), 0.003, 0.08))
        rho = float(np.clip(rng.normal(best["rho"], 0.10), 0.25, 0.98))
        xi = float(np.clip(rng.normal(best["xi"], 0.06), 0.01, 0.65))
        zeta = float(np.clip(rng.normal(best["zeta"], 0.06), 0.01, 0.70))
        eta = float(np.clip(rng.normal(best["eta"], 0.10), -0.70, 0.70))
        c = eval_one(omega, phi, beta, rho, xi, zeta, eta)
        if c["objective"] < best["objective"]:
            best = c

    return best


def odd_even_holdout_eval(
    y: np.ndarray,
    d: np.ndarray,
    t: Dict[str, float],
    nodes: List[int],
    tau: float,
    lambda_c: float,
    lambda_b: float,
    lambda_amp: float,
    model_kind: str,
    rng: np.random.Generator,
) -> Dict:
    odd_mask = (d.astype(int) % 2) == 1
    even_mask = ~odd_mask

    y_tr, d_tr = y[odd_mask], d[odd_mask]
    y_te, d_te = y[even_mask], d[even_mask]

    fit = fit_nonlinear_on_train(
        y_tr,
        d_tr,
        t,
        nodes,
        tau,
        lambda_c,
        lambda_b,
        lambda_amp,
        model_kind,
        rng,
        n_samples=200,
    )

    X_te = design_matrix(
        d_te,
        fit["omega"],
        fit["phi"],
        fit["beta"],
        fit["rho"],
        fit["xi"],
        fit["zeta"],
        fit["eta"],
        nodes,
        tau,
        model_kind,
    )

    # Refit linear coefficients on train only, then evaluate on holdout.
    X_tr = design_matrix(
        d_tr,
        fit["omega"],
        fit["phi"],
        fit["beta"],
        fit["rho"],
        fit["xi"],
        fit["zeta"],
        fit["eta"],
        nodes,
        tau,
        model_kind,
    )
    coef = fit_linear(X_tr, y_tr, model_kind, lambda_amp)
    yhat_te = X_te @ coef
    rmse_holdout = float(np.sqrt(np.mean((y_te - yhat_te) ** 2)))

    return {
        "rmse_holdout": rmse_holdout,
        "train_rmse": fit["train_rmse"],
        "canon_score": fit["canon_score"],
        "params": {k: float(fit[k]) for k in ["omega", "phi", "beta", "rho", "xi", "zeta", "eta"]},
    }


def summarize(rows: List[Dict]) -> Dict:
    hold = np.array([r["rmse_holdout"] for r in rows], dtype=float)
    train = np.array([r["train_rmse"] for r in rows], dtype=float)
    cs = np.array([r["canon_score"] for r in rows], dtype=float)
    return {
        "n": len(rows),
        "rmse_holdout_median": float(np.median(hold)),
        "rmse_train_median": float(np.median(train)),
        "generalization_gap_median": float(np.median(hold - train)),
        "canon_median": float(np.median(cs)),
    }


def main() -> None:
    d1739 = read_json("report_qw1739_signed_dynamic_micromodel_derivation.json")
    d1862 = read_json("report_qw1862_micromodel_canonical_compatibility_audit.json")
    d1894 = read_json("report_qw1894_amplitude_channel_reformulation.json")

    # Use strict test seeds from split 188000 (same as 1891/1894).
    all_rows = sorted(d1739.get("rows", []), key=lambda x: int(x["seed"]))

    def split_dataset(rows: List[Dict], seed: int = 188000) -> Dict[str, List[Dict]]:
        rng = np.random.default_rng(seed)
        by_n = {}
        for r in rows:
            by_n.setdefault(int(r["n_nodes"]), []).append(r)
        split = {"train": [], "val": [], "test": []}
        for n in sorted(by_n.keys()):
            g = list(by_n[n])
            idx = np.arange(len(g))
            rng.shuffle(idx)
            g = [g[i] for i in idx]
            n_train = min(8, len(g) - 4)
            n_val = min(3, len(g) - n_train - 1)
            split["train"].extend(g[:n_train])
            split["val"].extend(g[n_train : n_train + n_val])
            split["test"].extend(g[n_train + n_val :])
        return split

    split = split_dataset(all_rows, seed=188000)

    target = d1862.get("target_tuple", {})
    t = {
        "omega": float(target.get("omega", math.pi / 4.0)),
        "phi": float(target.get("phi", math.pi / 6.0)),
        "beta": float(target.get("beta", 0.01)),
    }

    nodes = [2, 5, 8, 11]
    tau = 2.5

    lc = 0.2
    lb = float(d1894.get("protocol", {}).get("selected_lambda_b", 0.35))
    la = float(d1894.get("protocol", {}).get("selected_lambda_amp", 0.01))

    rng = np.random.default_rng(189500)

    rows_base = []
    rows_amp = []

    for r in split["test"]:
        y = np.array([float(r["profile"][str(i)]) for i in range(1, 13)], dtype=float)
        d = np.arange(1, len(y) + 1, dtype=float)

        rows_base.append(
            odd_even_holdout_eval(
                y,
                d,
                t,
                nodes,
                tau,
                lambda_c=lc,
                lambda_b=lb,
                lambda_amp=0.0,
                model_kind="base",
                rng=rng,
            )
        )

        rows_amp.append(
            odd_even_holdout_eval(
                y,
                d,
                t,
                nodes,
                tau,
                lambda_c=lc,
                lambda_b=lb,
                lambda_amp=la,
                model_kind="amp",
                rng=rng,
            )
        )

    s_base = summarize(rows_base)
    s_amp = summarize(rows_amp)

    delta = {
        "holdout_rmse_gain": float(s_base["rmse_holdout_median"] - s_amp["rmse_holdout_median"]),
        "train_rmse_gain": float(s_base["rmse_train_median"] - s_amp["rmse_train_median"]),
        "gap_change": float(s_base["generalization_gap_median"] - s_amp["generalization_gap_median"]),
        "canon_gain": float(s_amp["canon_median"] - s_base["canon_median"]),
    }

    # Complexity check.
    n_obs = 12
    k_base = 9   # 7 nonlinear + 2 linear
    k_amp = 12   # 7 nonlinear + 5 linear
    complexity = {
        "n_obs": n_obs,
        "k_base": k_base,
        "k_amp": k_amp,
        "residual_dof_base": n_obs - k_base,
        "residual_dof_amp": n_obs - k_amp,
    }

    if complexity["residual_dof_amp"] <= 0:
        verdict = "AMPLITUDE_BRANCH_OVERPARAMETERIZED_NOT_ADMISSIBLE"
    elif delta["holdout_rmse_gain"] > 0 and delta["canon_gain"] >= 0:
        verdict = "JOINT_PARETO_GATE_PASS"
    elif delta["holdout_rmse_gain"] > 0:
        verdict = "JOINT_PARETO_GATE_PARTIAL"
    else:
        verdict = "JOINT_PARETO_GATE_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "evaluation": "odd-even intra-profile holdout",
            "split_seed": 188000,
            "lambda_c": lc,
            "lambda_b": lb,
            "lambda_amp": la,
        },
        "base_summary": s_base,
        "amp_summary": s_amp,
        "delta_amp_minus_base": delta,
        "complexity": complexity,
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1895: JOINT PARETO GATE UNDER CONSTRAINTS",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Holdout (odd-even) Summary",
        f"- base rmse_holdout median: {s_base['rmse_holdout_median']:.4f}",
        f"- amp  rmse_holdout median: {s_amp['rmse_holdout_median']:.4f}",
        f"- holdout rmse gain (amp): {delta['holdout_rmse_gain']:.4e}",
        f"- canon gain (amp): {delta['canon_gain']:.4e}",
        "",
        "## Complexity Check",
        f"- n_obs: {n_obs}",
        f"- k_base: {k_base}, residual dof: {complexity['residual_dof_base']}",
        f"- k_amp: {k_amp}, residual dof: {complexity['residual_dof_amp']}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1895] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1895] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
