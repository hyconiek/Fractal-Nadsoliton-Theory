#!/usr/bin/env python3
"""
QW-1941: Triple-sector shared model with AIC/BIC complexity control.

Compares:
- M0: strict derived model (no extra parameter),
- M1: one global shared scalar lambda (same lambda for mass, flavor, GW).

No per-sector parameter retuning is allowed in either model.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd
from scipy.linalg import polar


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1941_triple_sector_shared_complexity_aic_bic.json"
OUT_MD = ROOT / "RAPORT_QW1941_TRIPLE_SECTOR_SHARED_COMPLEXITY_AIC_BIC.md"


CKM_REF = np.array(
    [
        [0.97401, 0.22650, 0.00361],
        [0.22636, 0.97320, 0.04053],
        [0.00854, 0.03978, 0.99917],
    ],
    dtype=float,
)
PMNS_REF = np.array(
    [
        [0.821, 0.550, 0.150],
        [0.432, 0.582, 0.693],
        [0.378, 0.598, 0.707],
    ],
    dtype=float,
)

PARTICLES = [
    ("Top", 0.0, 173_000.0),
    ("Bottom", 7.0, 4_180.0),
    ("Tau", 9.0, 1_776.9),
    ("Charm", 9.0, 1_270.0),
    ("Muon", 14.0, 105.7),
    ("Electron", 24.0, 0.511),
]


def kernel_fn(d: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d + phi) / (1.0 + beta * (d**eta))


def cyclic_distance_matrix(q_left: np.ndarray, q_right: np.ndarray, modulus: int = 24) -> np.ndarray:
    dq = np.abs(q_left[:, None] - q_right[None, :]) % float(modulus)
    return np.minimum(dq, float(modulus) - dq)


def derive_shared_params(omega: float, phi: float, beta: float, eta: float) -> Dict[str, float]:
    d = np.arange(1.0, 13.0, dtype=float)
    k = np.abs(kernel_fn(d, omega, phi, beta, eta))
    w = k / max(np.sum(k), 1e-15)

    mean_d = float(np.sum(w * d))
    var_d = float(np.sum(w * (d - mean_d) ** 2))
    decay_ratio = float(k[3] / max(k[0], 1e-15))

    p_amp = float(np.clip(1.0 + 0.60 * np.tanh((mean_d - 2.0) / 2.0), 0.30, 2.20))
    r_dist = float(np.clip(np.tanh((var_d - 1.0) / 2.5), -1.20, 1.80))
    phase_scale = float(np.clip(1.0 + 0.70 * np.tanh(abs(phi)) + 0.25 * np.tanh(1.0 - decay_ratio), 0.00, 3.00))

    return {"p_amp": p_amp, "r_dist": r_dist, "phase_scale": phase_scale}


def flavor_prediction_abs(
    q_left: np.ndarray,
    q_right: np.ndarray,
    p_amp: float,
    r_dist: float,
    phase_scale: float,
    omega: float,
    phi: float,
    beta: float,
    eta: float,
) -> np.ndarray:
    n = len(q_left)
    d = 1.0 + cyclic_distance_matrix(q_left, q_right, modulus=24)
    k = kernel_fn(d, omega=omega, phi=phi, beta=beta, eta=eta)
    amp = np.sign(k) * ((np.abs(k) ** p_amp) * (d**r_dist))

    idx = np.arange(n, dtype=float)
    gap = np.abs(idx[:, None] - idx[None, :])
    sign = np.where((idx[:, None] - idx[None, :]) < 0.0, 1.0, -1.0)
    phase = np.exp(1j * (phi + phase_scale * omega * sign * gap))

    m = amp * phase
    u = polar(m)[0]
    return np.abs(u)


def flavor_metrics(params: Dict[str, float], kernel: Dict[str, float]) -> Dict[str, float]:
    q_up = np.array([0.0, 9.0, 14.0], dtype=float)
    q_down = np.array([7.0, 9.0, 14.0], dtype=float)
    q_nu = np.array([0.0, 1.0, 2.0], dtype=float)
    q_lep = np.array([24.0, 14.0, 9.0], dtype=float)

    ckm = flavor_prediction_abs(
        q_up,
        q_down,
        params["p_amp"],
        params["r_dist"],
        params["phase_scale"],
        kernel["omega"],
        kernel["phi"],
        kernel["beta"],
        kernel["eta"],
    )
    pmns = flavor_prediction_abs(
        q_nu,
        q_lep,
        params["p_amp"],
        params["r_dist"],
        params["phase_scale"],
        kernel["omega"],
        kernel["phi"],
        kernel["beta"],
        kernel["eta"],
    )

    ckm_rel = np.abs(ckm - CKM_REF) / np.clip(CKM_REF, 1e-12, None)
    pmns_rel = np.abs(pmns - PMNS_REF) / np.clip(PMNS_REF, 1e-12, None)
    return {
        "ckm_mean_rel_pct": float(100.0 * np.mean(ckm_rel)),
        "pmns_mean_rel_pct": float(100.0 * np.mean(pmns_rel)),
    }


def mass_metrics(gamma_scale: float, kernel: Dict[str, float]) -> Dict[str, float]:
    k1 = float(abs(kernel_fn(np.array([1.0]), kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])[0]))
    k4 = float(abs(kernel_fn(np.array([4.0]), kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])[0]))
    gamma_kernel = float(-4.0 * math.log(max(k4 / max(k1, 1e-15), 1e-15), 4.0) / 3.0)
    gamma_eff = float(np.clip(gamma_scale * gamma_kernel, 0.40, 2.60))

    errs = []
    for _, q, m_exp in PARTICLES:
        pred = 173_000.0 * (4.0 ** (-(gamma_eff * q / 4.0)))
        err = abs(pred - m_exp) / max(m_exp, 1e-15) * 100.0
        errs.append(err)
    return {
        "gamma_kernel": gamma_kernel,
        "gamma_eff": gamma_eff,
        "mean_rel_err_pct": float(np.mean(errs)),
        "max_rel_err_pct": float(np.max(errs)),
    }


def rank_auc_pos_gt_neg(pos: np.ndarray, neg: np.ndarray) -> float:
    y = np.concatenate([np.ones(len(pos), dtype=int), np.zeros(len(neg), dtype=int)])
    s = np.concatenate([pos, neg])
    n1 = len(pos)
    n0 = len(neg)
    order = np.argsort(s)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(s) + 1, dtype=float)
    rs = float(np.sum(ranks[y == 1]))
    return float((rs - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def gw_metrics(params: Dict[str, float], kernel: Dict[str, float], gw_cache: Dict[str, np.ndarray]) -> Dict[str, float]:
    d = np.array([1.0, 2.0, 3.0, 4.0], dtype=float)
    raw_w = (np.abs(kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])) ** params["p_amp"]) * (
        d**params["r_dist"]
    )
    w = raw_w / max(np.sum(raw_w), 1e-15)

    score = (
        w[0] * gw_cache["f_max"]
        + w[1] * gw_cache["f_mean"]
        + w[2] * gw_cache["f_0ms"]
        + w[3] * gw_cache["f_10ms"]
    )
    pair = gw_cache["pair"]

    s_hl = score[pair == "H1-L1"]
    s_hv = score[pair == "H1-V1"]
    s_lv = score[pair == "L1-V1"]
    s_ctrl = np.concatenate([s_hv, s_lv])

    q90 = float(np.quantile(s_ctrl, 0.90))
    auc = rank_auc_pos_gt_neg(s_hl, s_ctrl)
    adv = float(np.mean(s_hl > q90) - np.mean(s_ctrl > q90))
    sep = float(np.median(s_hl) - np.median(s_ctrl))
    ctrl_gap = float(abs(np.median(s_hv) - np.median(s_lv)))

    return {
        "auc_h1l1_vs_ctrl": auc,
        "adv_shared_minus_ctrl_q90": adv,
        "sep_median_h1l1_minus_ctrl": sep,
        "control_median_gap": ctrl_gap,
    }


def evaluate_lambda(
    lam: float,
    base_params: Dict[str, float],
    kernel: Dict[str, float],
    thresholds: Dict[str, float],
    gw_cache: Dict[str, np.ndarray],
) -> Dict[str, object]:
    params = {
        "p_amp": float(np.clip(base_params["p_amp"] * (lam**0.75), 0.20, 2.50)),
        "r_dist": float(np.clip(base_params["r_dist"] * lam, -1.50, 2.00)),
        "phase_scale": float(np.clip(1.0 + (base_params["phase_scale"] - 1.0) * lam, 0.00, 3.00)),
        "gamma_scale": float(np.clip(lam, 0.20, 1.80)),
    }

    mass = mass_metrics(params["gamma_scale"], kernel)
    flavor = flavor_metrics(params, kernel)
    gw = gw_metrics(params, kernel, gw_cache)

    flags = {
        "mass_mean_rel_pct_le_max": bool(mass["mean_rel_err_pct"] <= thresholds["mass_mean_rel_pct_max"]),
        "mass_max_rel_pct_le_max": bool(mass["max_rel_err_pct"] <= thresholds["mass_max_rel_pct_max"]),
        "ckm_mean_rel_pct_le_max": bool(flavor["ckm_mean_rel_pct"] <= thresholds["ckm_mean_rel_pct_max"]),
        "pmns_mean_rel_pct_le_max": bool(flavor["pmns_mean_rel_pct"] <= thresholds["pmns_mean_rel_pct_max"]),
        "gw_sep_ge_min": bool(gw["sep_median_h1l1_minus_ctrl"] >= thresholds["gw_sep_min"]),
        "gw_adv_ge_min": bool(gw["adv_shared_minus_ctrl_q90"] >= thresholds["gw_adv_min"]),
        "gw_auc_ge_min": bool(gw["auc_h1l1_vs_ctrl"] >= thresholds["gw_auc_min"]),
        "gw_control_gap_le_max": bool(gw["control_median_gap"] <= thresholds["gw_control_gap_max"]),
    }

    mass_loss = mass["mean_rel_err_pct"] / thresholds["mass_mean_rel_pct_max"] + mass["max_rel_err_pct"] / thresholds["mass_max_rel_pct_max"]
    flavor_loss = flavor["ckm_mean_rel_pct"] / thresholds["ckm_mean_rel_pct_max"] + flavor["pmns_mean_rel_pct"] / thresholds["pmns_mean_rel_pct_max"]
    gw_loss = (
        max(0.0, thresholds["gw_sep_min"] - gw["sep_median_h1l1_minus_ctrl"]) / max(thresholds["gw_sep_min"], 1e-12)
        + max(0.0, thresholds["gw_adv_min"] - gw["adv_shared_minus_ctrl_q90"]) / max(thresholds["gw_adv_min"], 1e-12)
        + max(0.0, thresholds["gw_auc_min"] - gw["auc_h1l1_vs_ctrl"]) / max(thresholds["gw_auc_min"], 1e-12)
        + max(0.0, gw["control_median_gap"] - thresholds["gw_control_gap_max"]) / max(thresholds["gw_control_gap_max"], 1e-12)
    )
    total_loss = float(mass_loss + flavor_loss + gw_loss)

    return {
        "lambda": float(lam),
        "params": params,
        "mass": mass,
        "flavor": flavor,
        "gw": gw,
        "flags": flags,
        "all_pass": bool(all(flags.values())),
        "loss_components": {
            "mass_loss": float(mass_loss),
            "flavor_loss": float(flavor_loss),
            "gw_loss": float(gw_loss),
        },
        "total_loss": total_loss,
    }


def info_criteria(loss: float, n_obs: int, k_params: int) -> Tuple[float, float]:
    # Proxy criterion based on normalized loss.
    rss_proxy = max((loss**2) * n_obs, 1e-12)
    aic = float(n_obs * math.log(rss_proxy / n_obs) + 2 * k_params)
    bic = float(n_obs * math.log(rss_proxy / n_obs) + k_params * math.log(n_obs))
    return aic, bic


def main() -> None:
    d1932 = json.loads((ROOT / "report_qw1932_physical_reparameterization_eta_scan.json").read_text(encoding="utf-8"))
    sel = d1932["selected"]
    kernel = {
        "omega": float(sel["fit"]["omega"]),
        "phi": float(sel["fit"]["phi"]),
        "beta": float(sel["fit"]["beta"]),
        "eta": float(sel["eta"]),
    }

    thresholds = {
        "mass_mean_rel_pct_max": 15.0,
        "mass_max_rel_pct_max": 75.0,
        "ckm_mean_rel_pct_max": 15.0,
        "pmns_mean_rel_pct_max": 15.0,
        "gw_sep_min": 0.0020,
        "gw_adv_min": 0.30,
        "gw_auc_min": 0.75,
        "gw_control_gap_max": 0.0025,
    }
    base_params = derive_shared_params(**kernel)

    df_gw = pd.read_csv(ROOT / "gw1831_window_features.csv")
    req = ["pair", "max_abs_corr", "mean_abs_corr", "corr_at_0ms", "corr_at_10ms"]
    miss = [c for c in req if c not in df_gw.columns]
    if miss:
        raise RuntimeError(f"Missing GW columns: {miss}")
    gw_cache = {
        "pair": df_gw["pair"].astype(str).to_numpy(),
        "f_max": df_gw["max_abs_corr"].to_numpy(dtype=float),
        "f_mean": df_gw["mean_abs_corr"].to_numpy(dtype=float),
        "f_0ms": df_gw["corr_at_0ms"].to_numpy(dtype=float),
        "f_10ms": df_gw["corr_at_10ms"].to_numpy(dtype=float),
    }

    n_obs = 28  # 6 mass + 9 CKM + 9 PMNS + 4 GW summary constraints

    m0 = evaluate_lambda(1.0, base_params, kernel, thresholds, gw_cache)
    m0_aic, m0_bic = info_criteria(m0["total_loss"], n_obs=n_obs, k_params=0)

    lambdas = np.linspace(0.50, 1.50, 101)
    best = None
    for lam in lambdas:
        cand = evaluate_lambda(float(lam), base_params, kernel, thresholds, gw_cache)
        if best is None or cand["total_loss"] < best["total_loss"]:
            best = cand

    m1_aic, m1_bic = info_criteria(best["total_loss"], n_obs=n_obs, k_params=1)
    delta_aic_m1_minus_m0 = float(m1_aic - m0_aic)
    delta_bic_m1_minus_m0 = float(m1_bic - m0_bic)

    if m0["all_pass"]:
        verdict = "TRIPLE_SECTOR_STRICT_DERIVED_PASS"
        chosen_model = "M0_STRICT_DERIVED"
        required_next = "LOCK_DERIVED_MODEL_FOR_FINAL_GATE"
    elif best["all_pass"] and delta_bic_m1_minus_m0 <= -6.0:
        verdict = "TRIPLE_SECTOR_ONE_PARAM_EXTENSION_SUPPORTED"
        chosen_model = "M1_ONE_SHARED_LAMBDA"
        required_next = "FREEZE_SHARED_LAMBDA_AND_RUN_FINAL_GATE"
    elif best["all_pass"]:
        verdict = "TRIPLE_SECTOR_PASS_BUT_COMPLEXITY_NOT_JUSTIFIED"
        chosen_model = "M1_ONE_SHARED_LAMBDA"
        required_next = "DO_NOT_PROMOTE_M1_UNTIL_COMPLEXITY_PENALTY_IS_RESOLVED"
    else:
        verdict = "TRIPLE_SECTOR_NOT_CLOSED_UNDER_SHARED_COMPLEXITY_CONTROL"
        chosen_model = "M0_STRICT_DERIVED" if m0["total_loss"] <= best["total_loss"] else "M1_ONE_SHARED_LAMBDA"
        required_next = "REDESIGN_SHARED_FLAVOR_MASS_LINK_UNDER_SINGLE_KERNEL_CONSTRAINT"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw1932_physical_reparameterization_eta_scan.json:selected",
        "kernel": kernel,
        "base_params": base_params,
        "thresholds": thresholds,
        "models": {
            "M0_strict_derived": {
                **m0,
                "aic_proxy": m0_aic,
                "bic_proxy": m0_bic,
                "k_params": 0,
            },
            "M1_one_shared_lambda_best": {
                **best,
                "aic_proxy": m1_aic,
                "bic_proxy": m1_bic,
                "k_params": 1,
            },
        },
        "model_comparison": {
            "delta_aic_m1_minus_m0": delta_aic_m1_minus_m0,
            "delta_bic_m1_minus_m0": delta_bic_m1_minus_m0,
        },
        "chosen_model": chosen_model,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1941: TRIPLE-SECTOR SHARED COMPLEXITY (AIC/BIC)",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Chosen model: **{chosen_model}**",
        "",
        "## M0 Strict Derived",
        f"- total_loss: {m0['total_loss']:.4f}",
        f"- all_pass: {m0['all_pass']}",
        f"- AIC/BIC proxy: {m0_aic:.4f}/{m0_bic:.4f}",
        "",
        "## M1 One Shared Lambda (best)",
        f"- lambda: {best['lambda']:.4f}",
        f"- total_loss: {best['total_loss']:.4f}",
        f"- all_pass: {best['all_pass']}",
        f"- AIC/BIC proxy: {m1_aic:.4f}/{m1_bic:.4f}",
        "",
        "## Comparison",
        f"- delta_aic (M1-M0): {delta_aic_m1_minus_m0:.4f}",
        f"- delta_bic (M1-M0): {delta_bic_m1_minus_m0:.4f}",
        "",
        "## Required Next Step",
        f"- {required_next}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1941] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1941] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1941] verdict={verdict} chosen={chosen_model}")


if __name__ == "__main__":
    main()
