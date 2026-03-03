#!/usr/bin/env python3
"""
QW-1962: Unified triad gate for best noncircular mass branch from QW-1961.

Strict protocol:
- one frozen kernel,
- one deterministic shared flavor/GW operator map (no fitting),
- noncircular mass chain from QW-1961 best_noncircular.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
from scipy.linalg import polar


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1962_noncircular_branch_unified_triad_gate.json"
OUT_MD = ROOT / "RAPORT_QW1962_NONCIRCULAR_BRANCH_UNIFIED_TRIAD_GATE.md"


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


MASS_EXP = {
    "Top": 173_000.0,
    "Bottom": 4_180.0,
    "Tau": 1_776.9,
    "Charm": 1_270.0,
    "Muon": 105.7,
    "Electron": 0.511,
}


def kernel_fn(d: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d + phi) / (1.0 + beta * (d**eta))


def cyclic_distance_matrix(q_left: np.ndarray, q_right: np.ndarray, modulus: int = 24) -> np.ndarray:
    dq = np.abs(q_left[:, None] - q_right[None, :]) % float(modulus)
    return np.minimum(dq, float(modulus) - dq)


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


def derive_shared_params(omega: float, phi: float, beta: float, eta: float) -> Dict[str, float]:
    # Deterministic map reused from QW-1943 family (no fitting).
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


def flavor_metrics(
    q_up: np.ndarray,
    q_down: np.ndarray,
    q_nu: np.ndarray,
    q_lep: np.ndarray,
    params: Dict[str, float],
    kernel: Dict[str, float],
) -> Dict[str, float]:
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


def mass_metrics(
    q_assign: Dict[str, float],
    gamma: float,
    delta_info: float,
) -> Dict[str, object]:
    sector = {
        "Top": 0.0,
        "Bottom": +1.0,
        "Charm": +1.0,
        "Tau": -1.0,
        "Muon": -1.0,
        "Electron": -1.0,
    }
    rows = []
    errs = []
    pred_tau = None
    pred_charm = None
    for name, exp_mev in MASS_EXP.items():
        q_base = float(q_assign[name])
        q_eff = q_base + sector[name] * delta_info
        pred = 173_000.0 * (4.0 ** (-(gamma * q_eff / 4.0)))
        err = abs(pred - exp_mev) / max(exp_mev, 1e-15) * 100.0
        rows.append(
            {
                "particle": name,
                "q_base": q_base,
                "q_eff": float(q_eff),
                "exp_mev": float(exp_mev),
                "pred_mev": float(pred),
                "rel_err_pct": float(err),
            }
        )
        errs.append(err)
        if name == "Tau":
            pred_tau = pred
        if name == "Charm":
            pred_charm = pred

    tau_charm_ratio_pred = float(pred_tau / max(pred_charm, 1e-15))
    tau_charm_ratio_exp = float(MASS_EXP["Tau"] / MASS_EXP["Charm"])
    tau_charm_ratio_err = abs(tau_charm_ratio_pred - tau_charm_ratio_exp) / tau_charm_ratio_exp * 100.0

    return {
        "rows": rows,
        "mean_rel_err_pct": float(np.mean(errs)),
        "max_rel_err_pct": float(np.max(errs)),
        "tau_charm_ratio_pred": tau_charm_ratio_pred,
        "tau_charm_ratio_exp": tau_charm_ratio_exp,
        "tau_charm_ratio_rel_err_pct": float(tau_charm_ratio_err),
    }


def gw_metrics(params: Dict[str, float], kernel: Dict[str, float]) -> Dict[str, float]:
    d = np.array([1.0, 2.0, 3.0, 4.0], dtype=float)
    raw_w = (np.abs(kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])) ** params["p_amp"]) * (
        d ** params["r_dist"]
    )
    w = raw_w / np.sum(raw_w)

    df_gw = pd.read_csv(ROOT / "gw1831_window_features.csv")
    score = (
        w[0] * df_gw["max_abs_corr"].to_numpy(dtype=float)
        + w[1] * df_gw["mean_abs_corr"].to_numpy(dtype=float)
        + w[2] * df_gw["corr_at_0ms"].to_numpy(dtype=float)
        + w[3] * df_gw["corr_at_10ms"].to_numpy(dtype=float)
    )
    pair = df_gw["pair"].astype(str).to_numpy()
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


def main() -> None:
    r1961 = json.loads((ROOT / "report_qw1961_noncircular_gamma_q_derivation_matrix.json").read_text(encoding="utf-8"))
    best = r1961["summary"]["best_noncircular"]
    q_assign_name = best["q_assignment"]
    q_assign = r1961["inputs"]["q_assignments"][q_assign_name]
    gamma = float(best["gamma_value"])
    delta_info = float(r1961["info_split_source"]["delta_info"])

    d1932 = json.loads((ROOT / "report_qw1932_physical_reparameterization_eta_scan.json").read_text(encoding="utf-8"))
    sel = d1932["selected"]
    kernel = {
        "omega": float(sel["fit"]["omega"]),
        "phi": float(sel["fit"]["phi"]),
        "beta": float(sel["fit"]["beta"]),
        "eta": float(sel["eta"]),
    }
    params = derive_shared_params(**kernel)

    q_up = np.array([q_assign["Top"], q_assign["Charm"], q_assign["Muon"]], dtype=float)
    q_down = np.array([q_assign["Bottom"], q_assign["Tau"], q_assign["Muon"]], dtype=float)
    q_nu = np.array([0.0, 1.0, 2.0], dtype=float)
    q_lep = np.array([q_assign["Electron"], q_assign["Muon"], q_assign["Tau"]], dtype=float)

    m = mass_metrics(q_assign=q_assign, gamma=gamma, delta_info=delta_info)
    f = flavor_metrics(q_up=q_up, q_down=q_down, q_nu=q_nu, q_lep=q_lep, params=params, kernel=kernel)
    g = gw_metrics(params=params, kernel=kernel)

    thresholds = {
        "mass_mean_rel_pct_max": 15.0,
        "mass_max_rel_pct_max": 35.0,
        "tau_charm_ratio_rel_err_pct_max": 20.0,
        "ckm_mean_rel_pct_max": 15.0,
        "pmns_mean_rel_pct_max": 15.0,
        "gw_sep_min": 0.0020,
        "gw_adv_min": 0.30,
        "gw_auc_min": 0.75,
        "gw_control_gap_max": 0.0025,
    }

    flags = {
        "mass_mean_rel_pct_le_max": bool(m["mean_rel_err_pct"] <= thresholds["mass_mean_rel_pct_max"]),
        "mass_max_rel_pct_le_max": bool(m["max_rel_err_pct"] <= thresholds["mass_max_rel_pct_max"]),
        "mass_tau_charm_ratio_err_le_max": bool(m["tau_charm_ratio_rel_err_pct"] <= thresholds["tau_charm_ratio_rel_err_pct_max"]),
        "ckm_mean_rel_pct_le_max": bool(f["ckm_mean_rel_pct"] <= thresholds["ckm_mean_rel_pct_max"]),
        "pmns_mean_rel_pct_le_max": bool(f["pmns_mean_rel_pct"] <= thresholds["pmns_mean_rel_pct_max"]),
        "gw_sep_ge_min": bool(g["sep_median_h1l1_minus_ctrl"] >= thresholds["gw_sep_min"]),
        "gw_adv_ge_min": bool(g["adv_shared_minus_ctrl_q90"] >= thresholds["gw_adv_min"]),
        "gw_auc_ge_min": bool(g["auc_h1l1_vs_ctrl"] >= thresholds["gw_auc_min"]),
        "gw_control_gap_le_max": bool(g["control_median_gap"] <= thresholds["gw_control_gap_max"]),
    }
    all_pass = bool(all(flags.values()))

    verdict = "NONCIRCULAR_UNIFIED_TRIAD_PASS" if all_pass else "NONCIRCULAR_UNIFIED_TRIAD_FAIL"
    required_next = (
        "LOCK_NONCIRCULAR_BRANCH_AND_PREPARE_EXTERNAL_CONFIRMATORY_PACKAGE"
        if all_pass
        else "KEEP_NONCIRCULAR_MASS_BRANCH_AND_REWORK_SHARED_FLAVOR_GW_OPERATOR"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_qw1961_best_noncircular": best,
        "kernel_source": "report_qw1932_physical_reparameterization_eta_scan.json:selected",
        "kernel": kernel,
        "deterministic_shared_params": params,
        "q_assignment_used": q_assign,
        "gamma_used": gamma,
        "delta_info_used": delta_info,
        "metrics": {
            "mass": m,
            "flavor": f,
            "gw": g,
        },
        "thresholds": thresholds,
        "flags": flags,
        "all_pass": all_pass,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1962: NONCIRCULAR BRANCH UNIFIED TRIAD GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Branch Input",
        f"- q_assignment: {q_assign_name}",
        f"- gamma_used: {gamma:.6f}",
        f"- delta_info_used: {delta_info:.6f}",
        (
            f"- shared params p/r/phase: "
            f"{params['p_amp']:.4f}/{params['r_dist']:.4f}/{params['phase_scale']:.4f}"
        ),
        "",
        "## Mass",
        (
            f"- mean/max rel err: {m['mean_rel_err_pct']:.3f}% / {m['max_rel_err_pct']:.3f}%"
        ),
        (
            f"- tau/charm ratio pred/exp/error: "
            f"{m['tau_charm_ratio_pred']:.4f}/{m['tau_charm_ratio_exp']:.4f}/"
            f"{m['tau_charm_ratio_rel_err_pct']:.3f}%"
        ),
        "",
        "## Flavor",
        f"- CKM mean rel err: {f['ckm_mean_rel_pct']:.3f}%",
        f"- PMNS mean rel err: {f['pmns_mean_rel_pct']:.3f}%",
        "",
        "## GW",
        (
            f"- auc/adv/sep/gap: "
            f"{g['auc_h1l1_vs_ctrl']:.4f}/"
            f"{g['adv_shared_minus_ctrl_q90']:.4f}/"
            f"{g['sep_median_h1l1_minus_ctrl']:.6f}/"
            f"{g['control_median_gap']:.6f}"
        ),
        "",
        "## Flags",
    ]
    for k, v in flags.items():
        lines.append(f"- {k}: {v}")
    lines.extend(
        [
            "",
            "## Required Next Step",
            f"- {required_next}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1962] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1962] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1962] verdict={verdict}")


if __name__ == "__main__":
    main()
