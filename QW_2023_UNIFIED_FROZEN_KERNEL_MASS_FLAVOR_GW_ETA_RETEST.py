#!/usr/bin/env python3
"""
QW-2023: Unified frozen-kernel test for mass + flavor + GW without sector retuning.

Core requirement:
- one frozen kernel tuple (omega, phi, beta, eta),
- one shared operator family across all sectors,
- no separate per-sector parameter resets.

This script evaluates:
1) a strict derived candidate (no fitting),
2) a single global fit candidate (one shared parameter vector for all sectors).
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
OUT_JSON = ROOT / "report_qw2023_unified_frozen_kernel_mass_flavor_gw_eta_retest.json"
OUT_MD = ROOT / "RAPORT_QW2023_UNIFIED_FROZEN_KERNEL_MASS_FLAVOR_GW_ETA_RETEST.md"


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
    gen_gap = np.abs(idx[:, None] - idx[None, :])
    sign = np.where((idx[:, None] - idx[None, :]) < 0.0, 1.0, -1.0)
    phase = np.exp(1j * (phi + phase_scale * omega * sign * gen_gap))

    m = amp * phase
    u = polar(m)[0]
    return np.abs(u)


def flavor_metrics(
    p_amp: float,
    r_dist: float,
    phase_scale: float,
    omega: float,
    phi: float,
    beta: float,
    eta: float,
) -> Dict[str, object]:
    q_up = np.array([0.0, 9.0, 14.0], dtype=float)
    q_down = np.array([7.0, 9.0, 14.0], dtype=float)
    q_nu = np.array([0.0, 1.0, 2.0], dtype=float)
    q_lep = np.array([24.0, 14.0, 9.0], dtype=float)

    ckm = flavor_prediction_abs(q_up, q_down, p_amp, r_dist, phase_scale, omega, phi, beta, eta)
    pmns = flavor_prediction_abs(q_nu, q_lep, p_amp, r_dist, phase_scale, omega, phi, beta, eta)

    ckm_rel = np.abs(ckm - CKM_REF) / np.clip(CKM_REF, 1e-12, None)
    pmns_rel = np.abs(pmns - PMNS_REF) / np.clip(PMNS_REF, 1e-12, None)

    return {
        "ckm_pred_abs": ckm,
        "pmns_pred_abs": pmns,
        "ckm_mean_rel_pct": float(100.0 * np.mean(ckm_rel)),
        "ckm_max_rel_pct": float(100.0 * np.max(ckm_rel)),
        "pmns_mean_rel_pct": float(100.0 * np.mean(pmns_rel)),
        "pmns_max_rel_pct": float(100.0 * np.max(pmns_rel)),
    }


def mass_metrics(
    p_amp: float,
    r_dist: float,
    phase_scale: float,
    gamma_scale: float,
    omega: float,
    phi: float,
    beta: float,
    eta: float,
) -> Dict[str, object]:
    # Particle table aligned with prior reports.
    particles = [
        ("Top", 0.0, 173_000.0, 0.0, 3.0),
        ("Bottom", 7.0, 4_180.0, 1.0, 3.0),
        ("Tau", 9.0, 1_776.9, -1.0, 3.0),
        ("Charm", 9.0, 1_270.0, 1.0, 2.0),
        ("Muon", 14.0, 105.7, -1.0, 2.0),
        ("Electron", 24.0, 0.511, -1.0, 1.0),
    ]

    d1, d4 = np.array([1.0, 4.0], dtype=float)
    k1 = float(abs(kernel_fn(d1, omega, phi, beta, eta)))
    k4 = float(abs(kernel_fn(d4, omega, phi, beta, eta)))
    gamma_kernel = float(-4.0 * math.log(max(k4 / max(k1, 1e-15), 1e-15), 4.0) / 3.0)
    gamma_eff = float(np.clip(gamma_scale * gamma_kernel, 0.40, 2.60))

    # Shared correction coefficients, deterministic from same shared parameters.
    c_q = float(np.tanh(0.90 * (p_amp - 1.0)))
    c_s = float(np.tanh(0.70 * r_dist))
    c_g = float(np.tanh(0.85 * (phase_scale - 1.0)))

    rows = []
    errs = []
    for name, q, m_exp, sector, gen in particles:
        base = 173_000.0 * (4.0 ** (-(gamma_eff * q / 4.0)))
        delta = c_q * (q / 24.0) + c_s * sector + c_g * (gen - 2.0)
        pred = float(base * math.exp(delta))
        err = float(abs(pred - m_exp) / max(m_exp, 1e-15) * 100.0)
        errs.append(err)
        rows.append(
            {
                "particle": name,
                "q": q,
                "sector": sector,
                "gen": gen,
                "exp_mev": m_exp,
                "pred_mev": pred,
                "rel_err_pct": err,
            }
        )

    return {
        "gamma_kernel": gamma_kernel,
        "gamma_eff": gamma_eff,
        "coefficients": {"c_q": c_q, "c_s": c_s, "c_g": c_g},
        "rows": rows,
        "mean_rel_err_pct": float(np.mean(errs)),
        "max_rel_err_pct": float(np.max(errs)),
    }


def gw_metrics(
    p_amp: float,
    r_dist: float,
    omega: float,
    phi: float,
    beta: float,
    eta: float,
    gw_cache: Dict[str, np.ndarray],
) -> Dict[str, object]:
    d = np.array([1.0, 2.0, 3.0, 4.0], dtype=float)
    raw_w = (np.abs(kernel_fn(d, omega, phi, beta, eta)) ** p_amp) * (d**r_dist)
    w = raw_w / np.sum(raw_w)

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
        "weights": {
            "w_max_abs_corr": float(w[0]),
            "w_mean_abs_corr": float(w[1]),
            "w_corr_at_0ms": float(w[2]),
            "w_corr_at_10ms": float(w[3]),
        },
        "summary": {
            "auc_h1l1_vs_ctrl": float(auc),
            "adv_shared_minus_ctrl_q90": float(adv),
            "sep_median_h1l1_minus_ctrl": float(sep),
            "control_median_gap": float(ctrl_gap),
        },
    }


def evaluate_candidate(
    p_amp: float,
    r_dist: float,
    phase_scale: float,
    gamma_scale: float,
    omega: float,
    phi: float,
    beta: float,
    eta: float,
    thresholds: Dict[str, float],
    gw_cache: Dict[str, np.ndarray],
) -> Dict[str, object]:
    fm = flavor_metrics(p_amp, r_dist, phase_scale, omega, phi, beta, eta)
    mm = mass_metrics(p_amp, r_dist, phase_scale, gamma_scale, omega, phi, beta, eta)
    gm = gw_metrics(p_amp, r_dist, omega, phi, beta, eta, gw_cache=gw_cache)

    flags = {
        "mass_mean_rel_pct_le_max": bool(mm["mean_rel_err_pct"] <= thresholds["mass_mean_rel_pct_max"]),
        "mass_max_rel_pct_le_max": bool(mm["max_rel_err_pct"] <= thresholds["mass_max_rel_pct_max"]),
        "ckm_mean_rel_pct_le_max": bool(fm["ckm_mean_rel_pct"] <= thresholds["ckm_mean_rel_pct_max"]),
        "pmns_mean_rel_pct_le_max": bool(fm["pmns_mean_rel_pct"] <= thresholds["pmns_mean_rel_pct_max"]),
        "gw_sep_ge_min": bool(gm["summary"]["sep_median_h1l1_minus_ctrl"] >= thresholds["gw_sep_min"]),
        "gw_adv_ge_min": bool(gm["summary"]["adv_shared_minus_ctrl_q90"] >= thresholds["gw_adv_min"]),
        "gw_auc_ge_min": bool(gm["summary"]["auc_h1l1_vs_ctrl"] >= thresholds["gw_auc_min"]),
        "gw_control_gap_le_max": bool(gm["summary"]["control_median_gap"] <= thresholds["gw_control_gap_max"]),
    }

    # Joint objective used only for one-time global shared fit.
    gw_pen = (
        30.0 * max(0.0, thresholds["gw_auc_min"] - gm["summary"]["auc_h1l1_vs_ctrl"])
        + 40.0 * max(0.0, thresholds["gw_adv_min"] - gm["summary"]["adv_shared_minus_ctrl_q90"])
        + 250.0 * max(0.0, thresholds["gw_sep_min"] - gm["summary"]["sep_median_h1l1_minus_ctrl"])
        + 120.0 * max(0.0, gm["summary"]["control_median_gap"] - thresholds["gw_control_gap_max"])
    )
    objective = (
        fm["ckm_mean_rel_pct"]
        + fm["pmns_mean_rel_pct"]
        + 0.06 * mm["mean_rel_err_pct"]
        + 0.01 * mm["max_rel_err_pct"]
        + gw_pen
    )

    return {
        "params": {
            "p_amp": float(p_amp),
            "r_dist": float(r_dist),
            "phase_scale": float(phase_scale),
            "gamma_scale": float(gamma_scale),
        },
        "mass": mm,
        "flavor": {
            "ckm_mean_rel_pct": fm["ckm_mean_rel_pct"],
            "ckm_max_rel_pct": fm["ckm_max_rel_pct"],
            "pmns_mean_rel_pct": fm["pmns_mean_rel_pct"],
            "pmns_max_rel_pct": fm["pmns_max_rel_pct"],
            "ckm_pred_abs": fm["ckm_pred_abs"].tolist(),
            "pmns_pred_abs": fm["pmns_pred_abs"].tolist(),
        },
        "gw": gm,
        "flags": flags,
        "all_pass": bool(all(flags.values())),
        "objective": float(objective),
    }


def main() -> None:
    d2021 = json.loads((ROOT / "report_qw2021_v2_eta_operator_beta_constraint_scan.json").read_text(encoding="utf-8"))
    sel = d2021["selected"]["fit"]
    omega = float(sel["omega"])
    phi = float(sel["phi"])
    beta = float(sel["beta"])
    eta = float(sel["eta"])

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

    # Strict derived (no fitting): deterministic shared parameters.
    derived = evaluate_candidate(
        p_amp=1.0,
        r_dist=0.0,
        phase_scale=1.0,
        gamma_scale=1.0,
        omega=omega,
        phi=phi,
        beta=beta,
        eta=eta,
        thresholds=thresholds,
        gw_cache=gw_cache,
    )

    # Single global fit (shared across all sectors).
    p_grid = np.linspace(0.45, 1.55, 11)
    r_grid = np.linspace(-0.60, 0.90, 11)
    ph_grid = np.linspace(0.30, 2.00, 11)
    gs_grid = np.linspace(0.45, 1.05, 11)

    best = None
    feasible_count = 0
    eval_count = 0

    for p_amp in p_grid:
        for r_dist in r_grid:
            for phase_scale in ph_grid:
                for gamma_scale in gs_grid:
                    cand = evaluate_candidate(
                        p_amp=float(p_amp),
                        r_dist=float(r_dist),
                        phase_scale=float(phase_scale),
                        gamma_scale=float(gamma_scale),
                        omega=omega,
                        phi=phi,
                        beta=beta,
                        eta=eta,
                        thresholds=thresholds,
                        gw_cache=gw_cache,
                    )
                    eval_count += 1
                    if cand["all_pass"]:
                        feasible_count += 1
                    if best is None or cand["objective"] < best["objective"]:
                        best = cand

    # Deterministic local refinement around coarse best.
    bp = best["params"]
    p_ref = np.linspace(bp["p_amp"] - 0.10, bp["p_amp"] + 0.10, 7)
    r_ref = np.linspace(bp["r_dist"] - 0.12, bp["r_dist"] + 0.12, 7)
    ph_ref = np.linspace(bp["phase_scale"] - 0.15, bp["phase_scale"] + 0.15, 7)
    gs_ref = np.linspace(bp["gamma_scale"] - 0.10, bp["gamma_scale"] + 0.10, 7)

    for p_amp in p_ref:
        for r_dist in r_ref:
            for phase_scale in ph_ref:
                for gamma_scale in gs_ref:
                    cand = evaluate_candidate(
                        p_amp=float(np.clip(p_amp, 0.20, 2.50)),
                        r_dist=float(np.clip(r_dist, -1.50, 2.00)),
                        phase_scale=float(np.clip(phase_scale, 0.00, 3.00)),
                        gamma_scale=float(np.clip(gamma_scale, 0.20, 1.50)),
                        omega=omega,
                        phi=phi,
                        beta=beta,
                        eta=eta,
                        thresholds=thresholds,
                        gw_cache=gw_cache,
                    )
                    eval_count += 1
                    if cand["all_pass"]:
                        feasible_count += 1
                    if cand["objective"] < best["objective"]:
                        best = cand

    if derived["all_pass"]:
        verdict = "UNIFIED_FROZEN_KERNEL_DERIVED_PASS"
        required_next = "LOCK_AND_EXTERNALIZE_CONFIRMATORY_PROTOCOL_WITHOUT_PARAMETER_UPDATES"
    elif best["all_pass"]:
        verdict = "UNIFIED_FROZEN_KERNEL_PASS_WITH_SINGLE_GLOBAL_FIT"
        required_next = "FREEZE_GLOBAL_VECTOR_AND_RUN_INDEPENDENT_CONFIRMATORY_REPLICATION"
    else:
        verdict = "UNIFIED_FROZEN_KERNEL_NOT_CLOSED_TRIPLE_SECTOR"
        required_next = "REWORK_SHARED_OPERATOR_FOR_FLAVOR_WHILE_KEEPING_SINGLE_KERNEL_CONSTRAINT"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw2021_v2_eta_operator_beta_constraint_scan.json:selected.fit",
        "kernel": {"omega": omega, "phi": phi, "beta": beta, "eta": eta},
        "protocol": {
            "single_shared_operator": True,
            "sector_specific_parameter_reset": False,
            "coarse_grid_sizes": {
                "p_amp": int(len(p_grid)),
                "r_dist": int(len(r_grid)),
                "phase_scale": int(len(ph_grid)),
                "gamma_scale": int(len(gs_grid)),
            },
            "total_evaluations": int(eval_count),
            "feasible_count_all_pass": int(feasible_count),
        },
        "thresholds": thresholds,
        "derived_candidate": derived,
        "joint_single_fit_best": best,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    d_mass = derived["mass"]
    d_fl = derived["flavor"]
    d_gw = derived["gw"]["summary"]
    b_mass = best["mass"]
    b_fl = best["flavor"]
    b_gw = best["gw"]["summary"]

    lines = [
        "# RAPORT QW-2023: UNIFIED FROZEN KERNEL MASS+FLAVOR+GW (NO SECTOR RETUNE)",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Kernel: omega={omega:.6f}, phi={phi:.6f}, beta={beta:.6f}, eta={eta:.2f}",
        f"- Verdict: **{verdict}**",
        f"- Total evaluations: {eval_count}",
        f"- Feasible (all-pass) count: {feasible_count}",
        "",
        "## Derived Candidate (No Fit)",
        (
            f"- params p/r/phase/gamma_scale: "
            f"{derived['params']['p_amp']:.3f}/{derived['params']['r_dist']:.3f}/"
            f"{derived['params']['phase_scale']:.3f}/{derived['params']['gamma_scale']:.3f}"
        ),
        (
            f"- mass mean/max rel%: {d_mass['mean_rel_err_pct']:.3f}/{d_mass['max_rel_err_pct']:.3f} | "
            f"gamma_eff={d_mass['gamma_eff']:.4f}"
        ),
        (
            f"- flavor CKM/PMNS mean rel%: "
            f"{d_fl['ckm_mean_rel_pct']:.3f}/{d_fl['pmns_mean_rel_pct']:.3f}"
        ),
        (
            f"- GW auc/adv/sep/gap: "
            f"{d_gw['auc_h1l1_vs_ctrl']:.4f}/{d_gw['adv_shared_minus_ctrl_q90']:.4f}/"
            f"{d_gw['sep_median_h1l1_minus_ctrl']:.6f}/{d_gw['control_median_gap']:.6f}"
        ),
        f"- all_pass: {derived['all_pass']}",
        "",
        "## Best Single Global Fit (Shared Across All Sectors)",
        (
            f"- params p/r/phase/gamma_scale: "
            f"{best['params']['p_amp']:.3f}/{best['params']['r_dist']:.3f}/"
            f"{best['params']['phase_scale']:.3f}/{best['params']['gamma_scale']:.3f}"
        ),
        (
            f"- mass mean/max rel%: {b_mass['mean_rel_err_pct']:.3f}/{b_mass['max_rel_err_pct']:.3f} | "
            f"gamma_eff={b_mass['gamma_eff']:.4f}"
        ),
        (
            f"- flavor CKM/PMNS mean rel%: "
            f"{b_fl['ckm_mean_rel_pct']:.3f}/{b_fl['pmns_mean_rel_pct']:.3f}"
        ),
        (
            f"- GW auc/adv/sep/gap: "
            f"{b_gw['auc_h1l1_vs_ctrl']:.4f}/{b_gw['adv_shared_minus_ctrl_q90']:.4f}/"
            f"{b_gw['sep_median_h1l1_minus_ctrl']:.6f}/{b_gw['control_median_gap']:.6f}"
        ),
        f"- all_pass: {best['all_pass']}",
        "",
        "## Required Next Step",
        f"- {required_next}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2023] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2023] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2023] verdict={verdict} feasible_count={feasible_count} evals={eval_count}")


if __name__ == "__main__":
    main()
