#!/usr/bin/env python3
"""
QW-2024: Eta-kernel isospin-split shared flavor rework scan.

Purpose:
- keep one frozen kernel (from QW-2021),
- test richer but still shared flavor dynamics (isospin/sector split),
- check if triple-sector blockers from QW-2023 can be reduced.
"""

from __future__ import annotations

import itertools
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2024_eta_kernel_isospin_flavor_rework_scan.json"
OUT_MD = ROOT / "RAPORT_QW2024_ETA_KERNEL_ISOSPIN_FLAVOR_REWORK_SCAN.md"


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


def gw_metrics(p_amp: float, r_dist: float, kernel: Dict[str, float], df_gw: pd.DataFrame) -> Dict[str, float]:
    d = np.array([1.0, 2.0, 3.0, 4.0], dtype=float)
    raw_w = (np.abs(kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])) ** p_amp) * (d**r_dist)
    w = raw_w / np.sum(raw_w)

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


def mass_metrics(
    p_amp: float,
    r_dist: float,
    phase_scale: float,
    gamma_scale: float,
    kernel: Dict[str, float],
) -> Dict[str, float]:
    particles = [
        ("Top", 0.0, 173_000.0, 0.0, 3.0),
        ("Bottom", 7.0, 4_180.0, 1.0, 3.0),
        ("Tau", 9.0, 1_776.9, -1.0, 3.0),
        ("Charm", 9.0, 1_270.0, 1.0, 2.0),
        ("Muon", 14.0, 105.7, -1.0, 2.0),
        ("Electron", 24.0, 0.511, -1.0, 1.0),
    ]

    d1, d4 = np.array([1.0]), np.array([4.0])
    k1 = float(abs(kernel_fn(d1, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])[0]))
    k4 = float(abs(kernel_fn(d4, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])[0]))
    gamma_kernel = float(-4.0 * math.log(max(k4 / max(k1, 1e-15), 1e-15), 4.0) / 3.0)
    gamma_eff = float(np.clip(gamma_scale * gamma_kernel, 0.40, 2.60))

    c_q = float(np.tanh(0.90 * (p_amp - 1.0)))
    c_s = float(np.tanh(0.70 * r_dist))
    c_g = float(np.tanh(0.85 * (phase_scale - 1.0)))

    errs = []
    for _, q, m_exp, sector, gen in particles:
        base = 173_000.0 * (4.0 ** (-(gamma_eff * q / 4.0)))
        delta = c_q * (q / 24.0) + c_s * sector + c_g * (gen - 2.0)
        pred = float(base * math.exp(delta))
        err = float(abs(pred - m_exp) / max(m_exp, 1e-15) * 100.0)
        errs.append(err)

    return {
        "mean_rel_err_pct": float(np.mean(errs)),
        "max_rel_err_pct": float(np.max(errs)),
        "gamma_eff": gamma_eff,
    }


def flavor_hamiltonian(
    q: np.ndarray,
    iso_tag: float,
    sector_tag: float,
    p_amp: float,
    r_dist: float,
    params: Dict[str, float],
    kernel: Dict[str, float],
) -> np.ndarray:
    n = len(q)
    d = 1.0 + cyclic_distance_matrix(q, q, modulus=24)
    k = kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])
    base = np.sign(k) * (np.abs(k) ** p_amp) * (d**r_dist)

    idx = np.arange(n, dtype=float)
    i_minus_j = idx[:, None] - idx[None, :]
    gap = np.abs(i_minus_j)
    q_diff = q[:, None] - q[None, :]
    near = np.exp(-params["rho_gap"] * gap)

    phase = (
        params["phase_q"] * q_diff
        + iso_tag * params["theta_iso"] * i_minus_j
        + sector_tag * params["theta_sector"] * i_minus_j
    )
    strength = params["lambda_mix"] * base * near

    re_raw = strength * np.cos(phase)
    im_raw = strength * np.sin(phase)

    re = 0.5 * (re_raw + re_raw.T)
    im_asym = 0.5 * (im_raw - im_raw.T)
    np.fill_diagonal(re, 0.0)
    np.fill_diagonal(im_asym, 0.0)

    idx_centered = idx - np.mean(idx)
    q_centered = q - np.mean(q)
    diag = (
        params["diag_q_coeff"] * q_centered
        + iso_tag * params["diag_iso"] * idx_centered
        + sector_tag * params["diag_sector"] * idx_centered
    )

    h = re + 1j * params["chi_im"] * im_asym + np.diag(diag)
    return 0.5 * (h + h.conj().T)


def flavor_metrics(
    q_up: np.ndarray,
    q_down: np.ndarray,
    q_nu: np.ndarray,
    q_lep: np.ndarray,
    p_amp: float,
    r_dist: float,
    params: Dict[str, float],
    kernel: Dict[str, float],
) -> Dict[str, float]:
    hu = flavor_hamiltonian(q_up, iso_tag=+1.0, sector_tag=+1.0, p_amp=p_amp, r_dist=r_dist, params=params, kernel=kernel)
    hd = flavor_hamiltonian(q_down, iso_tag=-1.0, sector_tag=+1.0, p_amp=p_amp, r_dist=r_dist, params=params, kernel=kernel)
    hn = flavor_hamiltonian(q_nu, iso_tag=+1.0, sector_tag=-1.0, p_amp=p_amp, r_dist=r_dist, params=params, kernel=kernel)
    hl = flavor_hamiltonian(q_lep, iso_tag=-1.0, sector_tag=-1.0, p_amp=p_amp, r_dist=r_dist, params=params, kernel=kernel)

    _, uu = np.linalg.eigh(hu)
    _, ud = np.linalg.eigh(hd)
    _, un = np.linalg.eigh(hn)
    _, ul = np.linalg.eigh(hl)

    ckm = np.abs(uu.conj().T @ ud)
    pmns = np.abs(un.conj().T @ ul)

    ckm_rel = np.abs(ckm - CKM_REF) / np.clip(CKM_REF, 1e-12, None)
    pmns_rel = np.abs(pmns - PMNS_REF) / np.clip(PMNS_REF, 1e-12, None)
    return {
        "ckm_mean_rel_pct": float(100.0 * np.mean(ckm_rel)),
        "pmns_mean_rel_pct": float(100.0 * np.mean(pmns_rel)),
        "ckm_max_rel_pct": float(100.0 * np.max(ckm_rel)),
        "pmns_max_rel_pct": float(100.0 * np.max(pmns_rel)),
    }


def main() -> None:
    d2021 = json.loads((ROOT / "report_qw2021_v2_eta_operator_beta_constraint_scan.json").read_text(encoding="utf-8"))
    kernel = d2021["selected"]["fit"]

    df_gw = pd.read_csv(ROOT / "gw1831_window_features.csv")

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

    phase_scale_fixed = 0.83
    gamma_scale_fixed = 0.7166666666666667

    q_up = np.array([0.0, 9.0, 14.0], dtype=float)
    q_down = np.array([7.0, 9.0, 14.0], dtype=float)
    q_lep = np.array([24.0, 14.0, 9.0], dtype=float)
    q_nu_space = list(itertools.permutations([0.0, 1.0, 2.0], 3))

    grid = {
        "p_amp": np.linspace(0.7, 1.7, 5).tolist(),
        "r_dist": np.linspace(-0.2, 0.8, 5).tolist(),
        "lambda_mix": [0.8, 1.2, 1.6],
        "rho_gap": [0.0, 0.6],
        "chi_im": [0.0, 0.3],
        "phase_q": [0.05, 0.2],
        "theta_iso": [0.0, 0.15, 0.3],
        "theta_sector": [0.0, 0.15],
        "diag_q_coeff": [0.0, 0.1],
        "diag_iso": [0.0, 0.08],
        "diag_sector": [0.0, 0.08],
    }

    total = (
        len(grid["p_amp"]) * len(grid["r_dist"]) * len(grid["lambda_mix"]) * len(grid["rho_gap"]) *
        len(grid["chi_im"]) * len(grid["phase_q"]) * len(grid["theta_iso"]) * len(grid["theta_sector"]) *
        len(grid["diag_q_coeff"]) * len(grid["diag_iso"]) * len(grid["diag_sector"]) * len(q_nu_space)
    )

    best = None
    pass_count = 0

    for p_amp in grid["p_amp"]:
        for r_dist in grid["r_dist"]:
            g = gw_metrics(float(p_amp), float(r_dist), kernel=kernel, df_gw=df_gw)
            m = mass_metrics(float(p_amp), float(r_dist), phase_scale=phase_scale_fixed, gamma_scale=gamma_scale_fixed, kernel=kernel)

            for lambda_mix in grid["lambda_mix"]:
                for rho_gap in grid["rho_gap"]:
                    for chi_im in grid["chi_im"]:
                        for phase_q in grid["phase_q"]:
                            for theta_iso in grid["theta_iso"]:
                                for theta_sector in grid["theta_sector"]:
                                    for diag_q_coeff in grid["diag_q_coeff"]:
                                        for diag_iso in grid["diag_iso"]:
                                            for diag_sector in grid["diag_sector"]:
                                                params = {
                                                    "lambda_mix": float(lambda_mix),
                                                    "rho_gap": float(rho_gap),
                                                    "chi_im": float(chi_im),
                                                    "phase_q": float(phase_q),
                                                    "theta_iso": float(theta_iso),
                                                    "theta_sector": float(theta_sector),
                                                    "diag_q_coeff": float(diag_q_coeff),
                                                    "diag_iso": float(diag_iso),
                                                    "diag_sector": float(diag_sector),
                                                }

                                                for q_nu in q_nu_space:
                                                    f = flavor_metrics(
                                                        q_up=q_up,
                                                        q_down=q_down,
                                                        q_nu=np.array(q_nu, dtype=float),
                                                        q_lep=q_lep,
                                                        p_amp=float(p_amp),
                                                        r_dist=float(r_dist),
                                                        params=params,
                                                        kernel=kernel,
                                                    )

                                                    flags = {
                                                        "mass_mean_rel_pct_le_max": bool(m["mean_rel_err_pct"] <= thresholds["mass_mean_rel_pct_max"]),
                                                        "mass_max_rel_pct_le_max": bool(m["max_rel_err_pct"] <= thresholds["mass_max_rel_pct_max"]),
                                                        "ckm_mean_rel_pct_le_max": bool(f["ckm_mean_rel_pct"] <= thresholds["ckm_mean_rel_pct_max"]),
                                                        "pmns_mean_rel_pct_le_max": bool(f["pmns_mean_rel_pct"] <= thresholds["pmns_mean_rel_pct_max"]),
                                                        "gw_sep_ge_min": bool(g["sep_median_h1l1_minus_ctrl"] >= thresholds["gw_sep_min"]),
                                                        "gw_adv_ge_min": bool(g["adv_shared_minus_ctrl_q90"] >= thresholds["gw_adv_min"]),
                                                        "gw_auc_ge_min": bool(g["auc_h1l1_vs_ctrl"] >= thresholds["gw_auc_min"]),
                                                        "gw_control_gap_le_max": bool(g["control_median_gap"] <= thresholds["gw_control_gap_max"]),
                                                    }
                                                    all_pass = bool(all(flags.values()))
                                                    if all_pass:
                                                        pass_count += 1

                                                    score = (
                                                        f["ckm_mean_rel_pct"] / thresholds["ckm_mean_rel_pct_max"]
                                                        + f["pmns_mean_rel_pct"] / thresholds["pmns_mean_rel_pct_max"]
                                                        + 0.06 * (m["mean_rel_err_pct"] / thresholds["mass_mean_rel_pct_max"])
                                                        + 0.02 * (m["max_rel_err_pct"] / thresholds["mass_max_rel_pct_max"])
                                                        + 30.0 * max(0.0, thresholds["gw_auc_min"] - g["auc_h1l1_vs_ctrl"])
                                                        + 35.0 * max(0.0, thresholds["gw_adv_min"] - g["adv_shared_minus_ctrl_q90"])
                                                        + 220.0 * max(0.0, thresholds["gw_sep_min"] - g["sep_median_h1l1_minus_ctrl"])
                                                        + 24.0 * max(0.0, g["control_median_gap"] - thresholds["gw_control_gap_max"])
                                                    )

                                                    row = {
                                                        "p_amp": float(p_amp),
                                                        "r_dist": float(r_dist),
                                                        "q_nu": [int(x) for x in q_nu],
                                                        "flavor_params": params,
                                                        "mass": m,
                                                        "flavor": f,
                                                        "gw": g,
                                                        "flags": flags,
                                                        "all_pass": all_pass,
                                                        "score": float(score),
                                                    }

                                                    if best is None or row["score"] < best["score"]:
                                                        best = row

    verdict = "ETA_KERNEL_ISOSPIN_FLAVOR_REWORK_PASS" if pass_count > 0 else "ETA_KERNEL_ISOSPIN_FLAVOR_REWORK_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw2021_v2_eta_operator_beta_constraint_scan.json:selected.fit",
        "kernel": kernel,
        "thresholds": thresholds,
        "search_space_size": int(total),
        "pass_count_all_flags": int(pass_count),
        "best_row": best,
        "verdict": verdict,
        "required_next_step": (
            "RERUN_TRIPLE_SECTOR_GATE_WITH_ISOSPIN_OPERATOR"
            if pass_count > 0
            else "ESCALATE_TO_MASS_OPERATOR_REFORMULATION_WITH_FLAVOR_REWORK"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    b = best
    lines = [
        "# RAPORT QW-2024: ETA-KERNEL ISOSPIN FLAVOR REWORK SCAN",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- search_space_size: {total}",
        f"- pass_count_all_flags: {pass_count}",
        "",
        "## Best Row",
        f"- p_amp/r_dist: {b['p_amp']:.3f}/{b['r_dist']:.3f}",
        f"- q_nu: {b['q_nu']}",
        f"- mass mean/max rel%: {b['mass']['mean_rel_err_pct']:.3f}/{b['mass']['max_rel_err_pct']:.3f}",
        f"- flavor CKM/PMNS mean rel%: {b['flavor']['ckm_mean_rel_pct']:.3f}/{b['flavor']['pmns_mean_rel_pct']:.3f}",
        f"- GW auc/adv/sep/gap: {b['gw']['auc_h1l1_vs_ctrl']:.4f}/{b['gw']['adv_shared_minus_ctrl_q90']:.4f}/{b['gw']['sep_median_h1l1_minus_ctrl']:.6f}/{b['gw']['control_median_gap']:.6f}",
        f"- all_pass: {b['all_pass']}",
        "",
        "## Required Next Step",
        f"- {out['required_next_step']}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2024] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2024] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2024] verdict={verdict} pass_count={pass_count}/{total}")


if __name__ == "__main__":
    main()
