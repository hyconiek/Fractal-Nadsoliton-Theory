#!/usr/bin/env python3
"""
QW-2028: Joint mass+flavor+GW scan on eta-kernel with GW kappa term.
"""

from __future__ import annotations

import itertools
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2028_joint_scan_with_gw_kappa_term.json"
OUT_MD = ROOT / "RAPORT_QW2028_JOINT_SCAN_WITH_GW_KAPPA_TERM.md"

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


def gw_metrics(p_amp: float, r_dist: float, kernel: Dict[str, float], df_gw: pd.DataFrame, kappa: float) -> Dict[str, float]:
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
    # Apply structural control correction.
    score = score.copy()
    score[pair == "H1-V1"] -= float(kappa)
    score[pair == "L1-V1"] += float(kappa)

    s_hl = score[pair == "H1-L1"]
    s_hv = score[pair == "H1-V1"]
    s_lv = score[pair == "L1-V1"]
    s_ctrl = np.concatenate([s_hv, s_lv])

    q90 = float(np.quantile(s_ctrl, 0.90))
    auc = rank_auc_pos_gt_neg(s_hl, s_ctrl)
    adv = float(np.mean(s_hl > q90) - np.mean(s_ctrl > q90))
    sep = float(np.median(s_hl) - np.median(s_ctrl))
    gap = float(abs(np.median(s_hv) - np.median(s_lv)))

    return {
        "auc_h1l1_vs_ctrl": auc,
        "adv_shared_minus_ctrl_q90": adv,
        "sep_median_h1l1_minus_ctrl": sep,
        "control_median_gap": gap,
    }


def mass_metrics(p_amp: float, r_dist: float, gamma_scale: float, coeffs: Dict[str, float], kernel: Dict[str, float]) -> Dict[str, float]:
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

    errs = []
    for _, q, m_exp, sector, gen in particles:
        xq = q / 24.0
        xg = gen - 2.0
        base = 173_000.0 * (4.0 ** (-(gamma_eff * q / 4.0)))
        delta = (
            coeffs["c_q"] * xq
            + coeffs["c_s"] * sector
            + coeffs["c_g"] * xg
            + coeffs["c_q2"] * (xq**2)
            + coeffs["c_sg"] * sector * xg
        )
        pred = float(base * math.exp(delta))
        err = float(abs(pred - m_exp) / max(m_exp, 1e-15) * 100.0)
        errs.append(err)

    return {"mean_rel_err_pct": float(np.mean(errs)), "max_rel_err_pct": float(np.max(errs))}


def flavor_hamiltonian(q: np.ndarray, iso_tag: float, sector_tag: float, p_amp: float, r_dist: float, params: Dict[str, float], kernel: Dict[str, float]) -> np.ndarray:
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


def flavor_metrics(q_up: np.ndarray, q_down: np.ndarray, q_nu: np.ndarray, q_lep: np.ndarray, p_amp: float, r_dist: float, params: Dict[str, float], kernel: Dict[str, float]) -> Dict[str, float]:
    hu = flavor_hamiltonian(q_up, +1.0, +1.0, p_amp, r_dist, params, kernel)
    hd = flavor_hamiltonian(q_down, -1.0, +1.0, p_amp, r_dist, params, kernel)
    hn = flavor_hamiltonian(q_nu, +1.0, -1.0, p_amp, r_dist, params, kernel)
    hl = flavor_hamiltonian(q_lep, -1.0, -1.0, p_amp, r_dist, params, kernel)

    _, uu = np.linalg.eigh(hu)
    _, ud = np.linalg.eigh(hd)
    _, un = np.linalg.eigh(hn)
    _, ul = np.linalg.eigh(hl)

    ckm = np.abs(uu.conj().T @ ud)
    pmns = np.abs(un.conj().T @ ul)

    ckm_rel = np.abs(ckm - CKM_REF) / np.clip(CKM_REF, 1e-12, None)
    pmns_rel = np.abs(pmns - PMNS_REF) / np.clip(PMNS_REF, 1e-12, None)
    return {"ckm_mean_rel_pct": float(100.0 * np.mean(ckm_rel)), "pmns_mean_rel_pct": float(100.0 * np.mean(pmns_rel))}


def main() -> None:
    kernel = json.loads((ROOT / "report_qw2021_v2_eta_operator_beta_constraint_scan.json").read_text(encoding="utf-8"))["selected"]["fit"]
    d2027 = json.loads((ROOT / "report_qw2027_gw_control_gap_structural_term_scan.json").read_text(encoding="utf-8"))
    kappa = float(d2027["selected"]["kappa"])

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

    q_up = np.array([0.0, 9.0, 14.0], dtype=float)
    q_down = np.array([7.0, 9.0, 14.0], dtype=float)
    q_lep = np.array([24.0, 14.0, 9.0], dtype=float)
    q_nu_space = list(itertools.permutations([0.0, 1.0, 2.0], 3))

    grid = {
        "p_amp": [0.7, 0.8],
        "r_dist": [-0.2, 0.0],
        "gamma_scale": [0.65],
        "c_q": [-0.8, -0.6],
        "lambda_mix": [0.8, 1.2],
        "rho_gap": [0.6],
        "chi_im": [0.0, 0.3],
        "phase_q": [0.2],
        "theta_iso": [0.3, 0.45],
        "theta_sector": [0.0, 0.15],
        "diag_q_coeff": [0.0, 0.1],
    }

    total = (
        len(grid["p_amp"]) * len(grid["r_dist"]) * len(grid["gamma_scale"]) * len(grid["c_q"]) *
        len(grid["lambda_mix"]) * len(grid["rho_gap"]) * len(grid["chi_im"]) * len(grid["phase_q"]) *
        len(grid["theta_iso"]) * len(grid["theta_sector"]) * len(grid["diag_q_coeff"]) * len(q_nu_space)
    )

    best = None
    pass_count = 0

    for p_amp in grid["p_amp"]:
        for r_dist in grid["r_dist"]:
            g = gw_metrics(float(p_amp), float(r_dist), kernel, df_gw, kappa=kappa)
            for gamma_scale in grid["gamma_scale"]:
                for c_q in grid["c_q"]:
                    coeffs = {"c_q": float(c_q), "c_s": 0.0, "c_g": 0.0, "c_q2": 0.0, "c_sg": 0.0}
                    m = mass_metrics(float(p_amp), float(r_dist), float(gamma_scale), coeffs, kernel)

                    for lambda_mix in grid["lambda_mix"]:
                        for rho_gap in grid["rho_gap"]:
                            for chi_im in grid["chi_im"]:
                                for phase_q in grid["phase_q"]:
                                    for theta_iso in grid["theta_iso"]:
                                        for theta_sector in grid["theta_sector"]:
                                            for diag_q_coeff in grid["diag_q_coeff"]:
                                                params = {
                                                    "lambda_mix": float(lambda_mix),
                                                    "rho_gap": float(rho_gap),
                                                    "chi_im": float(chi_im),
                                                    "phase_q": float(phase_q),
                                                    "theta_iso": float(theta_iso),
                                                    "theta_sector": float(theta_sector),
                                                    "diag_q_coeff": float(diag_q_coeff),
                                                    "diag_iso": 0.0,
                                                    "diag_sector": 0.0,
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
                                                        + m["mean_rel_err_pct"] / thresholds["mass_mean_rel_pct_max"]
                                                        + 0.25 * m["max_rel_err_pct"] / thresholds["mass_max_rel_pct_max"]
                                                        + 24.0 * max(0.0, g["control_median_gap"] - thresholds["gw_control_gap_max"])
                                                    )
                                                    row = {
                                                        "p_amp": float(p_amp),
                                                        "r_dist": float(r_dist),
                                                        "gamma_scale": float(gamma_scale),
                                                        "coeffs": coeffs,
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

    verdict = "JOINT_SCAN_WITH_GW_KAPPA_TERM_PASS" if pass_count > 0 else "JOINT_SCAN_WITH_GW_KAPPA_TERM_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw2021_v2_eta_operator_beta_constraint_scan.json:selected.fit",
        "gw_kappa_source": "report_qw2027_gw_control_gap_structural_term_scan.json:selected.kappa",
        "kappa": kappa,
        "thresholds": thresholds,
        "search_space_size": int(total),
        "pass_count_all_flags": int(pass_count),
        "best_row": best,
        "verdict": verdict,
        "required_next_step": (
            "FREEZE_AND_RERUN_STAGE_C_GATE"
            if pass_count > 0
            else "TARGET_CKM_BLOCKER_WITH_ADDITIONAL_SHARED_FLAVOR_BASIS"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    b = best
    lines = [
        "# RAPORT QW-2028: JOINT SCAN WITH GW KAPPA TERM",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- kappa: {kappa:.6f}",
        f"- search_space_size: {total}",
        f"- pass_count_all_flags: {pass_count}",
        "",
        "## Best Row",
        f"- p_amp/r_dist: {b['p_amp']:.3f}/{b['r_dist']:.3f}",
        f"- gamma_scale: {b['gamma_scale']:.3f}",
        f"- c_q: {b['coeffs']['c_q']:.3f}",
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

    print(f"[QW-2028] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2028] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2028] verdict={verdict} pass_count={pass_count}/{total}")


if __name__ == "__main__":
    main()
