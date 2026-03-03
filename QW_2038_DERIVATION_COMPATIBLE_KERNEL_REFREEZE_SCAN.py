#!/usr/bin/env python3
"""
QW-2038: Derivation-compatible kernel refreeze scan.

Goal:
- search for a kernel (omega, phi, beta, eta) with beta<=1.0
  that still passes mass+flavor+GW under the same shared context
  (no sector-specific retune of operators).
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
OUT_JSON = ROOT / "report_qw2038_derivation_compatible_kernel_refreeze_scan.json"
OUT_MD = ROOT / "RAPORT_QW2038_DERIVATION_COMPATIBLE_KERNEL_REFREEZE_SCAN.md"


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


def gw_metrics(
    p_amp: float,
    r_dist: float,
    kernel: Dict[str, float],
    df_gw: pd.DataFrame,
    kappa: float,
) -> Dict[str, float]:
    d = np.array([1.0, 2.0, 3.0, 4.0], dtype=float)
    raw_w = (np.abs(kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])) ** p_amp) * (
        d**r_dist
    )
    w = raw_w / np.sum(raw_w)

    score = (
        w[0] * df_gw["max_abs_corr"].to_numpy(dtype=float)
        + w[1] * df_gw["mean_abs_corr"].to_numpy(dtype=float)
        + w[2] * df_gw["corr_at_0ms"].to_numpy(dtype=float)
        + w[3] * df_gw["corr_at_10ms"].to_numpy(dtype=float)
    )

    pair = df_gw["pair"].astype(str).to_numpy()
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


def mass_metrics(
    p_amp: float,
    r_dist: float,
    gamma_scale: float,
    coeffs: Dict[str, float],
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

    q_diff = q[:, None] - q[None, :]
    q_abs = np.abs(q_diff) / 24.0
    base_amp = np.sign(k) * (np.abs(k) ** p_amp) * (d**r_dist)
    amp_mod = 1.0 + params["amp_qbias"] * (q_abs - float(np.mean(q_abs)))
    base = base_amp * amp_mod

    idx = np.arange(n, dtype=float)
    i_minus_j = idx[:, None] - idx[None, :]
    gap = np.abs(i_minus_j)
    near = np.exp(-params["rho_gap"] * gap)

    phase = (
        params["phase_q"] * q_diff
        + params["phase_q3"] * ((q_diff / 24.0) ** 3)
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
    return {
        "ckm_mean_rel_pct": float(100.0 * np.mean(ckm_rel)),
        "pmns_mean_rel_pct": float(100.0 * np.mean(pmns_rel)),
        "ckm_max_rel_pct": float(100.0 * np.max(ckm_rel)),
        "pmns_max_rel_pct": float(100.0 * np.max(pmns_rel)),
    }


def main() -> None:
    d2030 = json.loads((ROOT / "report_qw2030_final_stage_c_gate_combined_branch.json").read_text(encoding="utf-8"))
    d2028 = json.loads((ROOT / "report_qw2028_joint_scan_with_gw_kappa_term.json").read_text(encoding="utf-8"))
    d2029 = json.loads((ROOT / "report_qw2029_ckm_blocker_shared_flavor_basis_scan.json").read_text(encoding="utf-8"))
    d2027 = json.loads((ROOT / "report_qw2027_gw_control_gap_structural_term_scan.json").read_text(encoding="utf-8"))

    kernel0 = d2030["kernel"]
    best2028 = d2028["best_row"]
    best2029 = d2029["best_row"]
    kappa = float(d2027["selected"]["kappa"])

    p_amp = float(best2028["p_amp"])
    r_dist = float(best2028["r_dist"])
    gamma_scale = float(best2028["gamma_scale"])
    coeffs = best2028["coeffs"]
    q_nu = np.array(best2029["q_nu"], dtype=float)
    params = best2029["params"]

    q_up = np.array([0.0, 9.0, 14.0], dtype=float)
    q_down = np.array([7.0, 9.0, 14.0], dtype=float)
    q_lep = np.array([24.0, 14.0, 9.0], dtype=float)
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

    omega0 = float(kernel0["omega"])
    phi0 = float(kernel0["phi"])
    eta0 = float(kernel0["eta"])
    beta0 = float(kernel0["beta"])

    omega_grid = np.linspace(max(0.05, omega0 - 0.08), omega0 + 0.08, 11)
    phi_grid = np.linspace(phi0 - 0.30, phi0 + 0.30, 13)
    beta_grid = np.linspace(0.60, 1.00, 11)
    eta_grid = np.array([1.5, 1.6, 1.7, 1.8], dtype=float)
    total = int(len(omega_grid) * len(phi_grid) * len(beta_grid) * len(eta_grid))

    best = None
    pass_count = 0
    pass_rows = []

    for om, ph, be, et in itertools.product(omega_grid, phi_grid, beta_grid, eta_grid):
        kernel = {
            "omega": float(om),
            "phi": float(ph),
            "beta": float(be),
            "eta": float(et),
        }
        m = mass_metrics(p_amp=p_amp, r_dist=r_dist, gamma_scale=gamma_scale, coeffs=coeffs, kernel=kernel)
        f = flavor_metrics(
            q_up=q_up,
            q_down=q_down,
            q_nu=q_nu,
            q_lep=q_lep,
            p_amp=p_amp,
            r_dist=r_dist,
            params=params,
            kernel=kernel,
        )
        g = gw_metrics(p_amp=p_amp, r_dist=r_dist, kernel=kernel, df_gw=df_gw, kappa=kappa)

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
            m["mean_rel_err_pct"] / thresholds["mass_mean_rel_pct_max"]
            + f["ckm_mean_rel_pct"] / thresholds["ckm_mean_rel_pct_max"]
            + f["pmns_mean_rel_pct"] / thresholds["pmns_mean_rel_pct_max"]
            + max(0.0, (thresholds["gw_sep_min"] - g["sep_median_h1l1_minus_ctrl"]) / thresholds["gw_sep_min"])
            + max(0.0, (thresholds["gw_adv_min"] - g["adv_shared_minus_ctrl_q90"]) / thresholds["gw_adv_min"])
            + max(0.0, (thresholds["gw_auc_min"] - g["auc_h1l1_vs_ctrl"]) / thresholds["gw_auc_min"])
            + max(0.0, (g["control_median_gap"] - thresholds["gw_control_gap_max"]) / thresholds["gw_control_gap_max"])
            + 0.05 * abs(be - beta0)
            + 0.03 * abs(et - eta0)
        )

        row = {
            "kernel": kernel,
            "mass": m,
            "flavor": f,
            "gw": g,
            "flags": flags,
            "all_pass": all_pass,
            "score": float(score),
        }
        if all_pass and len(pass_rows) < 25:
            pass_rows.append(row)
        if best is None or row["score"] < best["score"]:
            best = row

    if pass_count > 0:
        verdict = "DERIVATION_COMPATIBLE_KERNEL_REFREEZE_PASS"
        readiness = "KERNEL_REFREEZE_CANDIDATE_AVAILABLE"
        required_next = "RUN_QW2030_STYLE_GATE_ON_NEW_CANDIDATE_KERNEL"
    else:
        verdict = "DERIVATION_COMPATIBLE_KERNEL_REFREEZE_FAIL"
        readiness = "NO_BETA_LE_1P0_KERNEL_PASSES_CURRENT_TRIAD"
        required_next = "EXTEND_MICRO_DYNAMICS_OR_OPERATOR_FOR_BETA_IDENTIFIABILITY"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "kernel0": "report_qw2030_final_stage_c_gate_combined_branch.json:kernel",
            "context_2028": "report_qw2028_joint_scan_with_gw_kappa_term.json:best_row",
            "context_2029": "report_qw2029_ckm_blocker_shared_flavor_basis_scan.json:best_row",
            "kappa": "report_qw2027_gw_control_gap_structural_term_scan.json:selected.kappa",
        },
        "fixed_context": {
            "p_amp": p_amp,
            "r_dist": r_dist,
            "gamma_scale": gamma_scale,
            "coeffs": coeffs,
            "q_nu": [float(x) for x in q_nu],
            "flavor_params": params,
            "kappa": kappa,
        },
        "grid": {
            "omega_grid": [float(x) for x in omega_grid],
            "phi_grid": [float(x) for x in phi_grid],
            "beta_grid": [float(x) for x in beta_grid],
            "eta_grid": [float(x) for x in eta_grid],
            "n_total": total,
        },
        "thresholds": thresholds,
        "pass_count": int(pass_count),
        "best_row": best,
        "pass_examples": pass_rows,
        "verdict": verdict,
        "readiness": readiness,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2038: DERIVATION-COMPATIBLE KERNEL REFREEZE SCAN",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Readiness: **{readiness}**",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {pass_count}/{total}",
        "",
        "## Best Row",
        (
            f"- kernel omega/phi/beta/eta: {best['kernel']['omega']:.6f} / {best['kernel']['phi']:.6f} / "
            f"{best['kernel']['beta']:.6f} / {best['kernel']['eta']:.6f}"
        ),
        (
            f"- mass mean/max rel%: {best['mass']['mean_rel_err_pct']:.3f} / "
            f"{best['mass']['max_rel_err_pct']:.3f}"
        ),
        (
            f"- flavor CKM/PMNS mean rel%: {best['flavor']['ckm_mean_rel_pct']:.3f} / "
            f"{best['flavor']['pmns_mean_rel_pct']:.3f}"
        ),
        (
            f"- GW auc/adv/sep/gap: {best['gw']['auc_h1l1_vs_ctrl']:.4f} / "
            f"{best['gw']['adv_shared_minus_ctrl_q90']:.4f} / "
            f"{best['gw']['sep_median_h1l1_minus_ctrl']:.6f} / "
            f"{best['gw']['control_median_gap']:.6f}"
        ),
        "",
        "## Required Next Step",
        f"- {required_next}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2038] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2038] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2038] readiness={readiness} verdict={verdict} pass={pass_count}/{total}")


if __name__ == "__main__":
    main()
