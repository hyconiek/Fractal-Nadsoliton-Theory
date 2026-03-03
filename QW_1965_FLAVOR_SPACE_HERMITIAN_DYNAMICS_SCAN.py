#!/usr/bin/env python3
"""
QW-1965: Shared Hermitian flavor-space dynamics scan on fixed noncircular mass branch.

Fixed:
- mass branch from QW-1962 (no edits to gamma/Q/delta_info),
- one frozen kernel.

Scanned (shared across all sectors):
- p_amp, r_dist (kernel envelope),
- chi_im (imaginary antisymmetric transport strength),
- signed_q_phase,
- diag_q_coeff.

Flavor construction:
- build one Hermitian H(q) per sector from the same rule,
- diagonalize H -> eigenvectors U,
- CKM = |U_u^\dagger U_d|, PMNS = |U_nu^\dagger U_lep|.
"""

from __future__ import annotations

import itertools
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1965_flavor_space_hermitian_dynamics_scan.json"
OUT_MD = ROOT / "RAPORT_QW1965_FLAVOR_SPACE_HERMITIAN_DYNAMICS_SCAN.md"


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


def flavor_hamiltonian(
    q: np.ndarray,
    p_amp: float,
    r_dist: float,
    chi_im: float,
    signed_q_phase: float,
    diag_q_coeff: float,
    omega: float,
    phi: float,
    beta: float,
    eta: float,
) -> np.ndarray:
    # Distances on cyclic topology.
    d = 1.0 + cyclic_distance_matrix(q, q, modulus=24)
    k = kernel_fn(d, omega=omega, phi=phi, beta=beta, eta=eta)
    amp = np.sign(k) * (np.abs(k) ** p_amp) * (d**r_dist)

    q_diff = q[:, None] - q[None, :]
    # Real symmetric carrier.
    h_re = amp * np.cos(signed_q_phase * q_diff)
    h_re = 0.5 * (h_re + h_re.T)
    # Imaginary antisymmetric transport term -> Hermitian i*A where A^T=-A.
    a_im = amp * np.sin(signed_q_phase * q_diff)
    a_im = 0.5 * (a_im - a_im.T)
    h_im = 1j * chi_im * a_im

    q_centered = q - np.mean(q)
    h_diag = np.diag(diag_q_coeff * q_centered)
    h = h_re + h_im + h_diag
    # Numerical guard for Hermiticity.
    return 0.5 * (h + h.conj().T)


def flavor_metrics(
    q_up: np.ndarray,
    q_down: np.ndarray,
    q_nu: np.ndarray,
    q_lep: np.ndarray,
    params: Dict[str, float],
    kernel: Dict[str, float],
) -> Dict[str, object]:
    hu = flavor_hamiltonian(q_up, **params, **kernel)
    hd = flavor_hamiltonian(q_down, **params, **kernel)
    hn = flavor_hamiltonian(q_nu, **params, **kernel)
    hl = flavor_hamiltonian(q_lep, **params, **kernel)

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
        "ckm_pred_abs": ckm.tolist(),
        "pmns_pred_abs": pmns.tolist(),
    }


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


def main() -> None:
    r1962 = json.loads((ROOT / "report_qw1962_noncircular_branch_unified_triad_gate.json").read_text(encoding="utf-8"))
    kernel = r1962["kernel"]
    q_assign = r1962["q_assignment_used"]
    mass_metrics_fixed = r1962["metrics"]["mass"]
    mass_flags = {
        "mass_mean_rel_pct_le_max": bool(r1962["flags"]["mass_mean_rel_pct_le_max"]),
        "mass_max_rel_pct_le_max": bool(r1962["flags"]["mass_max_rel_pct_le_max"]),
        "mass_tau_charm_ratio_err_le_max": bool(r1962["flags"]["mass_tau_charm_ratio_err_le_max"]),
    }

    # Keep same fixed branch mapping to avoid hidden retune.
    q_up = np.array([q_assign["Top"], q_assign["Charm"], q_assign["Muon"]], dtype=float)
    q_down = np.array([q_assign["Bottom"], q_assign["Tau"], q_assign["Muon"]], dtype=float)
    q_lep = np.array([q_assign["Electron"], q_assign["Muon"], q_assign["Tau"]], dtype=float)

    thresholds = {
        "ckm_mean_rel_pct_max": 15.0,
        "pmns_mean_rel_pct_max": 15.0,
        "gw_sep_min": 0.0020,
        "gw_adv_min": 0.30,
        "gw_auc_min": 0.75,
        "gw_control_gap_max": 0.0025,
    }

    grid = {
        "p_amp": [0.8, 1.2, 1.6],
        "r_dist": [-0.2, 0.4, 1.0],
        "chi_im": [0.0, 0.15, 0.30, 0.45],
        "signed_q_phase": [0.05, 0.15, 0.25, 0.35],
        "diag_q_coeff": [0.0, 0.1, 0.2],
    }
    q_nu_space = list(itertools.permutations([0, 1, 2, 3, 4], 3))
    n_total = (
        len(grid["p_amp"])
        * len(grid["r_dist"])
        * len(grid["chi_im"])
        * len(grid["signed_q_phase"])
        * len(grid["diag_q_coeff"])
        * len(q_nu_space)
    )

    df_gw = pd.read_csv(ROOT / "gw1831_window_features.csv")

    best = None
    pass_count = 0
    rows: List[Dict[str, object]] = []

    for p_amp in grid["p_amp"]:
        for r_dist in grid["r_dist"]:
            g = gw_metrics(p_amp=float(p_amp), r_dist=float(r_dist), kernel=kernel, df_gw=df_gw)
            for chi_im in grid["chi_im"]:
                for signed_q_phase in grid["signed_q_phase"]:
                    for diag_q_coeff in grid["diag_q_coeff"]:
                        params = {
                            "p_amp": float(p_amp),
                            "r_dist": float(r_dist),
                            "chi_im": float(chi_im),
                            "signed_q_phase": float(signed_q_phase),
                            "diag_q_coeff": float(diag_q_coeff),
                        }
                        for q_nu_tuple in q_nu_space:
                            q_nu = np.array(q_nu_tuple, dtype=float)
                            f = flavor_metrics(
                                q_up=q_up,
                                q_down=q_down,
                                q_nu=q_nu,
                                q_lep=q_lep,
                                params=params,
                                kernel=kernel,
                            )
                            flags = {
                                **mass_flags,
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
                                + 24.0 * max(0.0, g["control_median_gap"] - thresholds["gw_control_gap_max"])
                                + 30.0 * max(0.0, thresholds["gw_auc_min"] - g["auc_h1l1_vs_ctrl"])
                                + 35.0 * max(0.0, thresholds["gw_adv_min"] - g["adv_shared_minus_ctrl_q90"])
                                + 220.0 * max(0.0, thresholds["gw_sep_min"] - g["sep_median_h1l1_minus_ctrl"])
                            )
                            row = {
                                "params": params,
                                "q_nu": [int(x) for x in q_nu_tuple],
                                "flavor": f,
                                "gw": g,
                                "flags": flags,
                                "all_pass": all_pass,
                                "score": float(score),
                            }
                            rows.append(row)
                            if best is None or score < best["score"]:
                                best = row

    verdict = (
        "HERMITIAN_FLAVOR_DYNAMICS_PASS"
        if pass_count > 0
        else "HERMITIAN_FLAVOR_DYNAMICS_FAIL"
    )
    required_next = (
        "LOCK_HERMITIAN_OPERATOR_AND_RUN_STRICT_EXTERNAL_PACKAGE"
        if pass_count > 0
        else "FLAVOR_BLOCKER_IS_STRUCTURAL_BEYOND_CURRENT_SHARED_OPERATOR_CLASS"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_mass_branch": {
            "source_report": "report_qw1962_noncircular_branch_unified_triad_gate.json",
            "mass_metrics_fixed": mass_metrics_fixed,
            "mass_flags_fixed": mass_flags,
        },
        "kernel": kernel,
        "grid": {**grid, "q_nu_space_size": len(q_nu_space), "n_total": int(n_total)},
        "thresholds": thresholds,
        "summary": {
            "pass_count": int(pass_count),
            "best": best,
            "top10": sorted(rows, key=lambda x: x["score"])[:10],
        },
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    b = best
    lines = [
        "# RAPORT QW-1965: FLAVOR-SPACE HERMITIAN DYNAMICS SCAN",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {pass_count}/{n_total}",
        "",
        "## Fixed Mass Branch",
        (
            f"- mean/max/tau-charm err: "
            f"{mass_metrics_fixed['mean_rel_err_pct']:.3f}% / "
            f"{mass_metrics_fixed['max_rel_err_pct']:.3f}% / "
            f"{mass_metrics_fixed['tau_charm_ratio_rel_err_pct']:.3f}%"
        ),
        "",
        "## Best Candidate",
        (
            f"- params p/r/chi_im/signed_q_phase/diag_q_coeff: "
            f"{b['params']['p_amp']:.4f}/{b['params']['r_dist']:.4f}/"
            f"{b['params']['chi_im']:.4f}/{b['params']['signed_q_phase']:.4f}/"
            f"{b['params']['diag_q_coeff']:.4f}"
        ),
        f"- q_nu: {b['q_nu']}",
        (
            f"- flavor CKM/PMNS mean rel%: "
            f"{b['flavor']['ckm_mean_rel_pct']:.3f}/{b['flavor']['pmns_mean_rel_pct']:.3f}"
        ),
        (
            f"- GW auc/adv/sep/gap: "
            f"{b['gw']['auc_h1l1_vs_ctrl']:.4f}/"
            f"{b['gw']['adv_shared_minus_ctrl_q90']:.4f}/"
            f"{b['gw']['sep_median_h1l1_minus_ctrl']:.6f}/"
            f"{b['gw']['control_median_gap']:.6f}"
        ),
        f"- all_pass: {b['all_pass']}",
        "",
        "## Required Next Step",
        f"- {required_next}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1965] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1965] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1965] verdict={verdict} pass_count={pass_count}/{n_total}")


if __name__ == "__main__":
    main()
