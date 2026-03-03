#!/usr/bin/env python3
"""
QW-1967: Local refinement gate for isospin-split shared flavor dynamics.

Purpose:
- refine around near-threshold region found in QW-1966,
- keep fixed mass branch and single shared operator class,
- deterministic Monte Carlo + full q_nu permutation check.
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
OUT_JSON = ROOT / "report_qw1967_isospin_split_local_refinement_gate.json"
OUT_MD = ROOT / "RAPORT_QW1967_ISOSPIN_SPLIT_LOCAL_REFINEMENT_GATE.md"


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
    iso_tag: float,
    sector_tag: float,
    params: Dict[str, float],
    kernel: Dict[str, float],
) -> np.ndarray:
    n = len(q)
    d = 1.0 + cyclic_distance_matrix(q, q, modulus=24)
    k = kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])
    base = np.sign(k) * (np.abs(k) ** params["p_amp"]) * (d**params["r_dist"])

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
    params: Dict[str, float],
    kernel: Dict[str, float],
) -> Dict[str, float]:
    hu = flavor_hamiltonian(q_up, iso_tag=+1.0, sector_tag=+1.0, params=params, kernel=kernel)
    hd = flavor_hamiltonian(q_down, iso_tag=-1.0, sector_tag=+1.0, params=params, kernel=kernel)
    hn = flavor_hamiltonian(q_nu, iso_tag=+1.0, sector_tag=-1.0, params=params, kernel=kernel)
    hl = flavor_hamiltonian(q_lep, iso_tag=-1.0, sector_tag=-1.0, params=params, kernel=kernel)

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


def sample_params(rng: np.random.Generator) -> Dict[str, float]:
    return {
        "p_amp": float(rng.uniform(0.72, 0.95)),
        "r_dist": float(rng.uniform(-0.10, 1.00)),
        "lambda_mix": float(rng.uniform(0.60, 1.40)),
        "rho_gap": float(rng.uniform(0.20, 1.20)),
        "chi_im": float(rng.uniform(0.00, 0.40)),
        "phase_q": float(rng.uniform(0.12, 0.28)),
        "theta_iso": float(rng.uniform(0.18, 0.36)),
        "theta_sector": float(rng.uniform(0.05, 0.25)),
        "diag_q_coeff": float(rng.uniform(0.00, 0.08)),
        "diag_iso": float(rng.uniform(0.00, 0.06)),
        "diag_sector": float(rng.uniform(0.04, 0.12)),
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

    q_up = np.array([q_assign["Top"], q_assign["Charm"], q_assign["Muon"]], dtype=float)
    q_down = np.array([q_assign["Bottom"], q_assign["Tau"], q_assign["Muon"]], dtype=float)
    q_lep = np.array([q_assign["Electron"], q_assign["Muon"], q_assign["Tau"]], dtype=float)
    q_nu_space = list(itertools.permutations([0, 1, 2], 3))

    thresholds = {
        "ckm_mean_rel_pct_max": 15.0,
        "pmns_mean_rel_pct_max": 15.0,
        "gw_sep_min": 0.0020,
        "gw_adv_min": 0.30,
        "gw_auc_min": 0.75,
        "gw_control_gap_max": 0.0025,
    }

    n_param_samples = 15000
    n_total = n_param_samples * len(q_nu_space)
    rng = np.random.default_rng(1967)
    df_gw = pd.read_csv(ROOT / "gw1831_window_features.csv")

    best = None
    pass_count = 0
    rows: List[Dict[str, object]] = []

    for _ in range(n_param_samples):
        params = sample_params(rng)
        g = gw_metrics(params["p_amp"], params["r_dist"], kernel, df_gw)
        for q_nu_tuple in q_nu_space:
            q_nu = np.array(q_nu_tuple, dtype=float)
            f = flavor_metrics(q_up, q_down, q_nu, q_lep, params, kernel)

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
        "ISOSPIN_LOCAL_REFINEMENT_PASS"
        if pass_count > 0
        else "ISOSPIN_LOCAL_REFINEMENT_FAIL"
    )
    required_next = (
        "LOCK_REFINED_SHARED_OPERATOR_AND_RUN_EXTERNAL_CONFIRMATORY"
        if pass_count > 0
        else "NO_PASS_IN_LOCAL_REFINEMENT_NEED_NEW_FLAVOR_PRINCIPLE"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_mass_branch": {
            "source_report": "report_qw1962_noncircular_branch_unified_triad_gate.json",
            "mass_metrics_fixed": mass_metrics_fixed,
            "mass_flags_fixed": mass_flags,
        },
        "kernel": kernel,
        "search": {
            "mode": "deterministic_monte_carlo",
            "seed": 1967,
            "n_param_samples": int(n_param_samples),
            "q_nu_space_size": len(q_nu_space),
            "n_total": int(n_total),
            "ranges": {
                "p_amp": [0.72, 0.95],
                "r_dist": [-0.10, 1.00],
                "lambda_mix": [0.60, 1.40],
                "rho_gap": [0.20, 1.20],
                "chi_im": [0.00, 0.40],
                "phase_q": [0.12, 0.28],
                "theta_iso": [0.18, 0.36],
                "theta_sector": [0.05, 0.25],
                "diag_q_coeff": [0.00, 0.08],
                "diag_iso": [0.00, 0.06],
                "diag_sector": [0.04, 0.12],
            },
        },
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
        "# RAPORT QW-1967: ISOSPIN SPLIT LOCAL REFINEMENT GATE",
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
        f"- q_nu: {b['q_nu']}",
        f"- all_pass: {b['all_pass']}",
        "",
        "## Required Next Step",
        f"- {required_next}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1967] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1967] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1967] verdict={verdict} pass_count={pass_count}/{n_total}")


if __name__ == "__main__":
    main()
