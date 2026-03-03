#!/usr/bin/env python3
"""
QW-2029: Target CKM blocker with additional shared flavor basis.

Uses the best mass+GW-passing branch from QW-2028 and scans an extended
(shared) flavor phase basis to reduce CKM while preserving PMNS.
"""

from __future__ import annotations

import itertools
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2029_ckm_blocker_shared_flavor_basis_scan.json"
OUT_MD = ROOT / "RAPORT_QW2029_CKM_BLOCKER_SHARED_FLAVOR_BASIS_SCAN.md"


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

    # Extended shared phase basis (no per-sector free vectors).
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
    d2021 = json.loads((ROOT / "report_qw2021_v2_eta_operator_beta_constraint_scan.json").read_text(encoding="utf-8"))
    d2028 = json.loads((ROOT / "report_qw2028_joint_scan_with_gw_kappa_term.json").read_text(encoding="utf-8"))

    kernel = d2021["selected"]["fit"]
    b = d2028["best_row"]

    p_amp = float(b["p_amp"])
    r_dist = float(b["r_dist"])

    # Keep the already passing mass+GW context fixed.
    mass_context = b["mass"]
    gw_context = b["gw"]

    thresholds = {
        "ckm_mean_rel_pct_max": 15.0,
        "pmns_mean_rel_pct_max": 15.0,
    }

    q_up = np.array([0.0, 9.0, 14.0], dtype=float)
    q_down = np.array([7.0, 9.0, 14.0], dtype=float)
    q_lep = np.array([24.0, 14.0, 9.0], dtype=float)
    q_nu_space = list(itertools.permutations([0.0, 1.0, 2.0], 3))

    grid = {
        "lambda_mix": [0.6, 0.8, 1.0],
        "rho_gap": [0.4, 0.6, 0.8],
        "chi_im": [0.0, 0.15, 0.3],
        "phase_q": [0.15, 0.2, 0.25],
        "phase_q3": [-2.0, -1.0, 0.0, 1.0, 2.0],
        "theta_iso": [0.3, 0.45, 0.6],
        "theta_sector": [0.0, 0.15, 0.3],
        "diag_q_coeff": [0.0, 0.1, 0.2],
        "amp_qbias": [-0.4, 0.0, 0.4],
    }

    total = (
        len(grid["lambda_mix"]) * len(grid["rho_gap"]) * len(grid["chi_im"]) *
        len(grid["phase_q"]) * len(grid["phase_q3"]) * len(grid["theta_iso"]) *
        len(grid["theta_sector"]) * len(grid["diag_q_coeff"]) * len(grid["amp_qbias"]) * len(q_nu_space)
    )

    best = None
    pass_count = 0

    for q_nu in q_nu_space:
        for vals in itertools.product(
            grid["lambda_mix"],
            grid["rho_gap"],
            grid["chi_im"],
            grid["phase_q"],
            grid["phase_q3"],
            grid["theta_iso"],
            grid["theta_sector"],
            grid["diag_q_coeff"],
            grid["amp_qbias"],
        ):
            params = {
                "lambda_mix": float(vals[0]),
                "rho_gap": float(vals[1]),
                "chi_im": float(vals[2]),
                "phase_q": float(vals[3]),
                "phase_q3": float(vals[4]),
                "theta_iso": float(vals[5]),
                "theta_sector": float(vals[6]),
                "diag_q_coeff": float(vals[7]),
                "amp_qbias": float(vals[8]),
                "diag_iso": 0.0,
                "diag_sector": 0.0,
            }

            f = flavor_metrics(
                q_up=q_up,
                q_down=q_down,
                q_nu=np.array(q_nu, dtype=float),
                q_lep=q_lep,
                p_amp=p_amp,
                r_dist=r_dist,
                params=params,
                kernel=kernel,
            )

            flags = {
                "ckm_mean_rel_pct_le_max": bool(f["ckm_mean_rel_pct"] <= thresholds["ckm_mean_rel_pct_max"]),
                "pmns_mean_rel_pct_le_max": bool(f["pmns_mean_rel_pct"] <= thresholds["pmns_mean_rel_pct_max"]),
            }
            all_pass = bool(all(flags.values()))
            if all_pass:
                pass_count += 1

            score = (
                f["ckm_mean_rel_pct"] / thresholds["ckm_mean_rel_pct_max"]
                + f["pmns_mean_rel_pct"] / thresholds["pmns_mean_rel_pct_max"]
                + 0.05 * abs(params["phase_q3"])
                + 0.03 * abs(params["amp_qbias"])
            )

            row = {
                "q_nu": [int(x) for x in q_nu],
                "params": params,
                "flavor": f,
                "flags": flags,
                "all_pass": all_pass,
                "score": float(score),
            }

            if best is None or row["score"] < best["score"]:
                best = row

    verdict = "CKM_BLOCKER_SHARED_FLAVOR_BASIS_PASS" if pass_count > 0 else "CKM_BLOCKER_SHARED_FLAVOR_BASIS_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw2021_v2_eta_operator_beta_constraint_scan.json:selected.fit",
        "context_source": "report_qw2028_joint_scan_with_gw_kappa_term.json:best_row",
        "fixed_context": {
            "p_amp": p_amp,
            "r_dist": r_dist,
            "mass": mass_context,
            "gw": gw_context,
        },
        "thresholds": thresholds,
        "search_space_size": int(total),
        "pass_count_flavor": int(pass_count),
        "best_row": best,
        "verdict": verdict,
        "required_next_step": (
            "RUN_FINAL_STAGE_C_GATE_WITH_COMBINED_BRANCH"
            if pass_count > 0
            else "ESCALATE_TO_NONPERTURBATIVE_FLAVOR_BASIS_OR_CORE_REVISE"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2029: CKM BLOCKER SHARED FLAVOR BASIS SCAN",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- search_space_size: {total}",
        f"- pass_count_flavor: {pass_count}",
        "",
        "## Fixed Context",
        f"- p_amp/r_dist: {p_amp:.3f}/{r_dist:.3f}",
        f"- mass mean/max rel%: {mass_context['mean_rel_err_pct']:.3f}/{mass_context['max_rel_err_pct']:.3f}",
        f"- GW auc/adv/sep/gap: {gw_context['auc_h1l1_vs_ctrl']:.4f}/{gw_context['adv_shared_minus_ctrl_q90']:.4f}/{gw_context['sep_median_h1l1_minus_ctrl']:.6f}/{gw_context['control_median_gap']:.6f}",
        "",
        "## Best Flavor Row",
        f"- q_nu: {best['q_nu']}",
        f"- CKM/PMNS mean rel%: {best['flavor']['ckm_mean_rel_pct']:.3f}/{best['flavor']['pmns_mean_rel_pct']:.3f}",
        f"- all_pass: {best['all_pass']}",
        "",
        "## Required Next Step",
        f"- {out['required_next_step']}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2029] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2029] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2029] verdict={verdict} pass_count={pass_count}/{total}")


if __name__ == "__main__":
    main()
