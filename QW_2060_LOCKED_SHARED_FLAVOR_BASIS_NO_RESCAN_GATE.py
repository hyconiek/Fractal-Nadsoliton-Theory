#!/usr/bin/env python3
"""
QW-2060: Locked shared flavor basis gate (no rescan).

Purpose:
- run a transparent no-rescan check with a locked shared flavor basis
  (historically strongest from QW-2029 family),
- evaluate on current frozen kernel QW-2049 and fixed mass-chain QW-1961,
- isolate whether blocker is flavor or GW-control-gap.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2060_locked_shared_flavor_basis_no_rescan_gate.json"
OUT_MD = ROOT / "RAPORT_QW2060_LOCKED_SHARED_FLAVOR_BASIS_NO_RESCAN_GATE.md"


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


def mass_metrics(q_assign: Dict[str, float], gamma: float, delta_info: float) -> Dict[str, object]:
    sector = {
        "Top": 0.0,
        "Bottom": +1.0,
        "Charm": +1.0,
        "Tau": -1.0,
        "Muon": -1.0,
        "Electron": -1.0,
    }
    rows: List[Dict[str, float]] = []
    errs: List[float] = []
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
        errs.append(float(err))
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
        + iso_tag * params.get("diag_iso", 0.0) * idx_centered
        + sector_tag * params.get("diag_sector", 0.0) * idx_centered
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
) -> Dict[str, object]:
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
        "ckm_pred_abs": ckm.tolist(),
        "pmns_pred_abs": pmns.tolist(),
    }


def gw_metrics(kernel: Dict[str, float], p_amp: float, r_dist: float) -> Dict[str, float]:
    d = np.array([1.0, 2.0, 3.0, 4.0], dtype=float)
    raw = (np.abs(kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])) ** p_amp) * (d**r_dist)
    w = raw / np.sum(raw)

    df = pd.read_csv(ROOT / "gw1831_window_features.csv")
    score = (
        w[0] * df["max_abs_corr"].to_numpy(dtype=float)
        + w[1] * df["mean_abs_corr"].to_numpy(dtype=float)
        + w[2] * df["corr_at_0ms"].to_numpy(dtype=float)
        + w[3] * df["corr_at_10ms"].to_numpy(dtype=float)
    )
    pair = df["pair"].astype(str).to_numpy()
    s_hl = score[pair == "H1-L1"]
    s_hv = score[pair == "H1-V1"]
    s_lv = score[pair == "L1-V1"]
    s_ctrl = np.concatenate([s_hv, s_lv])
    q90 = float(np.quantile(s_ctrl, 0.90))
    return {
        "auc_h1l1_vs_ctrl": float(rank_auc_pos_gt_neg(s_hl, s_ctrl)),
        "adv_shared_minus_ctrl_q90": float(np.mean(s_hl > q90) - np.mean(s_ctrl > q90)),
        "sep_median_h1l1_minus_ctrl": float(np.median(s_hl) - np.median(s_ctrl)),
        "control_median_gap": float(abs(np.median(s_hv) - np.median(s_lv))),
    }


def main() -> None:
    d2049 = json.loads((ROOT / "report_qw2049_spectral_micro_stagec_intersection_gate.json").read_text(encoding="utf-8"))
    kernel = {k: float(v) for k, v in d2049["stagec_pool"]["selected_kernel"].items()}

    r1961 = json.loads((ROOT / "report_qw1961_noncircular_gamma_q_derivation_matrix.json").read_text(encoding="utf-8"))
    q_name = str(r1961["summary"]["best_noncircular"]["q_assignment"])
    q_assign = r1961["inputs"]["q_assignments"][q_name]
    gamma = float(r1961["summary"]["best_noncircular"]["gamma_value"])
    delta_info = float(r1961["info_split_source"]["delta_info"])

    q_up = np.array([q_assign["Top"], q_assign["Charm"], q_assign["Muon"]], dtype=float)
    q_down = np.array([q_assign["Bottom"], q_assign["Tau"], q_assign["Muon"]], dtype=float)
    q_lep = np.array([q_assign["Electron"], q_assign["Muon"], q_assign["Tau"]], dtype=float)
    q_nu = np.array([2.0, 1.0, 0.0], dtype=float)

    # Locked structural basis (no scan, no optimization in this run).
    p_amp = 0.7
    r_dist = -0.2
    params = {
        "lambda_mix": 0.8,
        "rho_gap": 0.4,
        "chi_im": 0.3,
        "phase_q": 0.25,
        "phase_q3": 0.0,
        "theta_iso": 0.6,
        "theta_sector": 0.3,
        "diag_q_coeff": 0.1,
        "amp_qbias": 0.4,
        "diag_iso": 0.0,
        "diag_sector": 0.0,
    }

    m = mass_metrics(q_assign=q_assign, gamma=gamma, delta_info=delta_info)
    f = flavor_metrics(q_up, q_down, q_nu, q_lep, p_amp=p_amp, r_dist=r_dist, params=params, kernel=kernel)
    g = gw_metrics(kernel, p_amp=p_amp, r_dist=r_dist)

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
        "run_has_no_rescan": True,
        "run_has_no_refit": True,
        "strict_first_principles_from_kernel_only": False,
    }
    pass_count = int(sum(1 for v in flags.values() if v))
    total_flags = int(len(flags))
    all_pass = bool(all(flags.values()))

    verdict = "LOCKED_SHARED_FLAVOR_BASIS_NO_RESCAN_PASS" if all_pass else "LOCKED_SHARED_FLAVOR_BASIS_NO_RESCAN_PARTIAL"
    required_next = (
        "PROMOTE_TO_STRICT_FIRST_PRINCIPLES_BY_DERIVING_LOCKED_COEFFICIENTS_FROM_KERNEL_INVARIANTS"
        if not all_pass
        else "FREEZE_AND_EXTERNALIZE"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw2049_spectral_micro_stagec_intersection_gate.json:stagec_pool.selected_kernel",
        "kernel": kernel,
        "mass_chain_source": {
            "file": "report_qw1961_noncircular_gamma_q_derivation_matrix.json",
            "q_assignment_name": q_name,
            "gamma_value": gamma,
            "delta_info": delta_info,
        },
        "locked_flavor_basis": {
            "p_amp": p_amp,
            "r_dist": r_dist,
            "q_nu": [int(x) for x in q_nu],
            "params": params,
            "origin": "historically strongest shared-basis family (QW-2029), fixed without rescan",
        },
        "metrics": {"mass": m, "flavor": f, "gw": g},
        "thresholds": thresholds,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "all_pass": all_pass,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2060: LOCKED SHARED FLAVOR BASIS NO-RESCAN GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {pass_count}/{total_flags}",
        "",
        "## Locked Basis",
        f"- p_amp/r_dist: {p_amp:.3f}/{r_dist:.3f}",
        f"- q_nu: {[int(x) for x in q_nu]}",
        "",
        "## Metrics",
        (
            f"- mass mean/max/tau-charm rel%: "
            f"{m['mean_rel_err_pct']:.3f}/{m['max_rel_err_pct']:.3f}/{m['tau_charm_ratio_rel_err_pct']:.3f}"
        ),
        f"- flavor CKM/PMNS mean rel%: {f['ckm_mean_rel_pct']:.3f}/{f['pmns_mean_rel_pct']:.3f}",
        (
            f"- GW auc/adv/sep/gap: "
            f"{g['auc_h1l1_vs_ctrl']:.4f}/{g['adv_shared_minus_ctrl_q90']:.4f}/"
            f"{g['sep_median_h1l1_minus_ctrl']:.6f}/{g['control_median_gap']:.6f}"
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

    print(f"[QW-2060] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2060] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2060] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
