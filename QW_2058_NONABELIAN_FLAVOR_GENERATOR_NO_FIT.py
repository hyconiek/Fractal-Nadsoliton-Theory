#!/usr/bin/env python3
"""
QW-2058: Nonabelian flavor generator gate (first-principles, no fit).

Strict protocol:
- frozen kernel from QW-2049,
- mass-chain inputs from QW-1961,
- deterministic nonabelian flavor generator in SU(3)-like basis,
- no optimization against CKM/PMNS and no sector retune.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2058_nonabelian_flavor_generator_no_fit.json"
OUT_MD = ROOT / "RAPORT_QW2058_NONABELIAN_FLAVOR_GENERATOR_NO_FIT.md"


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
    # Deterministic map from kernel invariants, reused for GW weights.
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


def gw_metrics(params: Dict[str, float], kernel: Dict[str, float]) -> Dict[str, float]:
    d = np.array([1.0, 2.0, 3.0, 4.0], dtype=float)
    raw_w = (np.abs(kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])) ** params["p_amp"]) * (
        d ** params["r_dist"]
    )
    w = raw_w / np.sum(raw_w)

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


def gell_mann_basis() -> Dict[str, np.ndarray]:
    l1 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex)
    l2 = np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], dtype=complex)
    l3 = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=complex)
    l4 = np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex)
    l5 = np.array([[0, 0, -1j], [0, 0, 0], [1j, 0, 0]], dtype=complex)
    l6 = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex)
    l7 = np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], dtype=complex)
    l8 = (1.0 / np.sqrt(3.0)) * np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]], dtype=complex)
    return {"l1": l1, "l2": l2, "l3": l3, "l4": l4, "l5": l5, "l6": l6, "l7": l7, "l8": l8}


def sector_generator(q_sector: np.ndarray, kernel: Dict[str, float]) -> Dict[str, object]:
    d = 1.0 + cyclic_distance_matrix(q_sector, q_sector, modulus=24)
    k = kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])

    i12 = float(np.abs(k[0, 1]) / (d[0, 1] ** (1.0 + abs(kernel["phi"]))))
    i13 = float(np.abs(k[0, 2]) / (d[0, 2] ** (1.0 + abs(kernel["phi"]))))
    i23 = float(np.abs(k[1, 2]) / (d[1, 2] ** (1.0 + abs(kernel["phi"]))))

    p12 = float(np.sign(k[0, 1]) * np.abs(np.sin(kernel["omega"] * d[0, 1] + kernel["phi"])))
    p13 = float(np.sign(k[0, 2]) * np.abs(np.sin(kernel["omega"] * d[0, 2] + kernel["phi"])))
    p23 = float(np.sign(k[1, 2]) * np.abs(np.sin(kernel["omega"] * d[1, 2] + kernel["phi"])))

    c_ord = float(np.tanh((q_sector[0] - 2.0 * q_sector[1] + q_sector[2]) / 24.0))
    c_mean = float((i12 + i13 + i23) / 3.0)

    # Deterministic coefficients in fixed nonabelian basis.
    a12 = i12 * (1.0 + 0.5 * c_ord)
    a13 = i13 * (1.0 - 0.5 * c_ord)
    a23 = i23 * (1.0 + 0.5 * np.sign(c_ord + 1e-12))

    b12 = p12 * i12 * (1.0 + abs(kernel["phi"]))
    b13 = p13 * i13 * (1.0 + abs(kernel["phi"]))
    b23 = p23 * i23 * (1.0 + abs(kernel["phi"]))

    d3 = c_ord * c_mean
    d8 = (i13 - i12) * np.tanh(kernel["eta"] - 1.0)

    basis = gell_mann_basis()
    g = (
        a12 * basis["l1"]
        + b12 * basis["l2"]
        + d3 * basis["l3"]
        + a13 * basis["l4"]
        + b13 * basis["l5"]
        + a23 * basis["l6"]
        + b23 * basis["l7"]
        + d8 * basis["l8"]
    )
    u = expm(1j * g)
    return {
        "invariants": {
            "i12": i12,
            "i13": i13,
            "i23": i23,
            "p12": p12,
            "p13": p13,
            "p23": p23,
            "c_ord": c_ord,
            "c_mean": c_mean,
        },
        "coefficients": {
            "a12": float(a12),
            "a13": float(a13),
            "a23": float(a23),
            "b12": float(b12),
            "b13": float(b13),
            "b23": float(b23),
            "d3": float(d3),
            "d8": float(d8),
        },
        "unitary": u,
    }


def flavor_nonabelian_metrics(
    q_up: np.ndarray,
    q_down: np.ndarray,
    q_nu: np.ndarray,
    q_lep: np.ndarray,
    kernel: Dict[str, float],
) -> Dict[str, object]:
    up = sector_generator(q_up, kernel)
    down = sector_generator(q_down, kernel)
    nu = sector_generator(q_nu, kernel)
    lep = sector_generator(q_lep, kernel)

    ckm = np.abs(up["unitary"].conj().T @ down["unitary"])
    pmns = np.abs(nu["unitary"].conj().T @ lep["unitary"])

    ckm_rel = np.abs(ckm - CKM_REF) / np.clip(CKM_REF, 1e-12, None)
    pmns_rel = np.abs(pmns - PMNS_REF) / np.clip(PMNS_REF, 1e-12, None)

    return {
        "ckm_pred_abs": ckm.tolist(),
        "pmns_pred_abs": pmns.tolist(),
        "ckm_mean_rel_pct": float(100.0 * np.mean(ckm_rel)),
        "pmns_mean_rel_pct": float(100.0 * np.mean(pmns_rel)),
        "sectors": {
            "up": {"invariants": up["invariants"], "coefficients": up["coefficients"]},
            "down": {"invariants": down["invariants"], "coefficients": down["coefficients"]},
            "nu": {"invariants": nu["invariants"], "coefficients": nu["coefficients"]},
            "lep": {"invariants": lep["invariants"], "coefficients": lep["coefficients"]},
        },
    }


def main() -> None:
    d2049 = json.loads((ROOT / "report_qw2049_spectral_micro_stagec_intersection_gate.json").read_text(encoding="utf-8"))
    kernel = {k: float(v) for k, v in d2049["stagec_pool"]["selected_kernel"].items()}

    r1961 = json.loads((ROOT / "report_qw1961_noncircular_gamma_q_derivation_matrix.json").read_text(encoding="utf-8"))
    best = r1961["summary"]["best_noncircular"]
    q_assign_name = str(best["q_assignment"])
    q_assign = r1961["inputs"]["q_assignments"][q_assign_name]
    gamma = float(best["gamma_value"])
    delta_info = float(r1961["info_split_source"]["delta_info"])

    params = derive_shared_params(kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])

    q_up = np.array([q_assign["Top"], q_assign["Charm"], q_assign["Muon"]], dtype=float)
    q_down = np.array([q_assign["Bottom"], q_assign["Tau"], q_assign["Muon"]], dtype=float)
    q_nu = np.array([0.0, 1.0, 2.0], dtype=float)
    q_lep = np.array([q_assign["Electron"], q_assign["Muon"], q_assign["Tau"]], dtype=float)

    m = mass_metrics(q_assign=q_assign, gamma=gamma, delta_info=delta_info)
    f = flavor_nonabelian_metrics(q_up=q_up, q_down=q_down, q_nu=q_nu, q_lep=q_lep, kernel=kernel)
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
        "protocol_no_fit": True,
        "protocol_no_sector_retune": True,
        "protocol_nonabelian_generator_explicit": True,
    }

    pass_count = int(sum(1 for v in flags.values() if v))
    total_flags = int(len(flags))
    all_pass = bool(pass_count == total_flags)

    verdict = "NONABELIAN_FIRST_PRINCIPLES_TRIAD_CLOSURE_PASS" if all_pass else "NONABELIAN_FIRST_PRINCIPLES_TRIAD_CLOSURE_FAIL"
    required_next = (
        "LOCK_NONABELIAN_GENERATOR_AND_RUN_INDEPENDENT_CONFIRMATORY_PACKAGE"
        if all_pass
        else "REVISE_FOUNDATIONAL_MAPPING_FROM_NADSOLITON_CHARACTERISTICS_TO_FLAVOR_GENERATOR"
    )

    gap_to_target = {
        "ckm_mean_rel_pct_over": float(max(0.0, f["ckm_mean_rel_pct"] - thresholds["ckm_mean_rel_pct_max"])),
        "pmns_mean_rel_pct_over": float(max(0.0, f["pmns_mean_rel_pct"] - thresholds["pmns_mean_rel_pct_max"])),
        "gw_control_gap_over": float(max(0.0, g["control_median_gap"] - thresholds["gw_control_gap_max"])),
    }

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw2049_spectral_micro_stagec_intersection_gate.json:stagec_pool.selected_kernel",
        "kernel": kernel,
        "mass_chain_source": {
            "file": "report_qw1961_noncircular_gamma_q_derivation_matrix.json",
            "q_assignment_name": q_assign_name,
            "gamma_value": gamma,
            "delta_info": delta_info,
        },
        "deterministic_shared_params_for_gw": params,
        "metrics": {"mass": m, "flavor_nonabelian": f, "gw": g},
        "thresholds": thresholds,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "all_pass": all_pass,
        "gap_to_target": gap_to_target,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2058: NONABELIAN FLAVOR GENERATOR (NO FIT)",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {pass_count}/{total_flags}",
        "",
        "## Frozen Inputs",
        (
            f"- kernel omega/phi/beta/eta: "
            f"{kernel['omega']:.6f}/{kernel['phi']:.6f}/{kernel['beta']:.6f}/{kernel['eta']:.6f}"
        ),
        f"- q_assignment (QW-1961): {q_assign_name}",
        f"- gamma: {gamma:.6f}",
        f"- delta_info: {delta_info:.6f}",
        "",
        "## Metrics",
        (
            f"- mass mean/max/tau-charm rel%: "
            f"{m['mean_rel_err_pct']:.3f}/{m['max_rel_err_pct']:.3f}/{m['tau_charm_ratio_rel_err_pct']:.3f}"
        ),
        (
            f"- flavor nonabelian CKM/PMNS mean rel%: "
            f"{f['ckm_mean_rel_pct']:.3f}/{f['pmns_mean_rel_pct']:.3f}"
        ),
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
            "## Gaps",
            f"- ckm_mean_rel_pct_over: {gap_to_target['ckm_mean_rel_pct_over']:.6f}",
            f"- pmns_mean_rel_pct_over: {gap_to_target['pmns_mean_rel_pct_over']:.6f}",
            f"- gw_control_gap_over: {gap_to_target['gw_control_gap_over']:.6f}",
            "",
            "## Required Next Step",
            f"- {required_next}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2058] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2058] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2058] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
