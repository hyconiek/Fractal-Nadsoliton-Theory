#!/usr/bin/env python3
"""
QW-2056: First-principles flavor operator family frontier (no fit).

Goal:
- test a finite set of deterministic, theoretically-motivated operator families,
- keep one frozen kernel (QW-2049),
- keep no fitting / no sector retuning,
- quantify whether CKM/PMNS closure is reachable in this family space.
"""

from __future__ import annotations

import itertools
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy.linalg import polar


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2056_first_principles_flavor_operator_family_frontier.json"
OUT_MD = ROOT / "RAPORT_QW2056_FIRST_PRINCIPLES_FLAVOR_OPERATOR_FAMILY_FRONTIER.md"


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


def flavor_prediction_abs(
    q_left: np.ndarray,
    q_right: np.ndarray,
    kernel: Dict[str, float],
    p_amp: float,
    r_dist: float,
    phase_scale: float,
    amp_mode: str,
    phase_mode: str,
) -> np.ndarray:
    n = len(q_left)
    d = 1.0 + cyclic_distance_matrix(q_left, q_right, modulus=24)
    k = kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])

    if amp_mode == "signed_pow":
        amp = np.sign(k) * ((np.abs(k) ** p_amp) * (d**r_dist))
    elif amp_mode == "diag_bias":
        amp = np.sign(k) * (np.abs(k) ** (1.2 * p_amp)) / (d ** abs(r_dist) + 1e-12)
    elif amp_mode == "local_exp":
        length = max(1e-9, 1.0 + 5.0 * float(np.mean(np.abs(k))))
        amp = np.sign(k) * np.exp(-(d - 1.0) / length) * (np.abs(k) ** (0.5 + 0.5 * p_amp))
    elif amp_mode == "kernel_only":
        amp = k
    else:
        raise ValueError(f"Unknown amp_mode: {amp_mode}")

    idx = np.arange(n, dtype=float)
    gap = np.abs(idx[:, None] - idx[None, :])
    sign = np.where((idx[:, None] - idx[None, :]) < 0.0, 1.0, -1.0)

    if phase_mode == "sgn_gap":
        phase = np.exp(1j * (kernel["phi"] + phase_scale * kernel["omega"] * sign * gap))
    elif phase_mode == "sgn_gap_pi":
        phase = np.exp(1j * (kernel["phi"] + phase_scale * np.pi * sign * gap / (max(n - 1, 1))))
    elif phase_mode == "diag_locked":
        phase = np.exp(1j * (kernel["phi"] + phase_scale * kernel["omega"] * (sign * gap + (idx[:, None] == idx[None, :]) * 0.5)))
    elif phase_mode == "kernel_sign":
        phase = np.exp(1j * (kernel["phi"] + phase_scale * kernel["omega"] * sign * np.sign(k) * gap))
    else:
        raise ValueError(f"Unknown phase_mode: {phase_mode}")

    u = polar(amp * phase)[0]
    return np.abs(u)


def rel_mean_pct(pred: np.ndarray, ref: np.ndarray) -> float:
    return float(100.0 * np.mean(np.abs(pred - ref) / np.clip(ref, 1e-12, None)))


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


def gw_metrics(kernel: Dict[str, float], p_amp: float, r_dist: float, gw_mode: str) -> Dict[str, float]:
    d = np.array([1.0, 2.0, 3.0, 4.0], dtype=float)
    k = np.abs(kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"]))
    if gw_mode == "std":
        raw_w = (k**p_amp) * (d**r_dist)
    elif gw_mode == "local":
        raw_w = (k ** (0.7 + 0.3 * p_amp)) / (d ** abs(r_dist) + 1e-12)
    elif gw_mode == "flat":
        raw_w = k + 1e-12
    else:
        raise ValueError(f"Unknown gw_mode: {gw_mode}")
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
    return {
        "auc_h1l1_vs_ctrl": float(rank_auc_pos_gt_neg(s_hl, s_ctrl)),
        "adv_shared_minus_ctrl_q90": float(np.mean(s_hl > q90) - np.mean(s_ctrl > q90)),
        "sep_median_h1l1_minus_ctrl": float(np.median(s_hl) - np.median(s_ctrl)),
        "control_median_gap": float(abs(np.median(s_hv) - np.median(s_lv))),
    }


def derive_family_params(kernel: Dict[str, float], family_name: str) -> Dict[str, float]:
    d = np.arange(1.0, 13.0, dtype=float)
    k = np.abs(kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"]))
    w = k / np.sum(k)
    mean_d = float(np.sum(w * d))
    var_d = float(np.sum(w * (d - mean_d) ** 2))
    decay = float(k[3] / max(k[0], 1e-15))

    if family_name == "legacy":
        return {
            "p_amp": float(1.0 + 0.60 * np.tanh((mean_d - 2.0) / 2.0)),
            "r_dist": float(np.tanh((var_d - 1.0) / 2.5)),
            "phase_scale": float(1.0 + 0.70 * np.tanh(abs(kernel["phi"])) + 0.25 * np.tanh(1.0 - decay)),
            "amp_mode": "signed_pow",
            "phase_mode": "sgn_gap",
            "gw_mode": "std",
        }
    if family_name == "locality":
        return {
            "p_amp": float(1.0 + 0.60 * np.tanh((mean_d - 2.0) / 2.0)),
            "r_dist": float(-abs(np.tanh((var_d - 1.0) / 2.5))),
            "phase_scale": float(1.0 + 0.50 * np.tanh(abs(kernel["phi"]))),
            "amp_mode": "diag_bias",
            "phase_mode": "sgn_gap_pi",
            "gw_mode": "local",
        }
    if family_name == "critical":
        return {
            "p_amp": 1.0,
            "r_dist": 0.0,
            "phase_scale": float(1.0 + 0.60 * np.tanh(abs(kernel["phi"]))),
            "amp_mode": "kernel_only",
            "phase_mode": "diag_locked",
            "gw_mode": "flat",
        }
    if family_name == "phase_sign":
        return {
            "p_amp": 1.2,
            "r_dist": -0.35,
            "phase_scale": float(1.0 + 0.80 * np.tanh(abs(kernel["phi"]))),
            "amp_mode": "local_exp",
            "phase_mode": "kernel_sign",
            "gw_mode": "local",
        }
    if family_name == "ultra_local":
        return {
            "p_amp": 1.35,
            "r_dist": -0.80,
            "phase_scale": float(1.0 + 0.45 * np.tanh(abs(kernel["phi"]))),
            "amp_mode": "diag_bias",
            "phase_mode": "sgn_gap_pi",
            "gw_mode": "local",
        }
    raise ValueError(f"Unknown family_name: {family_name}")


def qeff_from_mass(m_mev: float, gamma: float) -> float:
    return float(-(4.0 / gamma) * (np.log(m_mev / 173_000.0) / np.log(4.0)))


def main() -> None:
    d2049 = json.loads((ROOT / "report_qw2049_spectral_micro_stagec_intersection_gate.json").read_text(encoding="utf-8"))
    kernel = {k: float(v) for k, v in d2049["stagec_pool"]["selected_kernel"].items()}

    r1961 = json.loads((ROOT / "report_qw1961_noncircular_gamma_q_derivation_matrix.json").read_text(encoding="utf-8"))
    best = r1961["summary"]["best_noncircular"]
    q_assign_name = str(best["q_assignment"])
    q_assign = r1961["inputs"]["q_assignments"][q_assign_name]
    gamma = float(best["gamma_value"])
    delta_info = float(r1961["info_split_source"]["delta_info"])

    q_nu = np.array([0.0, 1.0, 2.0], dtype=float)
    q_lep = np.array([q_assign["Electron"], q_assign["Muon"], q_assign["Tau"]], dtype=float)

    q_proxy_up = np.array([q_assign["Top"], q_assign["Charm"], q_assign["Muon"]], dtype=float)
    q_proxy_down = np.array([q_assign["Bottom"], q_assign["Tau"], q_assign["Muon"]], dtype=float)

    # Quark-only diagnostic mapping from mass law inversion.
    mass_quark_mev = {"u": 2.16, "d": 4.67, "s": 93.4, "c": 1270.0, "b": 4180.0, "t": 173000.0}
    q_u = qeff_from_mass(mass_quark_mev["u"], gamma) - delta_info
    q_c = qeff_from_mass(mass_quark_mev["c"], gamma) - delta_info
    q_t = qeff_from_mass(mass_quark_mev["t"], gamma) - delta_info
    q_d = qeff_from_mass(mass_quark_mev["d"], gamma) - delta_info
    q_s = qeff_from_mass(mass_quark_mev["s"], gamma) - delta_info
    q_b = qeff_from_mass(mass_quark_mev["b"], gamma) - delta_info
    q_quark_up = np.array([q_u, q_c, q_t], dtype=float)
    q_quark_down = np.array([q_d, q_s, q_b], dtype=float)

    q_schemes = {
        "proxy_old": {"q_up": q_proxy_up, "q_down": q_proxy_down},
        "quark_mass_inversion": {"q_up": q_quark_up, "q_down": q_quark_down},
    }

    families = ["legacy", "locality", "critical", "phase_sign", "ultra_local"]

    thresholds = {
        "ckm_mean_rel_pct_max": 15.0,
        "pmns_mean_rel_pct_max": 15.0,
        "gw_sep_min": 0.0020,
        "gw_adv_min": 0.30,
        "gw_auc_min": 0.75,
        "gw_control_gap_max": 0.0025,
    }

    rows: List[Dict[str, object]] = []
    for family_name in families:
        fp = derive_family_params(kernel, family_name)
        gw = gw_metrics(kernel, fp["p_amp"], fp["r_dist"], fp["gw_mode"])
        pmns_pred = flavor_prediction_abs(
            q_nu,
            q_lep,
            kernel=kernel,
            p_amp=fp["p_amp"],
            r_dist=fp["r_dist"],
            phase_scale=fp["phase_scale"],
            amp_mode=fp["amp_mode"],
            phase_mode=fp["phase_mode"],
        )
        pmns_mean = rel_mean_pct(pmns_pred, PMNS_REF)

        for scheme_name, scheme in q_schemes.items():
            ckm_pred = flavor_prediction_abs(
                scheme["q_up"],
                scheme["q_down"],
                kernel=kernel,
                p_amp=fp["p_amp"],
                r_dist=fp["r_dist"],
                phase_scale=fp["phase_scale"],
                amp_mode=fp["amp_mode"],
                phase_mode=fp["phase_mode"],
            )
            ckm_mean = rel_mean_pct(ckm_pred, CKM_REF)

            best_perm_err = float("inf")
            best_perm = None
            for pu in itertools.permutations(range(3)):
                for pd in itertools.permutations(range(3)):
                    ckm_perm = flavor_prediction_abs(
                        scheme["q_up"][list(pu)],
                        scheme["q_down"][list(pd)],
                        kernel=kernel,
                        p_amp=fp["p_amp"],
                        r_dist=fp["r_dist"],
                        phase_scale=fp["phase_scale"],
                        amp_mode=fp["amp_mode"],
                        phase_mode=fp["phase_mode"],
                    )
                    err = rel_mean_pct(ckm_perm, CKM_REF)
                    if err < best_perm_err:
                        best_perm_err = float(err)
                        best_perm = {"perm_up": list(pu), "perm_down": list(pd)}

            flags = {
                "ckm_mean_rel_pct_le_max": bool(ckm_mean <= thresholds["ckm_mean_rel_pct_max"]),
                "pmns_mean_rel_pct_le_max": bool(pmns_mean <= thresholds["pmns_mean_rel_pct_max"]),
                "gw_sep_ge_min": bool(gw["sep_median_h1l1_minus_ctrl"] >= thresholds["gw_sep_min"]),
                "gw_adv_ge_min": bool(gw["adv_shared_minus_ctrl_q90"] >= thresholds["gw_adv_min"]),
                "gw_auc_ge_min": bool(gw["auc_h1l1_vs_ctrl"] >= thresholds["gw_auc_min"]),
                "gw_control_gap_le_max": bool(gw["control_median_gap"] <= thresholds["gw_control_gap_max"]),
            }
            pass_count = int(sum(1 for v in flags.values() if v))
            total_flags = int(len(flags))

            rows.append(
                {
                    "family": family_name,
                    "q_scheme": scheme_name,
                    "params": fp,
                    "metrics": {
                        "ckm_mean_rel_pct_fixed_order": ckm_mean,
                        "pmns_mean_rel_pct": pmns_mean,
                        "ckm_mean_rel_pct_perm_envelope": best_perm_err,
                        "perm_envelope_argmin": best_perm,
                        "gw": gw,
                    },
                    "flags": flags,
                    "pass_count": pass_count,
                    "total_flags": total_flags,
                    "all_pass": bool(pass_count == total_flags),
                }
            )

    def closure_score(row: Dict[str, object]) -> float:
        m = row["metrics"]
        g = m["gw"]
        return float(
            max(0.0, m["ckm_mean_rel_pct_fixed_order"] - thresholds["ckm_mean_rel_pct_max"])
            + max(0.0, m["pmns_mean_rel_pct"] - thresholds["pmns_mean_rel_pct_max"])
            + 1_000.0 * max(0.0, thresholds["gw_sep_min"] - g["sep_median_h1l1_minus_ctrl"])
            + 500.0 * max(0.0, thresholds["gw_adv_min"] - g["adv_shared_minus_ctrl_q90"])
            + 500.0 * max(0.0, thresholds["gw_auc_min"] - g["auc_h1l1_vs_ctrl"])
            + 100_000.0 * max(0.0, g["control_median_gap"] - thresholds["gw_control_gap_max"])
        )

    rows_sorted = sorted(rows, key=closure_score)
    any_all_pass = any(bool(r["all_pass"]) for r in rows_sorted)
    best_row = rows_sorted[0]
    best_flavor_fixed = min(rows, key=lambda r: float(r["metrics"]["ckm_mean_rel_pct_fixed_order"] + r["metrics"]["pmns_mean_rel_pct"]))
    best_ckm_fixed = min(rows, key=lambda r: float(r["metrics"]["ckm_mean_rel_pct_fixed_order"]))
    best_ckm_perm_env = min(rows, key=lambda r: float(r["metrics"]["ckm_mean_rel_pct_perm_envelope"]))

    verdict = (
        "FIRST_PRINCIPLES_OPERATOR_FAMILY_HAS_STRICT_CLOSURE_CANDIDATE"
        if any_all_pass
        else "FIRST_PRINCIPLES_OPERATOR_FAMILY_FRONTIER_FAILS_TO_CLOSE_FLAVOR"
    )
    required_next = (
        "FREEZE_BEST_STRICT_FAMILY_AND_RUN_FULL_TRIAD_GATE"
        if any_all_pass
        else "DERIVE_NONABELIAN_FLAVOR_GENERATOR_FROM_KERNEL_DYNAMICS_WITHOUT_NEW_FREE_PARAMETERS"
    )

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
        "thresholds": thresholds,
        "rows_sorted": rows_sorted,
        "best_row": best_row,
        "best_flavor_fixed": best_flavor_fixed,
        "best_ckm_fixed": best_ckm_fixed,
        "best_ckm_perm_envelope": best_ckm_perm_env,
        "any_all_pass": any_all_pass,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2056: FIRST-PRINCIPLES FLAVOR OPERATOR FAMILY FRONTIER",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- any_all_pass: {any_all_pass}",
        "",
        "## Best Row",
        f"- family: {best_row['family']}",
        f"- q_scheme: {best_row['q_scheme']}",
        f"- pass_count: {best_row['pass_count']}/{best_row['total_flags']}",
        (
            f"- CKM fixed / PMNS / CKM perm-envelope: "
            f"{best_row['metrics']['ckm_mean_rel_pct_fixed_order']:.3f} / "
            f"{best_row['metrics']['pmns_mean_rel_pct']:.3f} / "
            f"{best_row['metrics']['ckm_mean_rel_pct_perm_envelope']:.3f}"
        ),
        (
            f"- GW auc/adv/sep/gap: "
            f"{best_row['metrics']['gw']['auc_h1l1_vs_ctrl']:.4f}/"
            f"{best_row['metrics']['gw']['adv_shared_minus_ctrl_q90']:.4f}/"
            f"{best_row['metrics']['gw']['sep_median_h1l1_minus_ctrl']:.6f}/"
            f"{best_row['metrics']['gw']['control_median_gap']:.6f}"
        ),
        "",
        "## Best Flavor Fixed (CKM+PMNS)",
        (
            f"- family/q_scheme: {best_flavor_fixed['family']} / {best_flavor_fixed['q_scheme']} | "
            f"CKM={best_flavor_fixed['metrics']['ckm_mean_rel_pct_fixed_order']:.3f} | "
            f"PMNS={best_flavor_fixed['metrics']['pmns_mean_rel_pct']:.3f}"
        ),
        "",
        "## Best CKM (Fixed Order)",
        (
            f"- family/q_scheme: {best_ckm_fixed['family']} / {best_ckm_fixed['q_scheme']} | "
            f"CKM={best_ckm_fixed['metrics']['ckm_mean_rel_pct_fixed_order']:.3f} | "
            f"PMNS={best_ckm_fixed['metrics']['pmns_mean_rel_pct']:.3f}"
        ),
        "",
        "## Best CKM (Permutation Envelope, Diagnostic)",
        (
            f"- family/q_scheme: {best_ckm_perm_env['family']} / {best_ckm_perm_env['q_scheme']} | "
            f"CKM_perm_env={best_ckm_perm_env['metrics']['ckm_mean_rel_pct_perm_envelope']:.3f} | "
            f"PMNS={best_ckm_perm_env['metrics']['pmns_mean_rel_pct']:.3f}"
        ),
        "",
        "## Top 10 Rows",
    ]

    top = rows_sorted[:10]
    for i, row in enumerate(top, start=1):
        lines.append(
            (
                f"- {i}. {row['family']} | {row['q_scheme']} | pass {row['pass_count']}/{row['total_flags']} | "
                f"CKM={row['metrics']['ckm_mean_rel_pct_fixed_order']:.3f} | "
                f"PMNS={row['metrics']['pmns_mean_rel_pct']:.3f} | "
                f"CKM_perm_env={row['metrics']['ckm_mean_rel_pct_perm_envelope']:.3f}"
            )
        )

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

    print(f"[QW-2056] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2056] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2056] verdict={verdict} any_all_pass={any_all_pass}")


if __name__ == "__main__":
    main()
