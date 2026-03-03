#!/usr/bin/env python3
"""
QW-1944: Kernel-invariant flavor operator + GW check (no CKM/PMNS separate tuning).

Design:
- one deterministic map from kernel invariants -> operator params,
- one operator family used for CKM and PMNS,
- same shared params feed GW score synthesis.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1944_kernel_invariant_flavor_operator_and_gw.json"
OUT_MD = ROOT / "RAPORT_QW1944_KERNEL_INVARIANT_FLAVOR_OPERATOR_AND_GW.md"


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


def derive_kernel_invariant_params(kernel: Dict[str, float]) -> Dict[str, float]:
    d = np.arange(1.0, 13.0, dtype=float)
    k = np.abs(kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"]))
    w = k / max(np.sum(k), 1e-15)
    mean_d = float(np.sum(w * d))
    var_d = float(np.sum(w * (d - mean_d) ** 2))
    k1, k2 = float(k[0]), float(k[1])

    # Fixed deterministic mapping from invariants.
    p_amp = float(np.clip(0.25 + 0.30 * np.tanh((k1 - k2) / 0.20), 0.10, 1.20))
    r_dist = float(np.clip(-0.60 - 0.60 * np.tanh(var_d - 1.50), -1.50, 0.50))
    s_scale = float(np.clip(1.00 + 8.00 * abs(k1 - k2), 0.10, 20.00))

    return {
        "p_amp": p_amp,
        "r_dist": r_dist,
        "s_scale": s_scale,
        "mean_d": mean_d,
        "var_d": var_d,
        "abs_k1": k1,
        "abs_k2": k2,
    }


def sector_unitary(q: np.ndarray, params: Dict[str, float], kernel: Dict[str, float]) -> np.ndarray:
    q = np.asarray(q, dtype=float)
    d = 1.0 + cyclic_distance_matrix(q, q, modulus=24)
    n = len(q)

    a = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            amp = (abs(kernel_fn(np.array([d[i, j]]), kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])[0]) ** params["p_amp"]) * (
                d[i, j] ** params["r_dist"]
            )
            a[i, j] = (1.0 if i < j else -1.0) * params["s_scale"] * amp
    a = 0.5 * (a - a.T)
    u = expm(a)  # orthogonal/unitary from antisymmetric generator
    return u


def flavor_metrics(
    q_up: np.ndarray,
    q_down: np.ndarray,
    q_nu: np.ndarray,
    q_lep: np.ndarray,
    params: Dict[str, float],
    kernel: Dict[str, float],
) -> Dict[str, float]:
    u_up = sector_unitary(q_up, params=params, kernel=kernel)
    u_down = sector_unitary(q_down, params=params, kernel=kernel)
    u_nu = sector_unitary(q_nu, params=params, kernel=kernel)
    u_lep = sector_unitary(q_lep, params=params, kernel=kernel)

    ckm = np.abs(u_up.T @ u_down)
    pmns = np.abs(u_nu.T @ u_lep)

    ckm_rel = np.abs(ckm - CKM_REF) / np.clip(CKM_REF, 1e-12, None)
    pmns_rel = np.abs(pmns - PMNS_REF) / np.clip(PMNS_REF, 1e-12, None)
    return {
        "ckm_mean_rel_pct": float(100.0 * np.mean(ckm_rel)),
        "ckm_max_rel_pct": float(100.0 * np.max(ckm_rel)),
        "pmns_mean_rel_pct": float(100.0 * np.mean(pmns_rel)),
        "pmns_max_rel_pct": float(100.0 * np.max(pmns_rel)),
    }


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


def gw_metrics(params: Dict[str, float], kernel: Dict[str, float], gw_cache: Dict[str, np.ndarray]) -> Dict[str, float]:
    d = np.array([1.0, 2.0, 3.0, 4.0], dtype=float)
    raw_w = (np.abs(kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])) ** params["p_amp"]) * (
        d**params["r_dist"]
    )
    w = raw_w / max(np.sum(raw_w), 1e-15)

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
        "auc_h1l1_vs_ctrl": auc,
        "adv_shared_minus_ctrl_q90": adv,
        "sep_median_h1l1_minus_ctrl": sep,
        "control_median_gap": ctrl_gap,
    }


def scenario_eval(
    name: str,
    q_up: np.ndarray,
    q_down: np.ndarray,
    q_nu: np.ndarray,
    q_lep: np.ndarray,
    params: Dict[str, float],
    kernel: Dict[str, float],
    gw_cache: Dict[str, np.ndarray],
    thresholds: Dict[str, float],
) -> Dict[str, object]:
    f = flavor_metrics(q_up, q_down, q_nu, q_lep, params=params, kernel=kernel)
    g = gw_metrics(params=params, kernel=kernel, gw_cache=gw_cache)

    flags = {
        "ckm_mean_rel_pct_le_max": bool(f["ckm_mean_rel_pct"] <= thresholds["ckm_mean_rel_pct_max"]),
        "pmns_mean_rel_pct_le_max": bool(f["pmns_mean_rel_pct"] <= thresholds["pmns_mean_rel_pct_max"]),
        "gw_sep_ge_min": bool(g["sep_median_h1l1_minus_ctrl"] >= thresholds["gw_sep_min"]),
        "gw_adv_ge_min": bool(g["adv_shared_minus_ctrl_q90"] >= thresholds["gw_adv_min"]),
        "gw_auc_ge_min": bool(g["auc_h1l1_vs_ctrl"] >= thresholds["gw_auc_min"]),
        "gw_control_gap_le_max": bool(g["control_median_gap"] <= thresholds["gw_control_gap_max"]),
    }

    loss = (
        f["ckm_mean_rel_pct"] / thresholds["ckm_mean_rel_pct_max"]
        + f["pmns_mean_rel_pct"] / thresholds["pmns_mean_rel_pct_max"]
        + max(0.0, thresholds["gw_sep_min"] - g["sep_median_h1l1_minus_ctrl"]) / max(thresholds["gw_sep_min"], 1e-12)
        + max(0.0, thresholds["gw_adv_min"] - g["adv_shared_minus_ctrl_q90"]) / max(thresholds["gw_adv_min"], 1e-12)
        + max(0.0, thresholds["gw_auc_min"] - g["auc_h1l1_vs_ctrl"]) / max(thresholds["gw_auc_min"], 1e-12)
        + max(0.0, g["control_median_gap"] - thresholds["gw_control_gap_max"]) / max(thresholds["gw_control_gap_max"], 1e-12)
    )

    return {
        "name": name,
        "q_up": [int(x) for x in q_up],
        "q_down": [int(x) for x in q_down],
        "q_nu": [int(x) for x in q_nu],
        "q_lep": [int(x) for x in q_lep],
        "flavor": f,
        "gw": g,
        "flags": flags,
        "all_pass": bool(all(flags.values())),
        "loss": float(loss),
    }


def main() -> None:
    d1932 = json.loads((ROOT / "report_qw1932_physical_reparameterization_eta_scan.json").read_text(encoding="utf-8"))
    sel = d1932["selected"]
    kernel = {
        "omega": float(sel["fit"]["omega"]),
        "phi": float(sel["fit"]["phi"]),
        "beta": float(sel["fit"]["beta"]),
        "eta": float(sel["eta"]),
    }
    params = derive_kernel_invariant_params(kernel)

    thresholds = {
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

    scenarios = []

    # Strict baseline assignment.
    scenarios.append(
        scenario_eval(
            "baseline_assignment",
            np.array([0, 9, 14], dtype=float),
            np.array([7, 9, 14], dtype=float),
            np.array([0, 1, 2], dtype=float),
            np.array([24, 14, 9], dtype=float),
            params=params,
            kernel=kernel,
            gw_cache=gw_cache,
            thresholds=thresholds,
        )
    )

    # Optional audited candidates from QW-1943 (still one shared operator, no separate sector fit).
    p1943 = ROOT / "report_qw1943_topological_q_assignment_audit.json"
    if p1943.exists():
        d1943 = json.loads(p1943.read_text(encoding="utf-8"))
        for tag, key in [("audit_best_joint", "best_joint_row"), ("audit_best_flavor", "best_flavor_row")]:
            row = d1943.get(key)
            if row:
                scenarios.append(
                    scenario_eval(
                        tag,
                        np.array(row["q_up"], dtype=float),
                        np.array(row["q_down"], dtype=float),
                        np.array(row["q_nu"], dtype=float),
                        np.array(row["q_lep"], dtype=float),
                        params=params,
                        kernel=kernel,
                        gw_cache=gw_cache,
                        thresholds=thresholds,
                    )
                )

    best = min(scenarios, key=lambda s: s["loss"])
    baseline = next(s for s in scenarios if s["name"] == "baseline_assignment")

    if baseline["all_pass"]:
        verdict = "KERNEL_INVARIANT_OPERATOR_PASS_STRICT_BASELINE"
        selected_mode = "STRICT_BASELINE"
        required_next = "MERGE_WITH_HARD_MASS_FORMULA_AND_RUN_FINAL_GATE"
    elif best["all_pass"]:
        verdict = "KERNEL_INVARIANT_OPERATOR_PASS_AUDITED_ASSIGNMENT"
        selected_mode = "AUDITED_ASSIGNMENT"
        required_next = "RUN_COMPLEXITY_CONTROLLED_FINAL_GATE_WITH_AUDITED_ASSIGNMENT"
    else:
        verdict = "KERNEL_INVARIANT_OPERATOR_FAIL"
        selected_mode = "NONE"
        required_next = "REWORK_OPERATOR_STRUCTURE_AT_KERNEL_LEVEL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw1932_physical_reparameterization_eta_scan.json:selected",
        "kernel": kernel,
        "shared_params_from_kernel_invariants": params,
        "thresholds": thresholds,
        "scenarios": scenarios,
        "baseline": baseline,
        "best_scenario": best,
        "selected_mode": selected_mode,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1944: KERNEL-INVARIANT FLAVOR OPERATOR + GW",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Selected mode: **{selected_mode}**",
        "",
        "## Shared Params (Derived from Kernel Invariants)",
        f"- p_amp: {params['p_amp']:.6f}",
        f"- r_dist: {params['r_dist']:.6f}",
        f"- s_scale: {params['s_scale']:.6f}",
        "",
        "## Baseline Scenario",
        (
            f"- CKM/PMNS mean rel%: "
            f"{baseline['flavor']['ckm_mean_rel_pct']:.3f}/{baseline['flavor']['pmns_mean_rel_pct']:.3f}"
        ),
        (
            f"- GW auc/adv/sep/gap: "
            f"{baseline['gw']['auc_h1l1_vs_ctrl']:.4f}/{baseline['gw']['adv_shared_minus_ctrl_q90']:.4f}/"
            f"{baseline['gw']['sep_median_h1l1_minus_ctrl']:.6f}/{baseline['gw']['control_median_gap']:.6f}"
        ),
        f"- all_pass: {baseline['all_pass']}",
        "",
        "## Best Scenario",
        f"- name: {best['name']}",
        (
            f"- CKM/PMNS mean rel%: "
            f"{best['flavor']['ckm_mean_rel_pct']:.3f}/{best['flavor']['pmns_mean_rel_pct']:.3f}"
        ),
        (
            f"- GW auc/adv/sep/gap: "
            f"{best['gw']['auc_h1l1_vs_ctrl']:.4f}/{best['gw']['adv_shared_minus_ctrl_q90']:.4f}/"
            f"{best['gw']['sep_median_h1l1_minus_ctrl']:.6f}/{best['gw']['control_median_gap']:.6f}"
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

    print(f"[QW-1944] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1944] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1944] verdict={verdict} selected_mode={selected_mode}")


if __name__ == "__main__":
    main()

