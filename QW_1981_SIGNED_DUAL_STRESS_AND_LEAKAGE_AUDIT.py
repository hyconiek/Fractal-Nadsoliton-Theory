#!/usr/bin/env python3
"""
QW-1981: Stress and leakage audit for QW-1980 signed-dual strict candidate.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
IN_QW1970 = ROOT / "report_qw1970_structural_gw_control_term_gate.json"
IN_QW1980 = ROOT / "report_qw1980_signed_dual_basis_independent_gate.json"
OUT_JSON = ROOT / "report_qw1981_signed_dual_stress_and_leakage_audit.json"
OUT_MD = ROOT / "RAPORT_QW1981_SIGNED_DUAL_STRESS_AND_LEAKAGE_AUDIT.md"


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


def gw_metrics(s_hl: np.ndarray, s_hv: np.ndarray, s_lv: np.ndarray) -> Dict[str, float]:
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


def gw_flags(g: Dict[str, float], thr: Dict[str, float]) -> Dict[str, bool]:
    return {
        "gw_sep_ge_min": bool(g["sep_median_h1l1_minus_ctrl"] >= thr["gw_sep_min"]),
        "gw_adv_ge_min": bool(g["adv_shared_minus_ctrl_q90"] >= thr["gw_adv_min"]),
        "gw_auc_ge_min": bool(g["auc_h1l1_vs_ctrl"] >= thr["gw_auc_min"]),
        "gw_control_gap_le_max": bool(g["control_median_gap"] <= thr["gw_control_gap_max"]),
    }


def bootstrap_pass_rate(
    s_hl: np.ndarray,
    s_hv: np.ndarray,
    s_lv: np.ndarray,
    thr: Dict[str, float],
    n_boot: int,
    seed: int,
) -> float:
    rng = np.random.default_rng(seed)
    n_hl, n_hv, n_lv = len(s_hl), len(s_hv), len(s_lv)
    pass_count = 0
    for _ in range(n_boot):
        b_hl = s_hl[rng.integers(0, n_hl, size=n_hl, endpoint=False)]
        b_hv = s_hv[rng.integers(0, n_hv, size=n_hv, endpoint=False)]
        b_lv = s_lv[rng.integers(0, n_lv, size=n_lv, endpoint=False)]
        if all(gw_flags(gw_metrics(b_hl, b_hv, b_lv), thr).values()):
            pass_count += 1
    return float(pass_count / n_boot)


def block_resample(arr: np.ndarray, block_len: int, rng: np.random.Generator) -> np.ndarray:
    n = len(arr)
    idx = []
    while len(idx) < n:
        start = int(rng.integers(0, n))
        idx.extend([(start + k) % n for k in range(block_len)])
    return arr[np.array(idx[:n], dtype=int)]


def bootstrap_block_pass_rate(
    s_hl: np.ndarray,
    s_hv: np.ndarray,
    s_lv: np.ndarray,
    thr: Dict[str, float],
    n_boot: int,
    block_len: int,
    seed: int,
) -> float:
    rng = np.random.default_rng(seed)
    pass_count = 0
    for _ in range(n_boot):
        b_hl = block_resample(s_hl, block_len, rng)
        b_hv = block_resample(s_hv, block_len, rng)
        b_lv = block_resample(s_lv, block_len, rng)
        if all(gw_flags(gw_metrics(b_hl, b_hv, b_lv), thr).values()):
            pass_count += 1
    return float(pass_count / n_boot)


def build_scores(
    xi1: float,
    xi3: float,
    base_score: np.ndarray,
    c1: np.ndarray,
    c3: np.ndarray,
    pairs: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    s = base_score + xi1 * c1 + xi3 * c3
    return s[pairs == 0], s[pairs == 1], s[pairs == 2]


def verdict(real_iid_min: float, real_block_min: float, null_p90: float, adv_null_mean: float) -> Tuple[str, str]:
    if real_iid_min >= 0.95 and real_block_min >= 0.90 and null_p90 <= 0.25 and adv_null_mean <= 0.40:
        return (
            "SIGNED_DUAL_STRESS_PASS_STRONG",
            "PROMOTE_TO_EXTERNAL_CONFIRMATORY_PREP",
        )
    if real_iid_min >= 0.90 and real_block_min >= 0.85 and null_p90 <= 0.35:
        return (
            "SIGNED_DUAL_STRESS_PASS_PARTIAL",
            "RUN_FINAL_INDEPENDENT_NULL_REPLICATION_BEFORE_EXTERNAL",
        )
    return (
        "SIGNED_DUAL_STRESS_FAIL_OR_LEAKAGE_RISK",
        "DO_NOT_FREEZE_YET_AND_HARDEN_NULL_PROTOCOL",
    )


def main() -> None:
    r1970 = json.loads(IN_QW1970.read_text(encoding="utf-8"))
    r1980 = json.loads(IN_QW1980.read_text(encoding="utf-8"))
    thr = json.loads((ROOT / "report_qw1969_bootstrap_robust_recenter_search.json").read_text(encoding="utf-8"))["thresholds"]

    best = r1980["best"]
    xi1 = float(best["xi1"])
    xi3 = float(best["xi3"])

    kernel = r1970["fixed_components"]["kernel"]
    params = r1970["fixed_components"]["params"]

    df = pd.read_csv(ROOT / "gw1831_window_features.csv")
    pair_map = {"H1-L1": 0, "H1-V1": 1, "L1-V1": 2}
    pairs = df["pair"].map(pair_map).to_numpy(dtype=int)
    feats = df[["max_abs_corr", "mean_abs_corr", "corr_at_0ms", "corr_at_10ms"]].to_numpy(dtype=float)

    d = np.array([1.0, 2.0, 3.0, 4.0], dtype=float)
    kvec = np.cos(kernel["omega"] * d + kernel["phi"]) / (1.0 + kernel["beta"] * (d**kernel["eta"]))
    w_raw = (np.abs(kvec) ** params["p_amp"]) * (d**params["r_dist"])
    w = w_raw / np.sum(w_raw)
    base_score = feats @ w

    lag_s = df["best_lag_ms"].to_numpy(dtype=float) * 1e-3
    lag_phase_sin = np.sin(kernel["omega"] * lag_s + kernel["phi"])
    lag_phase_cos = np.cos(kernel["omega"] * lag_s + kernel["phi"])
    corr0 = df["corr_at_0ms"].to_numpy(dtype=float)
    corr10 = df["corr_at_10ms"].to_numpy(dtype=float)
    mean_abs = df["mean_abs_corr"].to_numpy(dtype=float)

    pair_sign = np.where(pairs == 1, 1.0, np.where(pairs == 2, -1.0, 0.0))
    ctrl_idx = np.where(pairs != 0)[0]
    n_ctrl = len(ctrl_idx)
    n_plus = int(np.sum(pairs == 1))

    c1_raw = pair_sign * lag_phase_sin * (corr0 - corr10)
    c3_raw = pair_sign * lag_phase_cos * (mean_abs + np.abs(corr0) + np.abs(corr10))
    c1 = c1_raw / max(float(np.std(c1_raw)), 1e-12)
    c3 = c3_raw / max(float(np.std(c3_raw)), 1e-12)

    s_hl, s_hv, s_lv = build_scores(xi1, xi3, base_score, c1, c3, pairs)
    g_det = gw_metrics(s_hl, s_hv, s_lv)
    det_flags = gw_flags(g_det, thr)

    # Real stability
    iid_rates = [bootstrap_pass_rate(s_hl, s_hv, s_lv, thr, n_boot=2000, seed=110000 + i) for i in range(8)]
    block_rates = [bootstrap_block_pass_rate(s_hl, s_hv, s_lv, thr, n_boot=1200, block_len=8, seed=111000 + i) for i in range(6)]

    # Random null topology
    rng_null = np.random.default_rng(112000)
    null_rates = []
    null_det_pass = 0
    for i in range(120):
        signs = np.array([1.0] * n_plus + [-1.0] * (n_ctrl - n_plus), dtype=float)
        rng_null.shuffle(signs)
        rand_sign = np.zeros_like(pair_sign)
        rand_sign[ctrl_idx] = signs
        c1n_raw = rand_sign * lag_phase_sin * (corr0 - corr10)
        c3n_raw = rand_sign * lag_phase_cos * (mean_abs + np.abs(corr0) + np.abs(corr10))
        c1n = c1n_raw / max(float(np.std(c1n_raw)), 1e-12)
        c3n = c3n_raw / max(float(np.std(c3n_raw)), 1e-12)
        nhl, nhv, nlv = build_scores(xi1, xi3, base_score, c1n, c3n, pairs)
        g_n = gw_metrics(nhl, nhv, nlv)
        if all(gw_flags(g_n, thr).values()):
            null_det_pass += 1
        rb = bootstrap_pass_rate(nhl, nhv, nlv, thr, n_boot=300, seed=113000 + i)
        null_rates.append(float(rb))

    null_mean = float(np.mean(null_rates))
    null_p90 = float(np.quantile(np.array(null_rates, dtype=float), 0.90))

    # Adversarial null: choose best sign among K candidates, then bootstrap.
    rng_adv = np.random.default_rng(114000)
    adv_rates = []
    for i in range(30):
        best_candidate = None
        for _ in range(30):
            signs = np.array([1.0] * n_plus + [-1.0] * (n_ctrl - n_plus), dtype=float)
            rng_adv.shuffle(signs)
            rand_sign = np.zeros_like(pair_sign)
            rand_sign[ctrl_idx] = signs
            c1n_raw = rand_sign * lag_phase_sin * (corr0 - corr10)
            c3n_raw = rand_sign * lag_phase_cos * (mean_abs + np.abs(corr0) + np.abs(corr10))
            c1n = c1n_raw / max(float(np.std(c1n_raw)), 1e-12)
            c3n = c3n_raw / max(float(np.std(c3n_raw)), 1e-12)
            nhl, nhv, nlv = build_scores(xi1, xi3, base_score, c1n, c3n, pairs)
            g_n = gw_metrics(nhl, nhv, nlv)
            det_score = (
                1.2 * (g_n["auc_h1l1_vs_ctrl"] - thr["gw_auc_min"])
                + 1.0 * (g_n["adv_shared_minus_ctrl_q90"] - thr["gw_adv_min"])
                + 120.0 * (g_n["sep_median_h1l1_minus_ctrl"] - thr["gw_sep_min"])
                + 180.0 * (thr["gw_control_gap_max"] - g_n["control_median_gap"])
            )
            if (best_candidate is None) or (det_score > best_candidate["det_score"]):
                best_candidate = {"nhl": nhl, "nhv": nhv, "nlv": nlv, "det_score": det_score}
        rb = bootstrap_pass_rate(best_candidate["nhl"], best_candidate["nhv"], best_candidate["nlv"], thr, n_boot=300, seed=115000 + i)
        adv_rates.append(float(rb))

    adv_mean = float(np.mean(adv_rates))
    adv_p90 = float(np.quantile(np.array(adv_rates, dtype=float), 0.90))

    real_iid_min = float(np.min(iid_rates))
    real_block_min = float(np.min(block_rates))
    v, nxt = verdict(real_iid_min, real_block_min, null_p90, adv_mean)

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": [IN_QW1970.name, IN_QW1980.name],
        "candidate": {
            "xi1": xi1,
            "xi3": xi3,
            "deterministic_gw": g_det,
            "deterministic_flags": det_flags,
        },
        "real_stability": {
            "iid_bootstrap_n": 2000,
            "iid_rates_8seeds": iid_rates,
            "iid_min": real_iid_min,
            "iid_median": float(np.median(iid_rates)),
            "iid_max": float(np.max(iid_rates)),
            "block_bootstrap_n": 1200,
            "block_len": 8,
            "block_rates_6seeds": block_rates,
            "block_min": real_block_min,
            "block_median": float(np.median(block_rates)),
            "block_max": float(np.max(block_rates)),
        },
        "null_random_topology": {
            "n_topologies": 120,
            "bootstrap_n_each": 300,
            "deterministic_pass_count": int(null_det_pass),
            "deterministic_pass_rate": float(null_det_pass / 120.0),
            "bootstrap_mean": null_mean,
            "bootstrap_p90": null_p90,
            "bootstrap_std": float(np.std(null_rates)),
        },
        "null_adversarial_topology": {
            "n_groups": 30,
            "k_candidates_per_group": 30,
            "bootstrap_n_each": 300,
            "bootstrap_mean": adv_mean,
            "bootstrap_p90": adv_p90,
            "bootstrap_std": float(np.std(adv_rates)),
        },
        "verdict": v,
        "required_next_step": nxt,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1981: SIGNED DUAL STRESS AND LEAKAGE AUDIT",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{v}**",
        "",
        "## Candidate",
        f"- xi1/xi3: {xi1:.6f}/{xi3:.6f}",
        (
            f"- deterministic GW auc/adv/sep/gap: "
            f"{g_det['auc_h1l1_vs_ctrl']:.4f}/"
            f"{g_det['adv_shared_minus_ctrl_q90']:.4f}/"
            f"{g_det['sep_median_h1l1_minus_ctrl']:.6f}/"
            f"{g_det['control_median_gap']:.6f}"
        ),
        "",
        "## Real Stability",
        (
            f"- IID min/median/max: "
            f"{100.0 * real_iid_min:.2f}%/"
            f"{100.0 * np.median(iid_rates):.2f}%/"
            f"{100.0 * np.max(iid_rates):.2f}%"
        ),
        (
            f"- Block min/median/max: "
            f"{100.0 * real_block_min:.2f}%/"
            f"{100.0 * np.median(block_rates):.2f}%/"
            f"{100.0 * np.max(block_rates):.2f}%"
        ),
        "",
        "## Null Leakage",
        (
            f"- Random null mean/p90: "
            f"{100.0 * null_mean:.2f}%/"
            f"{100.0 * null_p90:.2f}%"
        ),
        (
            f"- Adversarial null mean/p90: "
            f"{100.0 * adv_mean:.2f}%/"
            f"{100.0 * adv_p90:.2f}%"
        ),
        "",
        "## Required Next Step",
        f"- {nxt}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1981] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1981] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1981] verdict={v}")


if __name__ == "__main__":
    main()
