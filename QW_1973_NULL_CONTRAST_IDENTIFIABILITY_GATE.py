#!/usr/bin/env python3
"""
QW-1973: Null-contrast identifiability gate for two-term GW structure.

Question:
- does the two-term model produce high pass-rate specifically for physical topology,
  or is it similarly strong under randomized null topology?
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
IN_QW1971 = ROOT / "report_qw1971_two_term_structural_control_dynamics_gate.json"
IN_QW1970 = ROOT / "report_qw1970_structural_gw_control_term_gate.json"
OUT_JSON = ROOT / "report_qw1973_null_contrast_identifiability_gate.json"
OUT_MD = ROOT / "RAPORT_QW1973_NULL_CONTRAST_IDENTIFIABILITY_GATE.md"


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


def deterministic_score(g: Dict[str, float], thr: Dict[str, float]) -> float:
    return float(
        1.0 * (g["auc_h1l1_vs_ctrl"] - thr["gw_auc_min"])
        + 1.0 * (g["adv_shared_minus_ctrl_q90"] - thr["gw_adv_min"])
        + 100.0 * (g["sep_median_h1l1_minus_ctrl"] - thr["gw_sep_min"])
        + 120.0 * (thr["gw_control_gap_max"] - g["control_median_gap"])
    )


def verdict(real_boot: float, null_boot_mean: float, contrast: float) -> Tuple[str, str]:
    if real_boot >= 0.95 and null_boot_mean <= 0.70 and contrast >= 0.25:
        return (
            "IDENTIFIABLE_PHYSICAL_SIGNAL_LOCKABLE",
            "FREEZE_AND_PREPARE_EXTERNAL_CONFIRMATORY_EXECUTION",
        )
    if real_boot >= 0.90 and contrast >= 0.10:
        return (
            "PARTIAL_IDENTIFIABILITY",
            "ADD_STRICTER_PHYSICAL_CONSTRAINTS_BEFORE_EXTERNAL",
        )
    return (
        "NON_IDENTIFIABLE_OVERFLEXIBLE_STRUCTURE",
        "DO_NOT_CLAIM_CLOSURE_UNTIL_NULL_CONTRAST_IMPROVES",
    )


def main() -> None:
    r1971 = json.loads(IN_QW1971.read_text(encoding="utf-8"))
    r1970 = json.loads(IN_QW1970.read_text(encoding="utf-8"))
    thr = json.loads((ROOT / "report_qw1969_bootstrap_robust_recenter_search.json").read_text(encoding="utf-8"))["thresholds"]

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
    ctrl_mask = np.where(pairs == 0, 0.0, 1.0)
    ctrl_idx = np.where(pairs != 0)[0]
    n_ctrl = len(ctrl_idx)
    n_plus = int(np.sum(pairs == 1))

    c1_raw = pair_sign * lag_phase_sin * (corr0 - corr10)
    c2_raw = ctrl_mask * lag_phase_cos * (corr0 + corr10 + mean_abs)
    c1 = c1_raw / max(float(np.std(c1_raw)), 1e-12)
    c2 = c2_raw / max(float(np.std(c2_raw)), 1e-12)

    # Candidate region from QW-1972 boundary.
    xi1_grid = np.linspace(0.000275, 0.001075, 33)
    xi2_grid = np.linspace(-0.0018, -0.0002, 33)

    rng_det = np.random.default_rng(29730)
    candidates = []
    for xi1 in xi1_grid:
        for xi2 in xi2_grid:
            score = base_score + float(xi1) * c1 + float(xi2) * c2
            s_hl = score[pairs == 0]
            s_hv = score[pairs == 1]
            s_lv = score[pairs == 2]
            g = gw_metrics(s_hl, s_hv, s_lv)
            f = gw_flags(g, thr)
            if not all(f.values()):
                continue

            # Null deterministic pass-rate for this (xi1, xi2).
            null_pass = 0
            for _ in range(40):
                signs = np.array([1.0] * n_plus + [-1.0] * (n_ctrl - n_plus), dtype=float)
                rng_det.shuffle(signs)
                rand_sign = np.zeros_like(pair_sign)
                rand_sign[ctrl_idx] = signs
                c1n_raw = rand_sign * lag_phase_sin * (corr0 - corr10)
                c1n = c1n_raw / max(float(np.std(c1n_raw)), 1e-12)
                score_n = base_score + float(xi1) * c1n + float(xi2) * c2
                g_n = gw_metrics(score_n[pairs == 0], score_n[pairs == 1], score_n[pairs == 2])
                if all(gw_flags(g_n, thr).values()):
                    null_pass += 1
            null_det_rate = float(null_pass / 40.0)

            candidates.append(
                {
                    "xi1": float(xi1),
                    "xi2": float(xi2),
                    "gw_real": g,
                    "real_det_score": deterministic_score(g, thr),
                    "null_det_pass_rate_40": null_det_rate,
                    "scores_real": (s_hl, s_hv, s_lv),
                }
            )

    # Prefer low null-pass first, then strong real deterministic margin.
    candidates.sort(key=lambda x: (x["null_det_pass_rate_40"], -x["real_det_score"]))
    shortlist = candidates[:25]

    evaluated = []
    rng_null_boot = np.random.default_rng(30730)
    for i, c in enumerate(shortlist):
        xi1 = c["xi1"]
        xi2 = c["xi2"]
        s_hl, s_hv, s_lv = c["scores_real"]
        real_boot = bootstrap_pass_rate(s_hl, s_hv, s_lv, thr, n_boot=2000, seed=31000 + i)

        null_boot_rates = []
        for j in range(20):
            signs = np.array([1.0] * n_plus + [-1.0] * (n_ctrl - n_plus), dtype=float)
            rng_null_boot.shuffle(signs)
            rand_sign = np.zeros_like(pair_sign)
            rand_sign[ctrl_idx] = signs
            c1n_raw = rand_sign * lag_phase_sin * (corr0 - corr10)
            c1n = c1n_raw / max(float(np.std(c1n_raw)), 1e-12)
            score_n = base_score + float(xi1) * c1n + float(xi2) * c2
            rb = bootstrap_pass_rate(
                score_n[pairs == 0], score_n[pairs == 1], score_n[pairs == 2], thr, n_boot=500, seed=32000 + i * 20 + j
            )
            null_boot_rates.append(float(rb))

        null_boot_mean = float(np.mean(null_boot_rates))
        contrast = float(real_boot - null_boot_mean)
        evaluated.append(
            {
                "xi1": xi1,
                "xi2": xi2,
                "gw_real": c["gw_real"],
                "null_det_pass_rate_40": c["null_det_pass_rate_40"],
                "real_boot_pass_rate_2000": real_boot,
                "null_boot_pass_rate_mean_20x500": null_boot_mean,
                "null_boot_pass_rate_std_20x500": float(np.std(null_boot_rates)),
                "contrast_real_minus_null": contrast,
            }
        )

    evaluated.sort(
        key=lambda x: (
            x["contrast_real_minus_null"],
            x["real_boot_pass_rate_2000"],
            -x["null_boot_pass_rate_mean_20x500"],
        ),
        reverse=True,
    )
    best = evaluated[0]
    v, nxt = verdict(
        best["real_boot_pass_rate_2000"],
        best["null_boot_pass_rate_mean_20x500"],
        best["contrast_real_minus_null"],
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": [IN_QW1970.name, IN_QW1971.name],
        "grid": {
            "xi1_min": float(np.min(xi1_grid)),
            "xi1_max": float(np.max(xi1_grid)),
            "xi2_min": float(np.min(xi2_grid)),
            "xi2_max": float(np.max(xi2_grid)),
            "grid_size": int(len(xi1_grid) * len(xi2_grid)),
            "deterministic_pass_candidates": int(len(candidates)),
            "shortlist_size": int(len(shortlist)),
        },
        "best": best,
        "top10": evaluated[:10],
        "verdict": v,
        "required_next_step": nxt,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1973: NULL CONTRAST IDENTIFIABILITY GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{v}**",
        "",
        "## Best Candidate",
        f"- xi1/xi2: {best['xi1']:.6f}/{best['xi2']:.6f}",
        (
            f"- real GW auc/adv/sep/gap: "
            f"{best['gw_real']['auc_h1l1_vs_ctrl']:.4f}/"
            f"{best['gw_real']['adv_shared_minus_ctrl_q90']:.4f}/"
            f"{best['gw_real']['sep_median_h1l1_minus_ctrl']:.6f}/"
            f"{best['gw_real']['control_median_gap']:.6f}"
        ),
        f"- real bootstrap pass (2000): {100.0 * best['real_boot_pass_rate_2000']:.2f}%",
        f"- null bootstrap mean (20x500): {100.0 * best['null_boot_pass_rate_mean_20x500']:.2f}%",
        f"- contrast (real-null): {100.0 * best['contrast_real_minus_null']:.2f} pp",
        "",
        "## Required Next Step",
        f"- {nxt}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1973] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1973] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1973] verdict={v}")


if __name__ == "__main__":
    main()
