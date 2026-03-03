#!/usr/bin/env python3
"""
QW-1980: Signed dual-basis independent gate.

Hypothesis:
- previous lock/coupling constraints may over-restrict feasible region,
- two independent signed control terms can improve real-vs-null tradeoff.
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
IN_QW1979 = ROOT / "report_qw1979_real_null_pareto_feasibility_map.json"
OUT_JSON = ROOT / "report_qw1980_signed_dual_basis_independent_gate.json"
OUT_MD = ROOT / "RAPORT_QW1980_SIGNED_DUAL_BASIS_INDEPENDENT_GATE.md"


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


def verdict(strict_count: int, medium_count: int) -> Tuple[str, str]:
    if strict_count > 0:
        return (
            "SIGNED_DUAL_STRICT_FEASIBLE",
            "PROMOTE_STRICT_CANDIDATE_TO_FULL_EXTERNAL_PREP",
        )
    if medium_count > 0:
        return (
            "SIGNED_DUAL_MEDIUM_FEASIBLE",
            "ADD_FINAL_INFO_GUARD_AND_REPEAT_WITH_HIGHER_POWER",
        )
    return (
        "SIGNED_DUAL_NON_FEASIBLE",
        "OPERATOR_CLASS_STILL_INSUFFICIENT_FOR_STRICT_CLOSURE",
    )


def main() -> None:
    r1970 = json.loads(IN_QW1970.read_text(encoding="utf-8"))
    r1979 = json.loads(IN_QW1979.read_text(encoding="utf-8"))
    thr = json.loads((ROOT / "report_qw1969_bootstrap_robust_recenter_search.json").read_text(encoding="utf-8"))["thresholds"]

    kernel = r1970["fixed_components"]["kernel"]
    params = r1970["fixed_components"]["params"]
    center = r1979["best_conservative_contrast"]

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

    # Signed dual basis (independent coefficients).
    c1_raw = pair_sign * lag_phase_sin * (corr0 - corr10)
    c3_raw = pair_sign * lag_phase_cos * (mean_abs + np.abs(corr0) + np.abs(corr10))
    c1 = c1_raw / max(float(np.std(c1_raw)), 1e-12)
    c3 = c3_raw / max(float(np.std(c3_raw)), 1e-12)

    rng = np.random.default_rng(100800)
    n_random = 900
    xi1 = np.clip(rng.normal(float(center["xi1"]), 0.00025, size=n_random), 0.0001, 0.0022)
    xi3 = np.clip(rng.normal(0.00035, 0.00030, size=n_random), -0.0001, 0.0020)

    stage_a = []
    for i in range(n_random):
        a = float(xi1[i])
        b = float(xi3[i])
        score = base_score + a * c1 + b * c3
        s_hl = score[pairs == 0]
        s_hv = score[pairs == 1]
        s_lv = score[pairs == 2]
        g = gw_metrics(s_hl, s_hv, s_lv)
        if all(gw_flags(g, thr).values()):
            stage_a.append(
                {
                    "xi1": a,
                    "xi3": b,
                    "gw_real": g,
                    "scores_real": (s_hl, s_hv, s_lv),
                }
            )

    # Approximate ranking.
    stage_b = []
    rng_null = np.random.default_rng(101800)
    for i, row in enumerate(stage_a):
        a = row["xi1"]
        b = row["xi3"]
        s_hl, s_hv, s_lv = row["scores_real"]
        real_approx = bootstrap_pass_rate(s_hl, s_hv, s_lv, thr, n_boot=120, seed=102000 + i)

        null_rates = []
        for j in range(6):
            signs = np.array([1.0] * n_plus + [-1.0] * (n_ctrl - n_plus), dtype=float)
            rng_null.shuffle(signs)
            rand_sign = np.zeros_like(pair_sign)
            rand_sign[ctrl_idx] = signs
            c1n_raw = rand_sign * lag_phase_sin * (corr0 - corr10)
            c3n_raw = rand_sign * lag_phase_cos * (mean_abs + np.abs(corr0) + np.abs(corr10))
            c1n = c1n_raw / max(float(np.std(c1n_raw)), 1e-12)
            c3n = c3n_raw / max(float(np.std(c3n_raw)), 1e-12)
            score_n = base_score + a * c1n + b * c3n
            rb = bootstrap_pass_rate(
                score_n[pairs == 0], score_n[pairs == 1], score_n[pairs == 2], thr, n_boot=80, seed=103000 + i * 20 + j
            )
            null_rates.append(float(rb))

        null_mean = float(np.mean(null_rates))
        null_p90 = float(np.quantile(np.array(null_rates, dtype=float), 0.90))
        stage_b.append(
            {
                "xi1": a,
                "xi3": b,
                "gw_real": row["gw_real"],
                "real_approx": real_approx,
                "null_mean_approx": null_mean,
                "null_p90_approx": null_p90,
                "cons_approx": float(real_approx - null_p90),
                "scores_real": row["scores_real"],
            }
        )

    stage_b.sort(key=lambda x: (x["cons_approx"], x["real_approx"], -x["null_mean_approx"]), reverse=True)
    shortlist = [x for x in stage_b if x["real_approx"] >= 0.70][:10]

    # Full evaluation.
    final_rows = []
    rng_null_full = np.random.default_rng(104800)
    for i, row in enumerate(shortlist):
        a = row["xi1"]
        b = row["xi3"]
        s_hl, s_hv, s_lv = row["scores_real"]
        real_full = bootstrap_pass_rate(s_hl, s_hv, s_lv, thr, n_boot=1600, seed=105000 + i)

        null_full = []
        for j in range(20):
            signs = np.array([1.0] * n_plus + [-1.0] * (n_ctrl - n_plus), dtype=float)
            rng_null_full.shuffle(signs)
            rand_sign = np.zeros_like(pair_sign)
            rand_sign[ctrl_idx] = signs
            c1n_raw = rand_sign * lag_phase_sin * (corr0 - corr10)
            c3n_raw = rand_sign * lag_phase_cos * (mean_abs + np.abs(corr0) + np.abs(corr10))
            c1n = c1n_raw / max(float(np.std(c1n_raw)), 1e-12)
            c3n = c3n_raw / max(float(np.std(c3n_raw)), 1e-12)
            score_n = base_score + a * c1n + b * c3n
            rb = bootstrap_pass_rate(
                score_n[pairs == 0], score_n[pairs == 1], score_n[pairs == 2], thr, n_boot=200, seed=106000 + i * 30 + j
            )
            null_full.append(float(rb))

        null_mean = float(np.mean(null_full))
        null_p90 = float(np.quantile(np.array(null_full, dtype=float), 0.90))
        final_rows.append(
            {
                "xi1": a,
                "xi3": b,
                "gw_real": row["gw_real"],
                "real_boot_pass_rate_1600": real_full,
                "null_boot_pass_rate_mean_20x200": null_mean,
                "null_boot_pass_rate_p90_20x200": null_p90,
                "contrast_conservative": float(real_full - null_p90),
                "contrast_mean": float(real_full - null_mean),
            }
        )

    final_rows.sort(key=lambda x: (x["contrast_conservative"], x["real_boot_pass_rate_1600"], -x["null_boot_pass_rate_mean_20x200"]), reverse=True)

    strict = [r for r in final_rows if r["real_boot_pass_rate_1600"] >= 0.90 and r["null_boot_pass_rate_p90_20x200"] <= 0.45]
    medium = [r for r in final_rows if r["real_boot_pass_rate_1600"] >= 0.85 and r["null_boot_pass_rate_p90_20x200"] <= 0.50]

    v, nxt = verdict(len(strict), len(medium))
    best = final_rows[0] if final_rows else None

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": [IN_QW1970.name, IN_QW1979.name],
        "scan": {
            "n_random": int(n_random),
            "stage_a_pass": int(len(stage_a)),
            "stage_b_count": int(len(stage_b)),
            "shortlist_count": int(len(shortlist)),
            "final_count": int(len(final_rows)),
        },
        "feasibility_counts": {
            "strict_real_ge_0p90_nullp90_le_0p45": int(len(strict)),
            "medium_real_ge_0p85_nullp90_le_0p50": int(len(medium)),
        },
        "best": best,
        "top10": final_rows[:10],
        "verdict": v,
        "required_next_step": nxt,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1980: SIGNED DUAL BASIS INDEPENDENT GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{v}**",
        "",
        "## Feasibility Counts",
        f"- strict (real>=0.90, null_p90<=0.45): {len(strict)}",
        f"- medium (real>=0.85, null_p90<=0.50): {len(medium)}",
        "",
        "## Best Candidate",
        (
            f"- real/null_mean/null_p90/cons: "
            f"{100.0 * best['real_boot_pass_rate_1600']:.2f}%/"
            f"{100.0 * best['null_boot_pass_rate_mean_20x200']:.2f}%/"
            f"{100.0 * best['null_boot_pass_rate_p90_20x200']:.2f}%/"
            f"{100.0 * best['contrast_conservative']:.2f} pp"
        )
        if best
        else "- none",
        "",
        "## Required Next Step",
        f"- {nxt}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1980] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1980] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1980] verdict={v}")


if __name__ == "__main__":
    main()
