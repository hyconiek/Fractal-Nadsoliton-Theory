#!/usr/bin/env python3
"""
QW-1978: Worst-case null contrast frontier.

Objective:
- continue from QW-1977,
- search operator points that maximize conservative real-vs-null separation,
- evaluate not only null mean but also null upper quantile (hard nulls).
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
IN_QW1977 = ROOT / "report_qw1977_contrast_first_global_search_gate.json"
OUT_JSON = ROOT / "report_qw1978_worstcase_null_contrast_frontier.json"
OUT_MD = ROOT / "RAPORT_QW1978_WORSTCASE_NULL_CONTRAST_FRONTIER.md"


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


def verdict(real_rate: float, null_mean: float, null_p90: float, conservative_contrast: float) -> Tuple[str, str]:
    if real_rate >= 0.90 and null_p90 <= 0.45 and conservative_contrast >= 0.30:
        return (
            "WORSTCASE_FRONTIER_IDENTIFIABLE",
            "PROMOTE_TO_PRE_EXTERNAL_FREEZE_AND_PACKAGE",
        )
    if real_rate >= 0.80 and null_p90 <= 0.60 and conservative_contrast >= 0.20:
        return (
            "WORSTCASE_FRONTIER_PARTIAL",
            "ADD_ONE_INFORMATIONAL_GUARD_AND_RUN_QW1979",
        )
    return (
        "WORSTCASE_FRONTIER_NON_IDENTIFIABLE",
        "NO_CLOSURE_CLAIM_YET",
    )


def main() -> None:
    r1970 = json.loads(IN_QW1970.read_text(encoding="utf-8"))
    r1977 = json.loads(IN_QW1977.read_text(encoding="utf-8"))
    thr = json.loads((ROOT / "report_qw1969_bootstrap_robust_recenter_search.json").read_text(encoding="utf-8"))["thresholds"]

    kernel = r1970["fixed_components"]["kernel"]
    params = r1970["fixed_components"]["params"]
    base_best = r1977["best"]

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

    # QW-1971 basis.
    c1_raw = pair_sign * lag_phase_sin * (corr0 - corr10)
    c2_raw = ctrl_mask * lag_phase_cos * (corr0 + corr10 + mean_abs)
    c1 = c1_raw / max(float(np.std(c1_raw)), 1e-12)
    c2 = c2_raw / max(float(np.std(c2_raw)), 1e-12)

    rng = np.random.default_rng(80780)
    n_random = 900
    xi1_samples = np.clip(rng.normal(loc=float(base_best["xi1"]), scale=0.00020, size=n_random), 0.0002, 0.0015)
    xi2_samples = np.clip(rng.normal(loc=float(base_best["xi2"]), scale=0.00025, size=n_random), -0.0013, -0.00002)

    stage_a = []
    for i in range(n_random):
        xi1 = float(xi1_samples[i])
        xi2 = float(xi2_samples[i])
        score = base_score + xi1 * c1 + xi2 * c2
        s_hl = score[pairs == 0]
        s_hv = score[pairs == 1]
        s_lv = score[pairs == 2]
        g = gw_metrics(s_hl, s_hv, s_lv)
        if all(gw_flags(g, thr).values()):
            stage_a.append(
                {
                    "xi1": xi1,
                    "xi2": xi2,
                    "gw_real": g,
                    "scores_real": (s_hl, s_hv, s_lv),
                }
            )

    # Stage B: conservative approximate frontier.
    stage_b = []
    rng_null = np.random.default_rng(81780)
    for i, row in enumerate(stage_a):
        xi1 = row["xi1"]
        xi2 = row["xi2"]
        s_hl, s_hv, s_lv = row["scores_real"]
        real_approx = bootstrap_pass_rate(s_hl, s_hv, s_lv, thr, n_boot=180, seed=82000 + i)

        null_rates = []
        for j in range(8):
            signs = np.array([1.0] * n_plus + [-1.0] * (n_ctrl - n_plus), dtype=float)
            rng_null.shuffle(signs)
            rand_sign = np.zeros_like(pair_sign)
            rand_sign[ctrl_idx] = signs
            c1n_raw = rand_sign * lag_phase_sin * (corr0 - corr10)
            c1n = c1n_raw / max(float(np.std(c1n_raw)), 1e-12)
            score_n = base_score + xi1 * c1n + xi2 * c2
            rb = bootstrap_pass_rate(
                score_n[pairs == 0], score_n[pairs == 1], score_n[pairs == 2], thr, n_boot=90, seed=83000 + i * 20 + j
            )
            null_rates.append(float(rb))

        null_mean = float(np.mean(null_rates))
        null_p90 = float(np.quantile(np.array(null_rates, dtype=float), 0.90))
        conservative_contrast = float(real_approx - null_p90)
        stage_b.append(
            {
                "xi1": xi1,
                "xi2": xi2,
                "gw_real": row["gw_real"],
                "real_approx": real_approx,
                "null_mean_approx": null_mean,
                "null_p90_approx": null_p90,
                "conservative_contrast_approx": conservative_contrast,
                "scores_real": row["scores_real"],
            }
        )

    stage_b.sort(
        key=lambda x: (
            x["conservative_contrast_approx"],
            x["real_approx"],
            -x["null_mean_approx"],
        ),
        reverse=True,
    )

    shortlist = [x for x in stage_b if x["real_approx"] >= 0.70][:10]

    # Stage C: full evaluation.
    final_rows = []
    rng_null_full = np.random.default_rng(84780)
    for i, row in enumerate(shortlist):
        xi1 = row["xi1"]
        xi2 = row["xi2"]
        s_hl, s_hv, s_lv = row["scores_real"]
        real_full = bootstrap_pass_rate(s_hl, s_hv, s_lv, thr, n_boot=1600, seed=85000 + i)

        null_full = []
        for j in range(24):
            signs = np.array([1.0] * n_plus + [-1.0] * (n_ctrl - n_plus), dtype=float)
            rng_null_full.shuffle(signs)
            rand_sign = np.zeros_like(pair_sign)
            rand_sign[ctrl_idx] = signs
            c1n_raw = rand_sign * lag_phase_sin * (corr0 - corr10)
            c1n = c1n_raw / max(float(np.std(c1n_raw)), 1e-12)
            score_n = base_score + xi1 * c1n + xi2 * c2
            rb = bootstrap_pass_rate(
                score_n[pairs == 0], score_n[pairs == 1], score_n[pairs == 2], thr, n_boot=220, seed=86000 + i * 40 + j
            )
            null_full.append(float(rb))

        null_mean = float(np.mean(null_full))
        null_p90 = float(np.quantile(np.array(null_full, dtype=float), 0.90))
        conservative_contrast = float(real_full - null_p90)
        final_rows.append(
            {
                "xi1": xi1,
                "xi2": xi2,
                "gw_real": row["gw_real"],
                "real_boot_pass_rate_1600": real_full,
                "null_boot_pass_rate_mean_24x220": null_mean,
                "null_boot_pass_rate_p90_24x220": null_p90,
                "null_boot_pass_rate_std_24x220": float(np.std(null_full)),
                "conservative_contrast_real_minus_null_p90": conservative_contrast,
                "contrast_real_minus_null_mean": float(real_full - null_mean),
            }
        )

    final_rows.sort(
        key=lambda x: (
            x["conservative_contrast_real_minus_null_p90"],
            x["real_boot_pass_rate_1600"],
            -x["null_boot_pass_rate_mean_24x220"],
        ),
        reverse=True,
    )

    if not final_rows:
        out = {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "verdict": "WORSTCASE_FRONTIER_NO_FINAL_CANDIDATE",
            "required_next_step": "WIDEN_SEARCH_AND_REDUCE_BASIS_FLEXIBILITY",
            "scan": {
                "n_random": int(n_random),
                "stage_a_pass": int(len(stage_a)),
                "stage_b_count": int(len(stage_b)),
                "shortlist_count": int(len(shortlist)),
            },
        }
        OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
        OUT_MD.write_text(
            "# RAPORT QW-1978: WORSTCASE NULL CONTRAST FRONTIER\n\n- Verdict: **WORSTCASE_FRONTIER_NO_FINAL_CANDIDATE**\n",
            encoding="utf-8",
        )
        print(f"[QW-1978] Saved JSON: {OUT_JSON.name}")
        print(f"[QW-1978] Saved MD:   {OUT_MD.name}")
        print("[QW-1978] verdict=WORSTCASE_FRONTIER_NO_FINAL_CANDIDATE")
        return

    best = final_rows[0]
    v, nxt = verdict(
        best["real_boot_pass_rate_1600"],
        best["null_boot_pass_rate_mean_24x220"],
        best["null_boot_pass_rate_p90_24x220"],
        best["conservative_contrast_real_minus_null_p90"],
    )

    baseline = r1977["best"]
    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": [IN_QW1970.name, IN_QW1977.name],
        "scan": {
            "n_random": int(n_random),
            "stage_a_pass": int(len(stage_a)),
            "stage_b_count": int(len(stage_b)),
            "shortlist_count": int(len(shortlist)),
            "final_count": int(len(final_rows)),
        },
        "baseline_qw1977_best": {
            "real_boot_pass_rate_2500": float(baseline["real_boot_pass_rate_2500"]),
            "null_boot_pass_rate_mean_20x500": float(baseline["null_boot_pass_rate_mean_20x500"]),
            "contrast_real_minus_null": float(baseline["contrast_real_minus_null"]),
        },
        "best": best,
        "top10": final_rows[:10],
        "verdict": v,
        "required_next_step": nxt,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1978: WORSTCASE NULL CONTRAST FRONTIER",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{v}**",
        "",
        "## Baseline (QW-1977) vs QW-1978 best",
        (
            f"- baseline real/null_mean/contrast: "
            f"{100.0 * out['baseline_qw1977_best']['real_boot_pass_rate_2500']:.2f}% / "
            f"{100.0 * out['baseline_qw1977_best']['null_boot_pass_rate_mean_20x500']:.2f}% / "
            f"{100.0 * out['baseline_qw1977_best']['contrast_real_minus_null']:.2f} pp"
        ),
        (
            f"- QW-1978 best real/null_mean/null_p90/cons_contrast: "
            f"{100.0 * best['real_boot_pass_rate_1600']:.2f}% / "
            f"{100.0 * best['null_boot_pass_rate_mean_24x220']:.2f}% / "
            f"{100.0 * best['null_boot_pass_rate_p90_24x220']:.2f}% / "
            f"{100.0 * best['conservative_contrast_real_minus_null_p90']:.2f} pp"
        ),
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
        "",
        "## Required Next Step",
        f"- {nxt}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1978] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1978] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1978] verdict={v}")


if __name__ == "__main__":
    main()
