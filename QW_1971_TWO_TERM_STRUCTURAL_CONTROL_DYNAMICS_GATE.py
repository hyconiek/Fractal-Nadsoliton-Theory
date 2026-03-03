#!/usr/bin/env python3
"""
QW-1971: Two-term shared structural control dynamics gate for GW channel.

Objective:
- improve bootstrap triad pass-rate beyond QW-1970,
- keep frozen kernel + fixed mass/flavor branch,
- no sector retune.
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
OUT_JSON = ROOT / "report_qw1971_two_term_structural_control_dynamics_gate.json"
OUT_MD = ROOT / "RAPORT_QW1971_TWO_TERM_STRUCTURAL_CONTROL_DYNAMICS_GATE.md"


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


def gw_metrics_from_scores(s_hl: np.ndarray, s_hv: np.ndarray, s_lv: np.ndarray) -> Dict[str, float]:
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


def gw_flags(gw: Dict[str, float], thr: Dict[str, float]) -> Dict[str, bool]:
    return {
        "gw_sep_ge_min": bool(gw["sep_median_h1l1_minus_ctrl"] >= thr["gw_sep_min"]),
        "gw_adv_ge_min": bool(gw["adv_shared_minus_ctrl_q90"] >= thr["gw_adv_min"]),
        "gw_auc_ge_min": bool(gw["auc_h1l1_vs_ctrl"] >= thr["gw_auc_min"]),
        "gw_control_gap_le_max": bool(gw["control_median_gap"] <= thr["gw_control_gap_max"]),
    }


def deterministic_score(gw: Dict[str, float], thr: Dict[str, float]) -> float:
    # Higher is better.
    return float(
        1.2 * (gw["auc_h1l1_vs_ctrl"] - thr["gw_auc_min"])
        + 1.0 * (gw["adv_shared_minus_ctrl_q90"] - thr["gw_adv_min"])
        + 120.0 * (gw["sep_median_h1l1_minus_ctrl"] - thr["gw_sep_min"])
        + 180.0 * (thr["gw_control_gap_max"] - gw["control_median_gap"])
    )


def bootstrap_gw_pass_rate(
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
        g = gw_metrics_from_scores(b_hl, b_hv, b_lv)
        if all(gw_flags(g, thr).values()):
            pass_count += 1
    return float(pass_count / n_boot)


def verdict(best_rate: float, baseline_rate: float) -> Tuple[str, str]:
    if best_rate >= 0.95:
        return (
            "TWO_TERM_STRUCTURAL_LOCKABLE",
            "FREEZE_TWO_TERM_STRUCTURE_AND_RUN_TRUE_EXTERNAL_CONFIRMATORY",
        )
    if best_rate >= 0.80 and best_rate >= baseline_rate + 0.01:
        return (
            "TWO_TERM_STRUCTURAL_PARTIAL_SUCCESS",
            "PROMOTE_TWO_TERM_DYNAMICS_AND_RUN_STRICT_CROSS_CHECK",
        )
    return (
        "TWO_TERM_STRUCTURAL_INSUFFICIENT",
        "NEED_HIGHER_ORDER_SHARED_DYNAMICS_OR_NEW_OBSERVABLES",
    )


def main() -> None:
    r1970 = json.loads(IN_QW1970.read_text(encoding="utf-8"))
    kernel = r1970["fixed_components"]["kernel"]
    params = r1970["fixed_components"]["params"]
    thr = r1970["source_thresholds"] if "source_thresholds" in r1970 else json.loads(
        (ROOT / "report_qw1969_bootstrap_robust_recenter_search.json").read_text(encoding="utf-8")
    )["thresholds"]

    df = pd.read_csv(ROOT / "gw1831_window_features.csv")
    pair_map = {"H1-L1": 0, "H1-V1": 1, "L1-V1": 2}
    pairs = df["pair"].map(pair_map).to_numpy(dtype=int)
    feats = df[["max_abs_corr", "mean_abs_corr", "corr_at_0ms", "corr_at_10ms"]].to_numpy(dtype=float)

    # Base score from QW-1969 recentered params.
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

    # Term 1 (antisymmetric control mode), consistent with QW-1970.
    c1_raw = pair_sign * lag_phase_sin * (corr0 - corr10)
    c1 = c1_raw / max(float(np.std(c1_raw)), 1e-12)

    # Term 2 (common control mode), shared and kernel-phase-linked.
    c2_raw = ctrl_mask * lag_phase_cos * (corr0 + corr10 + mean_abs)
    c2 = c2_raw / max(float(np.std(c2_raw)), 1e-12)

    # Baseline from QW-1970 best.
    baseline_rate = float(r1970["best"]["boot_pass_rate_5000"])
    baseline_xi1 = float(r1970["best"]["xi"])

    xi1_grid = np.linspace(baseline_xi1 - 0.0010, baseline_xi1 + 0.0010, 81)
    xi2_grid = np.linspace(-0.0010, 0.0010, 81)

    rows = []
    for xi1 in xi1_grid:
        for xi2 in xi2_grid:
            score = base_score + float(xi1) * c1 + float(xi2) * c2
            s_hl = score[pairs == 0]
            s_hv = score[pairs == 1]
            s_lv = score[pairs == 2]
            gw = gw_metrics_from_scores(s_hl, s_hv, s_lv)
            flags = gw_flags(gw, thr)
            rows.append(
                {
                    "xi1": float(xi1),
                    "xi2": float(xi2),
                    "gw": gw,
                    "flags": flags,
                    "all_pass": bool(all(flags.values())),
                    "det_score": deterministic_score(gw, thr),
                    "scores": {"s_hl": s_hl, "s_hv": s_hv, "s_lv": s_lv},
                }
            )

    pass_rows = [r for r in rows if r["all_pass"]]
    if not pass_rows:
        out = {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "verdict": "TWO_TERM_STRUCTURAL_NO_DETERMINISTIC_PASS",
            "required_next_step": "CHANGE_TWO_TERM_BASIS",
        }
        OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
        OUT_MD.write_text(
            "# RAPORT QW-1971: TWO-TERM STRUCTURAL CONTROL DYNAMICS GATE\n\n- Verdict: **TWO_TERM_STRUCTURAL_NO_DETERMINISTIC_PASS**\n",
            encoding="utf-8",
        )
        print(f"[QW-1971] Saved JSON: {OUT_JSON.name}")
        print(f"[QW-1971] Saved MD:   {OUT_MD.name}")
        print("[QW-1971] verdict=TWO_TERM_STRUCTURAL_NO_DETERMINISTIC_PASS")
        return

    pass_rows.sort(key=lambda x: x["det_score"], reverse=True)
    top = pass_rows[:50]

    screened = []
    for i, r in enumerate(top):
        boot = bootstrap_gw_pass_rate(
            s_hl=r["scores"]["s_hl"],
            s_hv=r["scores"]["s_hv"],
            s_lv=r["scores"]["s_lv"],
            thr=thr,
            n_boot=2000,
            seed=21000 + i,
        )
        screened.append(
            {
                "xi1": r["xi1"],
                "xi2": r["xi2"],
                "gw": r["gw"],
                "det_score": r["det_score"],
                "boot_pass_rate_2000": float(boot),
            }
        )
    screened.sort(key=lambda x: x["boot_pass_rate_2000"], reverse=True)

    finalists = []
    for j, r in enumerate(screened[:5]):
        xi1 = r["xi1"]
        xi2 = r["xi2"]
        score = base_score + float(xi1) * c1 + float(xi2) * c2
        s_hl = score[pairs == 0]
        s_hv = score[pairs == 1]
        s_lv = score[pairs == 2]
        gw = gw_metrics_from_scores(s_hl, s_hv, s_lv)
        boot = bootstrap_gw_pass_rate(
            s_hl=s_hl,
            s_hv=s_hv,
            s_lv=s_lv,
            thr=thr,
            n_boot=5000,
            seed=22000 + j,
        )
        finalists.append(
            {
                "xi1": float(xi1),
                "xi2": float(xi2),
                "gw": gw,
                "boot_pass_rate_5000": float(boot),
            }
        )
    finalists.sort(key=lambda x: x["boot_pass_rate_5000"], reverse=True)
    best = finalists[0]

    # local neighborhood in (xi1, xi2).
    local_points = [
        (best["xi1"] - 0.0003, best["xi2"]),
        (best["xi1"], best["xi2"] - 0.0003),
        (best["xi1"], best["xi2"]),
        (best["xi1"], best["xi2"] + 0.0003),
        (best["xi1"] + 0.0003, best["xi2"]),
    ]
    local_rows = []
    for k, (xi1, xi2) in enumerate(local_points):
        score = base_score + float(xi1) * c1 + float(xi2) * c2
        s_hl = score[pairs == 0]
        s_hv = score[pairs == 1]
        s_lv = score[pairs == 2]
        gw = gw_metrics_from_scores(s_hl, s_hv, s_lv)
        flags = gw_flags(gw, thr)
        boot = bootstrap_gw_pass_rate(
            s_hl=s_hl,
            s_hv=s_hv,
            s_lv=s_lv,
            thr=thr,
            n_boot=3000,
            seed=23000 + k,
        )
        local_rows.append(
            {
                "xi1": float(xi1),
                "xi2": float(xi2),
                "gw": gw,
                "all_pass_deterministic": bool(all(flags.values())),
                "boot_pass_rate_3000": float(boot),
            }
        )

    v, nxt = verdict(best["boot_pass_rate_5000"], baseline_rate)

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_report": IN_QW1970.name,
        "baseline": {
            "xi1_from_qw1970": baseline_xi1,
            "bootstrap_pass_rate_5000": baseline_rate,
        },
        "grid": {
            "xi1_min": float(np.min(xi1_grid)),
            "xi1_max": float(np.max(xi1_grid)),
            "xi1_size": int(len(xi1_grid)),
            "xi2_min": float(np.min(xi2_grid)),
            "xi2_max": float(np.max(xi2_grid)),
            "xi2_size": int(len(xi2_grid)),
            "deterministic_pass_count": int(len(pass_rows)),
        },
        "screened_top10": screened[:10],
        "finalists": finalists,
        "best": best,
        "local_robustness": local_rows,
        "verdict": v,
        "required_next_step": nxt,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1971: TWO-TERM STRUCTURAL CONTROL DYNAMICS GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{v}**",
        "",
        "## Baseline vs Best Two-Term",
        f"- baseline (QW-1970) bootstrap pass 5000: {100.0 * baseline_rate:.2f}%",
        f"- best two-term bootstrap pass 5000: {100.0 * best['boot_pass_rate_5000']:.2f}%",
        f"- delta: {100.0 * (best['boot_pass_rate_5000'] - baseline_rate):.2f} pp",
        "",
        "## Best Two-Term Candidate",
        f"- xi1: {best['xi1']:.6f}",
        f"- xi2: {best['xi2']:.6f}",
        (
            f"- GW auc/adv/sep/gap: "
            f"{best['gw']['auc_h1l1_vs_ctrl']:.4f}/"
            f"{best['gw']['adv_shared_minus_ctrl_q90']:.4f}/"
            f"{best['gw']['sep_median_h1l1_minus_ctrl']:.6f}/"
            f"{best['gw']['control_median_gap']:.6f}"
        ),
        "",
        "## Local Robustness (xi1, xi2 neighborhood)",
        *[
            (
                f"- (xi1,xi2)=({r['xi1']:.6f},{r['xi2']:.6f}): "
                f"det_pass={r['all_pass_deterministic']}, "
                f"boot3000={100.0 * r['boot_pass_rate_3000']:.2f}%"
            )
            for r in local_rows
        ],
        "",
        "## Required Next Step",
        f"- {nxt}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1971] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1971] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1971] verdict={v}")


if __name__ == "__main__":
    main()
