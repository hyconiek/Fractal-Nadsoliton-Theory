#!/usr/bin/env python3
"""
QW-1975: Adaptive physics-constrained identifiability gate.

Changes vs QW-1974:
- keep kernel coupling lock xi2 = -rho * xi1,
- replace hard physical guards with soft penalties,
- preserve strict GW threshold gate,
- optimize for real-null contrast under constrained class.
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
IN_QW1973 = ROOT / "report_qw1973_null_contrast_identifiability_gate.json"
OUT_JSON = ROOT / "report_qw1975_adaptive_physics_constrained_gate.json"
OUT_MD = ROOT / "RAPORT_QW1975_ADAPTIVE_PHYSICS_CONSTRAINED_GATE.md"


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


def safe_corr(x: np.ndarray, y: np.ndarray) -> float:
    x0 = x - np.mean(x)
    y0 = y - np.mean(y)
    sx = float(np.std(x0))
    sy = float(np.std(y0))
    if sx < 1e-12 or sy < 1e-12:
        return 0.0
    return float(np.mean(x0 * y0) / (sx * sy))


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


def base_det_score(g: Dict[str, float], thr: Dict[str, float]) -> float:
    return float(
        1.1 * (g["auc_h1l1_vs_ctrl"] - thr["gw_auc_min"])
        + 1.0 * (g["adv_shared_minus_ctrl_q90"] - thr["gw_adv_min"])
        + 120.0 * (g["sep_median_h1l1_minus_ctrl"] - thr["gw_sep_min"])
        + 180.0 * (thr["gw_control_gap_max"] - g["control_median_gap"])
    )


def guard_penalty(control_balance_abs: float, hl_corr: float, ctrl_corr: float) -> float:
    # Soft physics penalties (0 when meeting target region).
    return float(
        280.0 * max(0.0, control_balance_abs - 0.0022)
        + 2.8 * max(0.0, 0.03 - hl_corr)
        + 2.0 * max(0.0, 0.01 - ctrl_corr)
    )


def verdict(best_real: float, best_null: float, best_contrast: float) -> Tuple[str, str]:
    if best_real >= 0.90 and best_null <= 0.35 and best_contrast >= 0.45:
        return (
            "ADAPTIVE_PHYSICS_IDENTIFIABLE",
            "PROMOTE_TO_PRE_EXTERNAL_FREEZE_PROTOCOL",
        )
    if best_real >= 0.85 and best_null <= 0.45 and best_contrast >= 0.35:
        return (
            "ADAPTIVE_PHYSICS_PARTIAL_IDENTIFIABILITY",
            "ADD_FINAL_GEOMETRIC_LOCK_AND_RUN_QW1976",
        )
    return (
        "ADAPTIVE_PHYSICS_STILL_NON_IDENTIFIABLE",
        "REWORK_OPERATOR_BASIS_BEFORE_CLOSURE_CLAIM",
    )


def main() -> None:
    r1971 = json.loads(IN_QW1971.read_text(encoding="utf-8"))
    r1970 = json.loads(IN_QW1970.read_text(encoding="utf-8"))
    r1973 = json.loads(IN_QW1973.read_text(encoding="utf-8"))
    thr = json.loads((ROOT / "report_qw1969_bootstrap_robust_recenter_search.json").read_text(encoding="utf-8"))["thresholds"]

    kernel = r1970["fixed_components"]["kernel"]
    params = r1970["fixed_components"]["params"]
    xi1_center = float(r1971["best"]["xi1"])

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

    rho = float(
        ((1.0 + abs(np.sin(kernel["phi"]))) / (1.0 + abs(np.cos(kernel["phi"]))))
        * ((1.0 + kernel["beta"]) / (1.0 + 0.5 * kernel["eta"]))
    )

    xi1_grid = np.linspace(max(0.0, xi1_center - 0.0006), xi1_center + 0.0010, 401)
    stage_rows = []
    gw_pass_count = 0

    for xi1 in xi1_grid:
        xi2 = -rho * float(xi1)
        score = base_score + float(xi1) * c1 + float(xi2) * c2
        s_hl = score[pairs == 0]
        s_hv = score[pairs == 1]
        s_lv = score[pairs == 2]
        g = gw_metrics(s_hl, s_hv, s_lv)
        f = gw_flags(g, thr)
        if not all(f.values()):
            continue
        gw_pass_count += 1

        cb = float(abs(np.median(s_hv) + np.median(s_lv)))
        hlc = safe_corr(s_hl, lag_phase_sin[pairs == 0])
        csc = safe_corr(
            np.concatenate([s_hv, s_lv]),
            np.concatenate([lag_phase_sin[pairs == 1], -lag_phase_sin[pairs == 2]]),
        )
        pen = guard_penalty(cb, hlc, csc)
        score_det = base_det_score(g, thr) - pen

        stage_rows.append(
            {
                "xi1": float(xi1),
                "xi2": float(xi2),
                "gw_real": g,
                "guard_metrics": {
                    "control_balance_abs": cb,
                    "hl_lag_corr": hlc,
                    "ctrl_signed_lag_corr": csc,
                    "guard_penalty": pen,
                },
                "det_score": score_det,
                "scores_real": (s_hl, s_hv, s_lv),
            }
        )

    stage_rows.sort(key=lambda x: x["det_score"], reverse=True)
    shortlist = stage_rows[:24]

    rng_null = np.random.default_rng(50750)
    evaluated = []
    for i, row in enumerate(shortlist):
        xi1 = row["xi1"]
        xi2 = row["xi2"]
        s_hl, s_hv, s_lv = row["scores_real"]
        real_boot = bootstrap_pass_rate(s_hl, s_hv, s_lv, thr, n_boot=1600, seed=51000 + i)

        null_boots = []
        for j in range(15):
            signs = np.array([1.0] * n_plus + [-1.0] * (n_ctrl - n_plus), dtype=float)
            rng_null.shuffle(signs)
            rand_sign = np.zeros_like(pair_sign)
            rand_sign[ctrl_idx] = signs
            c1n_raw = rand_sign * lag_phase_sin * (corr0 - corr10)
            c1n = c1n_raw / max(float(np.std(c1n_raw)), 1e-12)
            score_n = base_score + float(xi1) * c1n + float(xi2) * c2
            rb = bootstrap_pass_rate(
                score_n[pairs == 0], score_n[pairs == 1], score_n[pairs == 2], thr, n_boot=400, seed=52000 + i * 30 + j
            )
            null_boots.append(float(rb))

        null_mean = float(np.mean(null_boots))
        contrast = float(real_boot - null_mean)
        evaluated.append(
            {
                "xi1": xi1,
                "xi2": xi2,
                "gw_real": row["gw_real"],
                "guard_metrics": row["guard_metrics"],
                "real_boot_pass_rate_1600": real_boot,
                "null_boot_pass_rate_mean_15x400": null_mean,
                "null_boot_pass_rate_std_15x400": float(np.std(null_boots)),
                "contrast_real_minus_null": contrast,
            }
        )

    evaluated.sort(
        key=lambda x: (
            x["contrast_real_minus_null"],
            x["real_boot_pass_rate_1600"],
            -x["null_boot_pass_rate_mean_15x400"],
        ),
        reverse=True,
    )

    if not evaluated:
        out = {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "verdict": "ADAPTIVE_PHYSICS_NO_EVALUATED_CANDIDATE",
            "required_next_step": "REWORK_CONSTRAINT_LOCK_OR_OPERATOR",
            "kernel_coupling_lock": {"rho": rho},
            "scan": {"xi1_size": int(len(xi1_grid)), "gw_pass_count": int(gw_pass_count)},
        }
        OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
        OUT_MD.write_text(
            "# RAPORT QW-1975: ADAPTIVE PHYSICS CONSTRAINED GATE\n\n- Verdict: **ADAPTIVE_PHYSICS_NO_EVALUATED_CANDIDATE**\n",
            encoding="utf-8",
        )
        print(f"[QW-1975] Saved JSON: {OUT_JSON.name}")
        print(f"[QW-1975] Saved MD:   {OUT_MD.name}")
        print("[QW-1975] verdict=ADAPTIVE_PHYSICS_NO_EVALUATED_CANDIDATE")
        return

    best = evaluated[0]
    v, nxt = verdict(
        best["real_boot_pass_rate_1600"],
        best["null_boot_pass_rate_mean_15x400"],
        best["contrast_real_minus_null"],
    )

    baseline = r1973["best"]
    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": [IN_QW1970.name, IN_QW1971.name, IN_QW1973.name],
        "kernel_coupling_lock": {"formula": "xi2 = -rho * xi1", "rho": rho},
        "scan": {
            "xi1_min": float(np.min(xi1_grid)),
            "xi1_max": float(np.max(xi1_grid)),
            "xi1_size": int(len(xi1_grid)),
            "gw_pass_count": int(gw_pass_count),
            "shortlist_size": int(len(shortlist)),
            "evaluated_size": int(len(evaluated)),
        },
        "baseline_qw1973_best": {
            "real_boot_pass_rate_2000": float(baseline["real_boot_pass_rate_2000"]),
            "null_boot_pass_rate_mean_20x500": float(baseline["null_boot_pass_rate_mean_20x500"]),
            "contrast_real_minus_null": float(baseline["contrast_real_minus_null"]),
        },
        "best": best,
        "top10": evaluated[:10],
        "verdict": v,
        "required_next_step": nxt,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1975: ADAPTIVE PHYSICS CONSTRAINED GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{v}**",
        "",
        "## Coupling Lock",
        f"- xi2 = -rho * xi1, rho={rho:.6f}",
        "",
        "## Baseline (QW-1973) vs QW-1975 best",
        (
            f"- baseline real/null/contrast: "
            f"{100.0 * out['baseline_qw1973_best']['real_boot_pass_rate_2000']:.2f}% / "
            f"{100.0 * out['baseline_qw1973_best']['null_boot_pass_rate_mean_20x500']:.2f}% / "
            f"{100.0 * out['baseline_qw1973_best']['contrast_real_minus_null']:.2f} pp"
        ),
        (
            f"- QW-1975 best real/null/contrast: "
            f"{100.0 * best['real_boot_pass_rate_1600']:.2f}% / "
            f"{100.0 * best['null_boot_pass_rate_mean_15x400']:.2f}% / "
            f"{100.0 * best['contrast_real_minus_null']:.2f} pp"
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
        (
            f"- guard metrics: "
            f"penalty={best['guard_metrics']['guard_penalty']:.4f}, "
            f"control_balance_abs={best['guard_metrics']['control_balance_abs']:.6f}, "
            f"hl_lag_corr={best['guard_metrics']['hl_lag_corr']:.4f}, "
            f"ctrl_signed_lag_corr={best['guard_metrics']['ctrl_signed_lag_corr']:.4f}"
        ),
        "",
        "## Required Next Step",
        f"- {nxt}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1975] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1975] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1975] verdict={v}")


if __name__ == "__main__":
    main()
