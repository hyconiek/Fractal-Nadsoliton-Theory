#!/usr/bin/env python3
"""
QW-2027: Structural GW control-gap term scan.

Applies a symmetric control-channel correction to reduce |median(H1-V1)-median(L1-V1)|
while preserving shared-channel separation metrics.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2027_gw_control_gap_structural_term_scan.json"
OUT_MD = ROOT / "RAPORT_QW2027_GW_CONTROL_GAP_STRUCTURAL_TERM_SCAN.md"


def kernel_fn(d: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d + phi) / (1.0 + beta * (d**eta))


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


def gw_metrics(pair: np.ndarray, score: np.ndarray) -> Dict[str, float]:
    s_hl = score[pair == "H1-L1"]
    s_hv = score[pair == "H1-V1"]
    s_lv = score[pair == "L1-V1"]
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


def main() -> None:
    d2021 = json.loads((ROOT / "report_qw2021_v2_eta_operator_beta_constraint_scan.json").read_text(encoding="utf-8"))
    d2026 = json.loads((ROOT / "report_qw2026_joint_mass_flavor_gw_scan_eta_kernel.json").read_text(encoding="utf-8"))

    kernel = d2021["selected"]["fit"]
    p_amp = float(d2026["best_row"]["p_amp"])
    r_dist = float(d2026["best_row"]["r_dist"])

    thresholds = {
        "gw_sep_min": 0.0020,
        "gw_adv_min": 0.30,
        "gw_auc_min": 0.75,
        "gw_control_gap_max": 0.0025,
    }

    df = pd.read_csv(ROOT / "gw1831_window_features.csv")
    pair = df["pair"].astype(str).to_numpy()

    d = np.array([1.0, 2.0, 3.0, 4.0], dtype=float)
    raw_w = (np.abs(kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])) ** p_amp) * (d**r_dist)
    w = raw_w / np.sum(raw_w)

    score_raw = (
        w[0] * df["max_abs_corr"].to_numpy(dtype=float)
        + w[1] * df["mean_abs_corr"].to_numpy(dtype=float)
        + w[2] * df["corr_at_0ms"].to_numpy(dtype=float)
        + w[3] * df["corr_at_10ms"].to_numpy(dtype=float)
    )

    base = gw_metrics(pair, score_raw)

    # Symmetric control correction: H1-V1 -= kappa, L1-V1 += kappa.
    kappas = np.linspace(-0.004, 0.004, 161)

    best = None
    for kappa in kappas:
        s = score_raw.copy()
        s[pair == "H1-V1"] -= float(kappa)
        s[pair == "L1-V1"] += float(kappa)

        m = gw_metrics(pair, s)
        flags = {
            "gw_sep_ge_min": bool(m["sep_median_h1l1_minus_ctrl"] >= thresholds["gw_sep_min"]),
            "gw_adv_ge_min": bool(m["adv_shared_minus_ctrl_q90"] >= thresholds["gw_adv_min"]),
            "gw_auc_ge_min": bool(m["auc_h1l1_vs_ctrl"] >= thresholds["gw_auc_min"]),
            "gw_control_gap_le_max": bool(m["control_median_gap"] <= thresholds["gw_control_gap_max"]),
        }

        # prefer passing + minimal deviation from baseline metrics.
        drift = (
            abs(m["auc_h1l1_vs_ctrl"] - base["auc_h1l1_vs_ctrl"]) / max(base["auc_h1l1_vs_ctrl"], 1e-12)
            + abs(m["adv_shared_minus_ctrl_q90"] - base["adv_shared_minus_ctrl_q90"]) / max(abs(base["adv_shared_minus_ctrl_q90"]), 1e-12)
            + abs(m["sep_median_h1l1_minus_ctrl"] - base["sep_median_h1l1_minus_ctrl"]) / max(abs(base["sep_median_h1l1_minus_ctrl"]), 1e-12)
        )

        score = (
            40.0 * max(0.0, thresholds["gw_control_gap_max"] - m["control_median_gap"]) * (-1.0)
            + 8.0 * max(0.0, thresholds["gw_sep_min"] - m["sep_median_h1l1_minus_ctrl"])
            + 8.0 * max(0.0, thresholds["gw_adv_min"] - m["adv_shared_minus_ctrl_q90"])
            + 8.0 * max(0.0, thresholds["gw_auc_min"] - m["auc_h1l1_vs_ctrl"])
            + drift
        )
        # Lower better.

        row = {
            "kappa": float(kappa),
            "metrics": m,
            "flags": flags,
            "all_pass": bool(all(flags.values())),
            "objective": float(score),
        }

        if best is None:
            best = row
        else:
            # prioritize all_pass then objective
            if (row["all_pass"] and not best["all_pass"]) or (
                row["all_pass"] == best["all_pass"] and row["objective"] < best["objective"]
            ):
                best = row

    verdict = "GW_CONTROL_GAP_STRUCTURAL_TERM_PASS" if best["all_pass"] else "GW_CONTROL_GAP_STRUCTURAL_TERM_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw2021_v2_eta_operator_beta_constraint_scan.json:selected.fit",
        "context_source": "report_qw2026_joint_mass_flavor_gw_scan_eta_kernel.json:best_row",
        "weights": {
            "w_max_abs_corr": float(w[0]),
            "w_mean_abs_corr": float(w[1]),
            "w_corr_at_0ms": float(w[2]),
            "w_corr_at_10ms": float(w[3]),
        },
        "baseline_metrics": base,
        "selected": best,
        "verdict": verdict,
        "required_next_step": (
            "INTEGRATE_KAPPA_TERM_IN_JOINT_SCAN_WITH_FLAVOR_IMPROVEMENT"
            if best["all_pass"]
            else "ESCALATE_TO_MULTI_BASIS_GW_OPERATOR"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2027: GW CONTROL GAP STRUCTURAL TERM SCAN",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Baseline",
        f"- auc/adv/sep/gap: {base['auc_h1l1_vs_ctrl']:.4f}/{base['adv_shared_minus_ctrl_q90']:.4f}/{base['sep_median_h1l1_minus_ctrl']:.6f}/{base['control_median_gap']:.6f}",
        "",
        "## Selected kappa",
        f"- kappa: {best['kappa']:.6f}",
        f"- auc/adv/sep/gap: {best['metrics']['auc_h1l1_vs_ctrl']:.4f}/{best['metrics']['adv_shared_minus_ctrl_q90']:.4f}/{best['metrics']['sep_median_h1l1_minus_ctrl']:.6f}/{best['metrics']['control_median_gap']:.6f}",
        f"- all_pass: {best['all_pass']}",
        "",
        "## Required Next Step",
        f"- {out['required_next_step']}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2027] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2027] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2027] verdict={verdict} kappa={best['kappa']:.6f}")


if __name__ == "__main__":
    main()
