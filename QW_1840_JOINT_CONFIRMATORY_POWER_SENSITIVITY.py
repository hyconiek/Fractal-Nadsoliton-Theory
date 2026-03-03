#!/usr/bin/env python3
"""
QW-1840: Joint confirmatory power and margin-consumption sensitivity.

Purpose:
- quantify how likely the frozen QW-1839 joint gate is to pass under
  current empirical uncertainty (bootstrap over QW-1823 and QW-1836 folds),
- estimate robustness under controlled effect drift expressed as
  "consumed fraction of current margin to threshold".
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1840_joint_confirmatory_power_sensitivity.json"
OUT_MD = ROOT / "RAPORT_QW1840_JOINT_CONFIRMATORY_POWER_SENSITIVITY.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def ci95(p: float, n: int) -> Tuple[float, float]:
    if n <= 0:
        return (0.0, 0.0)
    se = math.sqrt(max(0.0, p * (1.0 - p) / n))
    lo = max(0.0, p - 1.96 * se)
    hi = min(1.0, p + 1.96 * se)
    return lo, hi


def adjust_pta_quantile_gain(
    base_gain: np.ndarray,
    frac_margin: float,
    mean_threshold: float,
    std_threshold: float,
) -> np.ndarray:
    """
    Drift model for PTA quantile gain:
    - consume fraction of mean margin by shifting mean downward,
    - consume fraction of std margin by inflating dispersion upward.
    """
    mu0 = float(np.mean(base_gain))
    sd0 = float(np.std(base_gain))

    mean_margin = mu0 - mean_threshold
    std_margin = std_threshold - sd0

    delta_mu = frac_margin * mean_margin
    mu_target = mu0 - delta_mu

    if sd0 <= 1e-15:
        scale = 1.0
    else:
        sd_target = max(1e-12, sd0 + frac_margin * std_margin)
        scale = sd_target / sd0

    centered = base_gain - mu0
    adjusted = centered * scale + mu_target
    return adjusted


def adjust_gw_arrays(
    auc: np.ndarray,
    adv: np.ndarray,
    gap: np.ndarray,
    frac_margin: float,
    auc_threshold: float,
    adv_threshold: float,
    gap_threshold: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Drift model for GW metrics:
    - beneficial metrics (AUC, ADV) move down by consumed margin,
    - cost metric (gap) moves up by consumed margin.
    """
    auc_margin = float(np.mean(auc)) - auc_threshold
    adv_margin = float(np.mean(adv)) - adv_threshold
    gap_margin = gap_threshold - float(np.mean(gap))

    auc_adj = auc - frac_margin * auc_margin
    adv_adj = adv - frac_margin * adv_margin
    gap_adj = gap + frac_margin * gap_margin
    return auc_adj, adv_adj, gap_adj


def bootstrap_joint(
    gain: np.ndarray,
    auc: np.ndarray,
    adv: np.ndarray,
    gap: np.ndarray,
    n_rep: int,
    n_fold: int,
    n_mc: int,
    seed: int,
    th_mean_gain: float,
    th_prob_pos: float,
    th_std_max: float,
    th_auc_min: float,
    th_adv_min: float,
    th_gap_max: float,
    th_prob_adv_pos: float,
) -> Dict:
    rng = np.random.default_rng(seed)

    idx_rep = rng.integers(0, gain.size, size=(n_mc, n_rep))
    s_gain = gain[idx_rep]

    pta_mean = s_gain.mean(axis=1)
    pta_prob_pos = (s_gain > 0.0).mean(axis=1)
    pta_std = s_gain.std(axis=1)

    pta_pass = (
        (pta_mean >= th_mean_gain)
        & (pta_prob_pos >= th_prob_pos)
        & (pta_std <= th_std_max)
    )

    idx_fold = rng.integers(0, auc.size, size=(n_mc, n_fold))
    s_auc = auc[idx_fold]
    s_adv = adv[idx_fold]
    s_gap = gap[idx_fold]

    gw_auc = s_auc.mean(axis=1)
    gw_adv = s_adv.mean(axis=1)
    gw_gap = s_gap.mean(axis=1)
    gw_prob_adv_pos = (s_adv > 0.0).mean(axis=1)

    gw_pass = (
        (gw_auc >= th_auc_min)
        & (gw_adv >= th_adv_min)
        & (gw_gap <= th_gap_max)
        & (gw_prob_adv_pos >= th_prob_adv_pos)
    )

    joint_pass = pta_pass & gw_pass

    p_pta = float(np.mean(pta_pass))
    p_gw = float(np.mean(gw_pass))
    p_joint = float(np.mean(joint_pass))

    pta_lo, pta_hi = ci95(p_pta, n_mc)
    gw_lo, gw_hi = ci95(p_gw, n_mc)
    joint_lo, joint_hi = ci95(p_joint, n_mc)

    return {
        "n_mc": n_mc,
        "n_rep": n_rep,
        "n_fold": n_fold,
        "pass_probability": {
            "pta": p_pta,
            "gw": p_gw,
            "joint": p_joint,
        },
        "pass_probability_ci95": {
            "pta": [pta_lo, pta_hi],
            "gw": [gw_lo, gw_hi],
            "joint": [joint_lo, joint_hi],
        },
        "metric_means": {
            "pta_mean_quantile_gain": float(np.mean(pta_mean)),
            "pta_prob_positive": float(np.mean(pta_prob_pos)),
            "pta_std_quantile_gain": float(np.mean(pta_std)),
            "gw_mean_auc": float(np.mean(gw_auc)),
            "gw_mean_adv": float(np.mean(gw_adv)),
            "gw_mean_gap": float(np.mean(gw_gap)),
            "gw_prob_adv_positive": float(np.mean(gw_prob_adv_pos)),
        },
    }


def main() -> None:
    d1839 = load_json("report_qw1839_joint_confirmatory_prereg_protocol.json")
    d1823 = load_json("report_qw1823_quantile_score_oos.json")
    d1836 = load_json("report_qw1836_gw_control_calibrated_objective.json")

    if d1839.get("verdict") != "JOINT_CONFIRMATORY_PREREG_FROZEN":
        raise RuntimeError("QW-1839 prereg protocol is not frozen.")

    pta_t = d1839["protocol"]["pta_protocol"]["thresholds"]
    gw_t = d1839["protocol"]["gw_protocol"]["thresholds"]

    rep = d1823.get("replications", [])
    folds = d1836.get("folds", [])
    if not rep or not folds:
        raise RuntimeError("Missing replications/folds in input reports.")

    gain = np.array([float(x["quantile_gain_m2_minus_m2e"]) for x in rep], dtype=float)
    auc = np.array([float(x["cal_auc"]) for x in folds], dtype=float)
    adv = np.array([float(x["cal_adv"]) for x in folds], dtype=float)
    gap = np.array([float(x["cal_control_gap"]) for x in folds], dtype=float)

    n_rep = len(gain)
    n_fold = len(auc)

    margins = {
        "pta_mean_margin": float(np.mean(gain) - pta_t["mean_quantile_gain_m2_minus_m2e_min"]),
        "pta_std_margin": float(pta_t["std_quantile_gain_m2_minus_m2e_max"] - np.std(gain)),
        "gw_auc_margin": float(np.mean(auc) - gw_t["calibrated_mean_auc_min"]),
        "gw_adv_margin": float(np.mean(adv) - gw_t["calibrated_mean_adv_min"]),
        "gw_gap_margin": float(gw_t["calibrated_mean_control_gap_max"] - np.mean(gap)),
    }

    scenario_fracs = [0.0, 0.5, 0.8, 1.0, 1.2]
    scenario_results = []

    for i, frac in enumerate(scenario_fracs):
        gain_adj = adjust_pta_quantile_gain(
            base_gain=gain,
            frac_margin=frac,
            mean_threshold=pta_t["mean_quantile_gain_m2_minus_m2e_min"],
            std_threshold=pta_t["std_quantile_gain_m2_minus_m2e_max"],
        )
        auc_adj, adv_adj, gap_adj = adjust_gw_arrays(
            auc=auc,
            adv=adv,
            gap=gap,
            frac_margin=frac,
            auc_threshold=gw_t["calibrated_mean_auc_min"],
            adv_threshold=gw_t["calibrated_mean_adv_min"],
            gap_threshold=gw_t["calibrated_mean_control_gap_max"],
        )

        boot = bootstrap_joint(
            gain=gain_adj,
            auc=auc_adj,
            adv=adv_adj,
            gap=gap_adj,
            n_rep=n_rep,
            n_fold=n_fold,
            n_mc=30000,
            seed=21925 + i,
            th_mean_gain=pta_t["mean_quantile_gain_m2_minus_m2e_min"],
            th_prob_pos=pta_t["prob_quantile_gain_positive_min"],
            th_std_max=pta_t["std_quantile_gain_m2_minus_m2e_max"],
            th_auc_min=gw_t["calibrated_mean_auc_min"],
            th_adv_min=gw_t["calibrated_mean_adv_min"],
            th_gap_max=gw_t["calibrated_mean_control_gap_max"],
            th_prob_adv_pos=gw_t["calibrated_prob_adv_positive_min"],
        )

        scenario_results.append(
            {
                "scenario": f"margin_consumption_{int(frac * 100)}pct",
                "frac_margin_consumed": frac,
                **boot,
            }
        )

    # robustness scan for the joint gate at high resolution
    scan_fracs = np.linspace(0.0, 1.4, 29)
    scan = []
    for i, frac in enumerate(scan_fracs):
        gain_adj = adjust_pta_quantile_gain(
            base_gain=gain,
            frac_margin=float(frac),
            mean_threshold=pta_t["mean_quantile_gain_m2_minus_m2e_min"],
            std_threshold=pta_t["std_quantile_gain_m2_minus_m2e_max"],
        )
        auc_adj, adv_adj, gap_adj = adjust_gw_arrays(
            auc=auc,
            adv=adv,
            gap=gap,
            frac_margin=float(frac),
            auc_threshold=gw_t["calibrated_mean_auc_min"],
            adv_threshold=gw_t["calibrated_mean_adv_min"],
            gap_threshold=gw_t["calibrated_mean_control_gap_max"],
        )
        boot = bootstrap_joint(
            gain=gain_adj,
            auc=auc_adj,
            adv=adv_adj,
            gap=gap_adj,
            n_rep=n_rep,
            n_fold=n_fold,
            n_mc=10000,
            seed=32925 + i,
            th_mean_gain=pta_t["mean_quantile_gain_m2_minus_m2e_min"],
            th_prob_pos=pta_t["prob_quantile_gain_positive_min"],
            th_std_max=pta_t["std_quantile_gain_m2_minus_m2e_max"],
            th_auc_min=gw_t["calibrated_mean_auc_min"],
            th_adv_min=gw_t["calibrated_mean_adv_min"],
            th_gap_max=gw_t["calibrated_mean_control_gap_max"],
            th_prob_adv_pos=gw_t["calibrated_prob_adv_positive_min"],
        )
        scan.append(
            {
                "frac_margin_consumed": float(frac),
                "joint_pass_probability": boot["pass_probability"]["joint"],
            }
        )

    frac_joint_ge_90 = max(
        (x["frac_margin_consumed"] for x in scan if x["joint_pass_probability"] >= 0.90),
        default=0.0,
    )
    frac_joint_ge_50 = max(
        (x["frac_margin_consumed"] for x in scan if x["joint_pass_probability"] >= 0.50),
        default=0.0,
    )

    nominal = scenario_results[0]
    verdict = "JOINT_CONFIRMATORY_HIGH_POWER" if nominal["pass_probability"]["joint"] >= 0.9 else "JOINT_CONFIRMATORY_POWER_LIMITED"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "qw1839_protocol_sha256": d1839.get("protocol_sha256"),
            "n_replications_pta": n_rep,
            "n_folds_gw": n_fold,
        },
        "threshold_margins": margins,
        "scenario_results": scenario_results,
        "robustness_scan": {
            "scan_points": scan,
            "max_margin_consumption_with_joint_power_ge_0p90": float(frac_joint_ge_90),
            "max_margin_consumption_with_joint_power_ge_0p50": float(frac_joint_ge_50),
        },
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1840: JOINT CONFIRMATORY POWER SENSITIVITY",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Protocol SHA256 (QW-1839): `{out['inputs']['qw1839_protocol_sha256']}`",
        f"- Rozmiary bootstrap: PTA n={n_rep}, GW n={n_fold}",
        "",
        "## Marginesy do progow (stan referencyjny)",
        f"- PTA mean margin: {margins['pta_mean_margin']:.6f}",
        f"- PTA std margin: {margins['pta_std_margin']:.6f}",
        f"- GW AUC margin: {margins['gw_auc_margin']:.6f}",
        f"- GW ADV margin: {margins['gw_adv_margin']:.6f}",
        f"- GW GAP margin: {margins['gw_gap_margin']:.6f}",
        "",
        "## Wyniki scenariuszy (consumed margin)",
    ]

    for row in scenario_results:
        p = row["pass_probability"]
        ci = row["pass_probability_ci95"]
        lines += [
            (
                f"- {row['scenario']}: p(PTA)={p['pta']:.3f} "
                f"[{ci['pta'][0]:.3f},{ci['pta'][1]:.3f}], "
                f"p(GW)={p['gw']:.3f} [{ci['gw'][0]:.3f},{ci['gw'][1]:.3f}], "
                f"p(JOINT)={p['joint']:.3f} [{ci['joint'][0]:.3f},{ci['joint'][1]:.3f}]"
            )
        ]

    lines += [
        "",
        "## Robustnosc joint gate",
        f"- Max consumed margin z p(JOINT)>=0.90: {frac_joint_ge_90:.2f}",
        f"- Max consumed margin z p(JOINT)>=0.50: {frac_joint_ge_50:.2f}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1840] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1840] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
