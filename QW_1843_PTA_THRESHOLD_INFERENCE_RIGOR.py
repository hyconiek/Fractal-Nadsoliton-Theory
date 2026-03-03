#!/usr/bin/env python3
"""
QW-1843: PTA threshold inference rigor audit.

Checks inferential strength for frozen QW-1839 PTA thresholds using QW-1823
replication-level quantile gains.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1843_pta_threshold_inference_rigor.json"
OUT_MD = ROOT / "RAPORT_QW1843_PTA_THRESHOLD_INFERENCE_RIGOR.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def binom_tail(n: int, k: int, p: float) -> float:
    s = 0.0
    for i in range(k, n + 1):
        s += math.comb(n, i) * (p**i) * ((1.0 - p) ** (n - i))
    return float(s)


def lower_bound_one_sided(n: int, k: int, alpha: float = 0.05) -> float:
    """Exact one-sided lower confidence bound for Bernoulli p.

    Returns p_low such that P[X>=k | p_low] = alpha.
    """
    if k <= 0:
        return 0.0
    if k >= n:
        return float(alpha ** (1.0 / n))

    lo, hi = 0.0, 1.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        t = binom_tail(n, k, mid)
        if t < alpha:
            lo = mid
        else:
            hi = mid
    return float(0.5 * (lo + hi))


def bootstrap_ci(arr: np.ndarray, stat: str, n_mc: int, seed: int) -> Tuple[float, float, float]:
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(arr), size=(n_mc, len(arr)))
    samp = arr[idx]

    if stat == "mean":
        vals = np.mean(samp, axis=1)
    elif stat == "std":
        vals = np.std(samp, axis=1, ddof=1)
    else:
        raise ValueError(stat)

    return float(np.mean(vals)), float(np.quantile(vals, 0.025)), float(np.quantile(vals, 0.975))


def main() -> None:
    d1823 = load("report_qw1823_quantile_score_oos.json")
    d1839 = load("report_qw1839_joint_confirmatory_prereg_protocol.json")

    t = d1839["protocol"]["pta_protocol"]["thresholds"]
    thr_mean = float(t["mean_quantile_gain_m2_minus_m2e_min"])
    thr_prob = float(t["prob_quantile_gain_positive_min"])
    thr_std = float(t["std_quantile_gain_m2_minus_m2e_max"])

    rep = d1823.get("replications", [])
    gain = np.array([float(x["quantile_gain_m2_minus_m2e"]) for x in rep], dtype=float)

    n = int(len(gain))
    k_pos = int(np.sum(gain > 0.0))

    mean_obs = float(np.mean(gain))
    std_obs = float(np.std(gain, ddof=1))
    prob_obs = float(k_pos / n)

    mean_boot_mu, mean_ci_lo, mean_ci_hi = bootstrap_ci(gain, "mean", n_mc=100000, seed=21927)
    std_boot_mu, std_ci_lo, std_ci_hi = bootstrap_ci(gain, "std", n_mc=100000, seed=21928)

    pval_prob_vs_thr = binom_tail(n=n, k=k_pos, p=thr_prob)
    lower95_prob = lower_bound_one_sided(n=n, k=k_pos, alpha=0.05)

    # Minimal n (assuming all-positive outcomes) to reject H0: p <= thr_prob at alpha=0.05
    n_min_all_positive = None
    for n_try in range(1, 500):
        if (thr_prob**n_try) <= 0.05:
            n_min_all_positive = n_try
            break

    pass_mean = mean_obs >= thr_mean and mean_ci_lo >= thr_mean
    pass_std = std_obs <= thr_std and std_ci_hi <= thr_std
    pass_prob_point = prob_obs >= thr_prob
    pass_prob_inferential = pval_prob_vs_thr <= 0.05 and lower95_prob >= thr_prob

    if pass_mean and pass_std and pass_prob_inferential:
        verdict = "PTA_THRESHOLD_INFERENCE_STRONG"
    elif pass_mean and pass_std and pass_prob_point and not pass_prob_inferential:
        verdict = "PTA_POINT_PASS_BUT_PROBABILITY_UNDERPOWERED"
    else:
        verdict = "PTA_THRESHOLD_INFERENCE_INSUFFICIENT"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_replications": n,
        "observed": {
            "mean_quantile_gain": mean_obs,
            "std_quantile_gain": std_obs,
            "prob_quantile_gain_positive": prob_obs,
            "k_positive": k_pos,
        },
        "thresholds": {
            "mean_min": thr_mean,
            "std_max": thr_std,
            "prob_min": thr_prob,
        },
        "bootstrap": {
            "mean": {
                "bootstrap_mean": mean_boot_mu,
                "ci95": [mean_ci_lo, mean_ci_hi],
            },
            "std": {
                "bootstrap_mean": std_boot_mu,
                "ci95": [std_ci_lo, std_ci_hi],
            },
        },
        "probability_inference": {
            "pvalue_h0_p_le_threshold": pval_prob_vs_thr,
            "lower95_one_sided_for_p": lower95_prob,
            "n_min_all_positive_for_alpha0p05_vs_threshold": n_min_all_positive,
        },
        "pass_flags": {
            "mean_threshold_with_ci": bool(pass_mean),
            "std_threshold_with_ci": bool(pass_std),
            "prob_threshold_point": bool(pass_prob_point),
            "prob_threshold_inferential": bool(pass_prob_inferential),
        },
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1843: PTA THRESHOLD INFERENCE RIGOR",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- n replications: {n}",
        f"- Observed mean/std/prob+: {mean_obs:.6f} / {std_obs:.6f} / {prob_obs:.3f}",
        f"- Thresholds mean>=/std<=/prob>=: {thr_mean:.3f} / {thr_std:.3f} / {thr_prob:.3f}",
        f"- Mean CI95: [{mean_ci_lo:.6f}, {mean_ci_hi:.6f}]",
        f"- Std CI95: [{std_ci_lo:.6f}, {std_ci_hi:.6f}]",
        (
            "- Prob-inference: "
            f"p-value(H0:p<=thr)={pval_prob_vs_thr:.6f}, "
            f"one-sided lower95={lower95_prob:.6f}"
        ),
        (
            "- n_min(all-positive) for alpha=0.05 vs prob-threshold: "
            f"{n_min_all_positive}"
        ),
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- mean_threshold_with_ci: {pass_mean}",
        f"- std_threshold_with_ci: {pass_std}",
        f"- prob_threshold_point: {pass_prob_point}",
        f"- prob_threshold_inferential: {pass_prob_inferential}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1843] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1843] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
