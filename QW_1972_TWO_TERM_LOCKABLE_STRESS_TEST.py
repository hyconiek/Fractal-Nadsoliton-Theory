#!/usr/bin/env python3
"""
QW-1972: Stress test for QW-1971 two-term structural lock candidate.

Checks:
- multi-seed IID bootstrap stability,
- block bootstrap stability (temporal dependence proxy),
- boundary scan around xi1/xi2 optimum,
- null topology test (randomized control-sign structure).
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
OUT_JSON = ROOT / "report_qw1972_two_term_lockable_stress_test.json"
OUT_MD = ROOT / "RAPORT_QW1972_TWO_TERM_LOCKABLE_STRESS_TEST.md"


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


def bootstrap_iid_pass_rate(
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
    idx: List[int] = []
    while len(idx) < n:
        start = int(rng.integers(0, n))
        block = [(start + i) % n for i in range(block_len)]
        idx.extend(block)
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


def deterministic_score(g: Dict[str, float], thr: Dict[str, float]) -> float:
    return float(
        1.0 * (g["auc_h1l1_vs_ctrl"] - thr["gw_auc_min"])
        + 1.0 * (g["adv_shared_minus_ctrl_q90"] - thr["gw_adv_min"])
        + 100.0 * (g["sep_median_h1l1_minus_ctrl"] - thr["gw_sep_min"])
        + 140.0 * (thr["gw_control_gap_max"] - g["control_median_gap"])
    )


def verdict(iid_min: float, block_min: float, null_det_pass_rate: float) -> Tuple[str, str]:
    if iid_min >= 0.98 and block_min >= 0.95 and null_det_pass_rate <= 0.10:
        return (
            "TWO_TERM_LOCK_CANDIDATE_PASSES_STRESS",
            "PROCEED_TO_TRUE_EXTERNAL_CONFIRMATORY_PACKAGE",
        )
    if iid_min >= 0.90 and block_min >= 0.85:
        return (
            "TWO_TERM_LOCK_CANDIDATE_PARTIAL_STRESS_PASS",
            "RUN_MORE_STRINGENT_NULLS_AND_EXTERNAL_PREP",
        )
    return (
        "TWO_TERM_LOCK_CANDIDATE_STRESS_FAIL",
        "REWORK_SHARED_CONTROL_DYNAMICS_BEFORE_EXTERNAL",
    )


def main() -> None:
    r1971 = json.loads(IN_QW1971.read_text(encoding="utf-8"))
    r1970 = json.loads(IN_QW1970.read_text(encoding="utf-8"))
    thr = json.loads((ROOT / "report_qw1969_bootstrap_robust_recenter_search.json").read_text(encoding="utf-8"))["thresholds"]

    best_xi1 = float(r1971["best"]["xi1"])
    best_xi2 = float(r1971["best"]["xi2"])
    baseline_rate = float(r1971["baseline"]["bootstrap_pass_rate_5000"])

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

    c1_raw = pair_sign * lag_phase_sin * (corr0 - corr10)
    c2_raw = ctrl_mask * lag_phase_cos * (corr0 + corr10 + mean_abs)
    c1 = c1_raw / max(float(np.std(c1_raw)), 1e-12)
    c2 = c2_raw / max(float(np.std(c2_raw)), 1e-12)

    score_best = base_score + best_xi1 * c1 + best_xi2 * c2
    s_hl_best = score_best[pairs == 0]
    s_hv_best = score_best[pairs == 1]
    s_lv_best = score_best[pairs == 2]
    gw_best = gw_metrics(s_hl_best, s_hv_best, s_lv_best)
    det_flags_best = gw_flags(gw_best, thr)

    # 1) Multi-seed IID bootstrap.
    iid_seeds = [24000 + i for i in range(10)]
    iid_rates = [
        bootstrap_iid_pass_rate(s_hl_best, s_hv_best, s_lv_best, thr, n_boot=3000, seed=s)
        for s in iid_seeds
    ]

    # 2) Block bootstrap.
    block_seeds = [25000 + i for i in range(8)]
    block_rates = [
        bootstrap_block_pass_rate(
            s_hl_best,
            s_hv_best,
            s_lv_best,
            thr,
            n_boot=2000,
            block_len=8,
            seed=s,
        )
        for s in block_seeds
    ]

    # 3) Boundary scan around best (ensure not narrow boundary artifact).
    xi1_grid = np.linspace(best_xi1 - 0.0004, best_xi1 + 0.0004, 33)
    xi2_grid = np.linspace(-0.0018, -0.0002, 33)
    pass_rows = []
    for xi1 in xi1_grid:
        for xi2 in xi2_grid:
            score = base_score + float(xi1) * c1 + float(xi2) * c2
            s_hl = score[pairs == 0]
            s_hv = score[pairs == 1]
            s_lv = score[pairs == 2]
            g = gw_metrics(s_hl, s_hv, s_lv)
            f = gw_flags(g, thr)
            if all(f.values()):
                pass_rows.append(
                    {
                        "xi1": float(xi1),
                        "xi2": float(xi2),
                        "gw": g,
                        "det_score": deterministic_score(g, thr),
                        "scores": (s_hl, s_hv, s_lv),
                    }
                )
    pass_rows.sort(key=lambda x: x["det_score"], reverse=True)
    top_boundary = pass_rows[:5]
    boundary_boot = []
    for i, r in enumerate(top_boundary):
        boot = bootstrap_iid_pass_rate(
            r["scores"][0], r["scores"][1], r["scores"][2], thr, n_boot=3000, seed=26000 + i
        )
        boundary_boot.append(
            {
                "xi1": r["xi1"],
                "xi2": r["xi2"],
                "gw": r["gw"],
                "boot_pass_rate_3000": float(boot),
            }
        )

    # 4) Null topology test: randomize control-sign structure in c1.
    ctrl_idx = np.where(pairs != 0)[0]
    n_ctrl = len(ctrl_idx)
    n_plus = int(np.sum(pairs == 1))
    null_det_pass = 0
    null_rows = []
    rng = np.random.default_rng(27000)
    for i in range(200):
        signs = np.array([1.0] * n_plus + [-1.0] * (n_ctrl - n_plus), dtype=float)
        rng.shuffle(signs)
        rand_sign = np.zeros_like(pair_sign)
        rand_sign[ctrl_idx] = signs

        c1n_raw = rand_sign * lag_phase_sin * (corr0 - corr10)
        c1n = c1n_raw / max(float(np.std(c1n_raw)), 1e-12)
        score_n = base_score + best_xi1 * c1n + best_xi2 * c2
        s_hl_n = score_n[pairs == 0]
        s_hv_n = score_n[pairs == 1]
        s_lv_n = score_n[pairs == 2]
        g_n = gw_metrics(s_hl_n, s_hv_n, s_lv_n)
        f_n = gw_flags(g_n, thr)
        det_pass = bool(all(f_n.values()))
        if det_pass:
            null_det_pass += 1

        if i < 20:
            boot_n = bootstrap_iid_pass_rate(s_hl_n, s_hv_n, s_lv_n, thr, n_boot=1000, seed=28000 + i)
            null_rows.append(
                {
                    "iter": int(i),
                    "det_pass": det_pass,
                    "gw": g_n,
                    "boot_pass_rate_1000": float(boot_n),
                }
            )

    null_det_pass_rate = float(null_det_pass / 200.0)

    iid_min = float(np.min(iid_rates))
    block_min = float(np.min(block_rates))
    v, nxt = verdict(iid_min, block_min, null_det_pass_rate)

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": [IN_QW1970.name, IN_QW1971.name],
        "best_candidate": {
            "xi1": best_xi1,
            "xi2": best_xi2,
            "gw": gw_best,
            "deterministic_flags": det_flags_best,
            "baseline_bootstrap_pass_5000_from_qw1971": baseline_rate,
        },
        "iid_bootstrap_stress": {
            "seeds": iid_seeds,
            "n_boot_each": 3000,
            "rates": iid_rates,
            "min_rate": iid_min,
            "median_rate": float(np.median(iid_rates)),
            "max_rate": float(np.max(iid_rates)),
        },
        "block_bootstrap_stress": {
            "seeds": block_seeds,
            "n_boot_each": 2000,
            "block_len": 8,
            "rates": block_rates,
            "min_rate": block_min,
            "median_rate": float(np.median(block_rates)),
            "max_rate": float(np.max(block_rates)),
        },
        "boundary_scan": {
            "xi1_min": float(np.min(xi1_grid)),
            "xi1_max": float(np.max(xi1_grid)),
            "xi2_min": float(np.min(xi2_grid)),
            "xi2_max": float(np.max(xi2_grid)),
            "deterministic_pass_count": int(len(pass_rows)),
            "top_bootstrap": boundary_boot,
        },
        "null_topology_test": {
            "n_null": 200,
            "deterministic_pass_count": int(null_det_pass),
            "deterministic_pass_rate": null_det_pass_rate,
            "sample_bootstrap_rows": null_rows,
        },
        "verdict": v,
        "required_next_step": nxt,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1972: TWO-TERM LOCKABLE STRESS TEST",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{v}**",
        "",
        "## Best Candidate (from QW-1971)",
        f"- xi1/xi2: {best_xi1:.6f}/{best_xi2:.6f}",
        (
            f"- GW auc/adv/sep/gap: "
            f"{gw_best['auc_h1l1_vs_ctrl']:.4f}/"
            f"{gw_best['adv_shared_minus_ctrl_q90']:.4f}/"
            f"{gw_best['sep_median_h1l1_minus_ctrl']:.6f}/"
            f"{gw_best['control_median_gap']:.6f}"
        ),
        "",
        "## IID Bootstrap Stress",
        (
            f"- min/median/max pass-rate: "
            f"{100.0 * iid_min:.2f}% / "
            f"{100.0 * np.median(iid_rates):.2f}% / "
            f"{100.0 * np.max(iid_rates):.2f}%"
        ),
        "",
        "## Block Bootstrap Stress",
        (
            f"- min/median/max pass-rate: "
            f"{100.0 * block_min:.2f}% / "
            f"{100.0 * np.median(block_rates):.2f}% / "
            f"{100.0 * np.max(block_rates):.2f}%"
        ),
        "",
        "## Null Topology Test",
        f"- deterministic pass-rate under null topology: {100.0 * null_det_pass_rate:.2f}%",
        "",
        "## Required Next Step",
        f"- {nxt}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1972] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1972] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1972] verdict={v}")


if __name__ == "__main__":
    main()
