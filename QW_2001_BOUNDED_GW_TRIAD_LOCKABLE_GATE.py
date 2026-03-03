#!/usr/bin/env python3
"""
QW-2001: Bounded-GW triad lockable gate.

Uses:
- frozen kernel + shared mass/flavor branch from QW-1967/1969,
- bounded GW operator from QW-2000,
- no sector retune.

Goal:
- verify whether triad is now lockable under bootstrap and local robustness,
- compare against QW-1970 baseline bootstrap rate.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

from QW_1983_FOLD_ROBUST_OPERATOR_SEARCH import ROOT
from QW_1993_GLOBAL_SCORE_COMPRESSION_HARD_LOCK_SEARCH import build_components
from QW_1999_BOUNDED_COUPLING_FOLD2_GUARDED_SEARCH import score_bounded


IN_QW1967 = ROOT / "report_qw1967_isospin_split_local_refinement_gate.json"
IN_QW1969 = ROOT / "report_qw1969_bootstrap_robust_recenter_search.json"
IN_QW1970 = ROOT / "report_qw1970_structural_gw_control_term_gate.json"
IN_QW2000 = ROOT / "report_qw2000_bounded_coupling_deep_audit.json"

OUT_JSON = ROOT / "report_qw2001_bounded_gw_triad_lockable_gate.json"
OUT_MD = ROOT / "RAPORT_QW2001_BOUNDED_GW_TRIAD_LOCKABLE_GATE.md"


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


def triad_flags(mass_flags: Dict[str, bool], flavor: Dict[str, float], gw: Dict[str, float], thr: Dict[str, float]) -> Dict[str, bool]:
    return {
        **mass_flags,
        "ckm_mean_rel_pct_le_max": bool(flavor["ckm_mean_rel_pct"] <= thr["ckm_mean_rel_pct_max"]),
        "pmns_mean_rel_pct_le_max": bool(flavor["pmns_mean_rel_pct"] <= thr["pmns_mean_rel_pct_max"]),
        "gw_sep_ge_min": bool(gw["sep_median_h1l1_minus_ctrl"] >= thr["gw_sep_min"]),
        "gw_adv_ge_min": bool(gw["adv_shared_minus_ctrl_q90"] >= thr["gw_adv_min"]),
        "gw_auc_ge_min": bool(gw["auc_h1l1_vs_ctrl"] >= thr["gw_auc_min"]),
        "gw_control_gap_le_max": bool(gw["control_median_gap"] <= thr["gw_control_gap_max"]),
    }


def bootstrap_triad_pass_rate(
    s_hl: np.ndarray,
    s_hv: np.ndarray,
    s_lv: np.ndarray,
    n_boot: int,
    thr: Dict[str, float],
    flavor: Dict[str, float],
    mass_flags: Dict[str, bool],
    seed: int,
) -> float:
    rng = np.random.default_rng(seed)
    n_hl, n_hv, n_lv = len(s_hl), len(s_hv), len(s_lv)
    pass_count = 0

    flavor_flags = {
        "ckm_mean_rel_pct_le_max": bool(flavor["ckm_mean_rel_pct"] <= thr["ckm_mean_rel_pct_max"]),
        "pmns_mean_rel_pct_le_max": bool(flavor["pmns_mean_rel_pct"] <= thr["pmns_mean_rel_pct_max"]),
    }

    for _ in range(n_boot):
        b_hl = s_hl[rng.integers(0, n_hl, size=n_hl, endpoint=False)]
        b_hv = s_hv[rng.integers(0, n_hv, size=n_hv, endpoint=False)]
        b_lv = s_lv[rng.integers(0, n_lv, size=n_lv, endpoint=False)]
        g = gw_metrics_from_scores(b_hl, b_hv, b_lv)
        flags = {
            **mass_flags,
            **flavor_flags,
            "gw_sep_ge_min": bool(g["sep_median_h1l1_minus_ctrl"] >= thr["gw_sep_min"]),
            "gw_adv_ge_min": bool(g["adv_shared_minus_ctrl_q90"] >= thr["gw_adv_min"]),
            "gw_auc_ge_min": bool(g["auc_h1l1_vs_ctrl"] >= thr["gw_auc_min"]),
            "gw_control_gap_le_max": bool(g["control_median_gap"] <= thr["gw_control_gap_max"]),
        }
        if all(flags.values()):
            pass_count += 1
    return float(pass_count / n_boot)


def local_deterministic_pass_rate(
    base: np.ndarray,
    c1: np.ndarray,
    c3: np.ndarray,
    c4: np.ndarray,
    pairs: np.ndarray,
    center: Dict[str, float],
    thr: Dict[str, float],
    flavor: Dict[str, float],
    mass_flags: Dict[str, bool],
    radius: float,
    n_samples: int,
    seed: int,
) -> float:
    # Parameter bounds inherited from bounded-GW search family.
    bounds = {
        "xi1": (0.0001, 0.0032),
        "xi3": (-0.0002, 0.0032),
        "xi4": (-0.0032, 0.0032),
        "gamma_c": (0.75, 1.02),
        "kappa_t": (0.90, 2.20),
    }

    rng = np.random.default_rng(seed)
    pass_count = 0
    for _ in range(n_samples):
        p = {}
        for k, (lo, hi) in bounds.items():
            span = hi - lo
            v = float(center[k]) + float(rng.uniform(-radius, radius) * span)
            p[k] = float(min(hi, max(lo, v)))

        s = score_bounded(base, c1, c3, c4, p["xi1"], p["xi3"], p["xi4"], p["gamma_c"], p["kappa_t"])
        s_hl = s[pairs == 0]
        s_hv = s[pairs == 1]
        s_lv = s[pairs == 2]
        g = gw_metrics_from_scores(s_hl, s_hv, s_lv)
        flags = triad_flags(mass_flags=mass_flags, flavor=flavor, gw=g, thr=thr)
        if all(flags.values()):
            pass_count += 1

    return float(pass_count / n_samples)


def verdict(boot5000_min: float, local_rates: Dict[str, float], deterministic_all_pass: bool) -> str:
    if (
        deterministic_all_pass
        and boot5000_min >= 0.95
        and local_rates["r01"] >= 0.95
        and local_rates["r02"] >= 0.90
        and local_rates["r05"] >= 0.70
    ):
        return "BOUNDED_GW_TRIAD_LOCKABLE_PASS"
    if deterministic_all_pass and boot5000_min >= 0.80 and local_rates["r01"] >= 0.90 and local_rates["r02"] >= 0.80:
        return "BOUNDED_GW_TRIAD_PARTIAL_ROBUSTNESS"
    return "BOUNDED_GW_TRIAD_NOT_LOCKABLE"


def main() -> None:
    r1967 = json.loads(IN_QW1967.read_text(encoding="utf-8"))
    r1969 = json.loads(IN_QW1969.read_text(encoding="utf-8"))
    r1970 = json.loads(IN_QW1970.read_text(encoding="utf-8"))
    r2000 = json.loads(IN_QW2000.read_text(encoding="utf-8"))

    kernel = r1969["frozen_kernel"]
    params = r1969["best_recentered_candidate"]["params"]
    thr = r1969["thresholds"]
    mass_flags = r1967["source_mass_branch"]["mass_flags_fixed"]
    flavor = r1967["summary"]["best"]["flavor"]

    c = r2000["candidate"]
    center = {
        "xi1": float(c["xi1"]),
        "xi3": float(c["xi3"]),
        "xi4": float(c["xi4"]),
        "gamma_c": float(c["gamma_c"]),
        "kappa_t": float(c["kappa_t"]),
    }

    df = pd.read_csv(ROOT / "gw1831_window_features.csv")
    base, c1, c3, c4, pairs = build_components(df, kernel, params)

    s = score_bounded(base, c1, c3, c4, center["xi1"], center["xi3"], center["xi4"], center["gamma_c"], center["kappa_t"])
    s_hl = s[pairs == 0]
    s_hv = s[pairs == 1]
    s_lv = s[pairs == 2]
    gw = gw_metrics_from_scores(s_hl, s_hv, s_lv)
    flags = triad_flags(mass_flags=mass_flags, flavor=flavor, gw=gw, thr=thr)

    boot2500 = [
        bootstrap_triad_pass_rate(s_hl, s_hv, s_lv, n_boot=2500, thr=thr, flavor=flavor, mass_flags=mass_flags, seed=200100 + i)
        for i in range(3)
    ]
    boot5000 = [
        bootstrap_triad_pass_rate(s_hl, s_hv, s_lv, n_boot=5000, thr=thr, flavor=flavor, mass_flags=mass_flags, seed=200200 + i)
        for i in range(5)
    ]
    boot10000 = [
        bootstrap_triad_pass_rate(
            s_hl,
            s_hv,
            s_lv,
            n_boot=10000,
            thr=thr,
            flavor=flavor,
            mass_flags=mass_flags,
            seed=200300,
        )
    ]

    local_rates = {
        "r01": local_deterministic_pass_rate(
            base,
            c1,
            c3,
            c4,
            pairs,
            center,
            thr,
            flavor,
            mass_flags,
            radius=0.01,
            n_samples=800,
            seed=200401,
        ),
        "r02": local_deterministic_pass_rate(
            base,
            c1,
            c3,
            c4,
            pairs,
            center,
            thr,
            flavor,
            mass_flags,
            radius=0.02,
            n_samples=800,
            seed=200402,
        ),
        "r05": local_deterministic_pass_rate(
            base,
            c1,
            c3,
            c4,
            pairs,
            center,
            thr,
            flavor,
            mass_flags,
            radius=0.05,
            n_samples=800,
            seed=200405,
        ),
    }

    boot5000_min = float(np.min(np.array(boot5000, dtype=float)))
    deterministic_all_pass = bool(all(flags.values()))
    v = verdict(boot5000_min=boot5000_min, local_rates=local_rates, deterministic_all_pass=deterministic_all_pass)

    baseline_boot5000 = float(r1970["best"]["boot_pass_rate_5000"])

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": [
            IN_QW1967.name,
            IN_QW1969.name,
            IN_QW1970.name,
            IN_QW2000.name,
        ],
        "frozen_components": {
            "kernel": kernel,
            "shared_params": params,
            "gw_operator": center,
        },
        "thresholds": thr,
        "deterministic": {
            "flavor": flavor,
            "gw": gw,
            "flags": flags,
            "all_pass": deterministic_all_pass,
        },
        "bootstrap": {
            "triad_pass_rate_2500_seeds": boot2500,
            "triad_pass_rate_5000_seeds": boot5000,
            "triad_pass_rate_10000_seeds": boot10000,
            "triad_pass_rate_5000_min": boot5000_min,
            "triad_pass_rate_5000_mean": float(np.mean(np.array(boot5000, dtype=float))),
        },
        "local_deterministic_pass_rates": local_rates,
        "baseline_comparison": {
            "q1970_boot_pass_rate_5000": baseline_boot5000,
            "improvement_abs": float(np.mean(np.array(boot5000, dtype=float)) - baseline_boot5000),
        },
        "verdict": v,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2001: BOUNDED GW TRIAD LOCKABLE GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{v}**",
        "",
        "## Deterministic Triad",
        f"- all_pass: {deterministic_all_pass}",
        f"- ckm_mean_rel_pct: {flavor['ckm_mean_rel_pct']:.3f}%",
        f"- pmns_mean_rel_pct: {flavor['pmns_mean_rel_pct']:.3f}%",
        f"- gw auc/adv/sep/gap: {gw['auc_h1l1_vs_ctrl']:.4f} / {gw['adv_shared_minus_ctrl_q90']:.4f} / {gw['sep_median_h1l1_minus_ctrl']:.6f} / {gw['control_median_gap']:.6f}",
        "",
        "## Bootstrap Triad Pass Rate",
        f"- 2500 seeds: {', '.join(f'{x:.4f}' for x in boot2500)}",
        f"- 5000 seeds: {', '.join(f'{x:.4f}' for x in boot5000)}",
        f"- 10000 seed: {', '.join(f'{x:.4f}' for x in boot10000)}",
        f"- 5000 min/mean: {boot5000_min:.4f} / {np.mean(np.array(boot5000, dtype=float)):.4f}",
        "",
        "## Local Deterministic Pass Rate",
        f"- radius 1%: {local_rates['r01']:.4f}",
        f"- radius 2%: {local_rates['r02']:.4f}",
        f"- radius 5%: {local_rates['r05']:.4f}",
        "",
        "## Baseline Comparison",
        f"- QW-1970 boot_pass_rate_5000: {baseline_boot5000:.4f}",
        f"- Improvement abs: {out['baseline_comparison']['improvement_abs']:.4f}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2001] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2001] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2001] verdict={v}")


if __name__ == "__main__":
    main()
