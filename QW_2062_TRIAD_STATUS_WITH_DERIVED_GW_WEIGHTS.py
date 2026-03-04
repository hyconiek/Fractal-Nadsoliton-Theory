#!/usr/bin/env python3
"""
QW-2062: Triad status with derived GW weights (no scan).

Combines:
- mass+flavor from QW-2060 locked shared flavor basis,
- deterministic GW weights derived from kernel invariants (no rescan),
to check whether physical thresholds are jointly passed.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2062_triad_status_with_derived_gw_weights.json"
OUT_MD = ROOT / "RAPORT_QW2062_TRIAD_STATUS_WITH_DERIVED_GW_WEIGHTS.md"


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


def derived_gw_weights(kernel: Dict[str, float]) -> Dict[str, float]:
    # Deterministic map from kernel invariants (no optimization).
    d = np.arange(1.0, 13.0, dtype=float)
    kv = np.abs(kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"]))
    decay_ratio = float(kv[3] / max(kv[0], 1e-15))

    w_max = 0.0 if decay_ratio < 0.2 else 0.05
    w_mean = float(np.clip(0.55 + 0.20 * np.tanh((1.0 - decay_ratio) - 0.4), 0.4, 0.8))
    w_0ms = float(np.clip(0.15 + 0.10 * np.tanh(abs(kernel["phi"]) + kernel["omega"]), 0.05, 0.35))
    w_10ms = max(0.0, 1.0 - (w_max + w_mean + w_0ms))

    w = np.array([w_max, w_mean, w_0ms, w_10ms], dtype=float)
    w = w / max(np.sum(w), 1e-12)
    return {
        "w_max_abs_corr": float(w[0]),
        "w_mean_abs_corr": float(w[1]),
        "w_corr_at_0ms": float(w[2]),
        "w_corr_at_10ms": float(w[3]),
        "decay_ratio": decay_ratio,
    }


def gw_metrics(weights: Dict[str, float]) -> Dict[str, float]:
    df = pd.read_csv(ROOT / "gw1831_window_features.csv")
    score = (
        weights["w_max_abs_corr"] * df["max_abs_corr"].to_numpy(dtype=float)
        + weights["w_mean_abs_corr"] * df["mean_abs_corr"].to_numpy(dtype=float)
        + weights["w_corr_at_0ms"] * df["corr_at_0ms"].to_numpy(dtype=float)
        + weights["w_corr_at_10ms"] * df["corr_at_10ms"].to_numpy(dtype=float)
    )
    pair = df["pair"].astype(str).to_numpy()
    s_hl = score[pair == "H1-L1"]
    s_hv = score[pair == "H1-V1"]
    s_lv = score[pair == "L1-V1"]
    s_ctrl = np.concatenate([s_hv, s_lv])
    q90 = float(np.quantile(s_ctrl, 0.90))
    return {
        "auc_h1l1_vs_ctrl": float(rank_auc_pos_gt_neg(s_hl, s_ctrl)),
        "adv_shared_minus_ctrl_q90": float(np.mean(s_hl > q90) - np.mean(s_ctrl > q90)),
        "sep_median_h1l1_minus_ctrl": float(np.median(s_hl) - np.median(s_ctrl)),
        "control_median_gap": float(abs(np.median(s_hv) - np.median(s_lv))),
    }


def main() -> None:
    r2060 = json.loads((ROOT / "report_qw2060_locked_shared_flavor_basis_no_rescan_gate.json").read_text(encoding="utf-8"))
    r2049 = json.loads((ROOT / "report_qw2049_spectral_micro_stagec_intersection_gate.json").read_text(encoding="utf-8"))
    kernel = {k: float(v) for k, v in r2049["stagec_pool"]["selected_kernel"].items()}

    mass = r2060["metrics"]["mass"]
    flavor = r2060["metrics"]["flavor"]

    w = derived_gw_weights(kernel)
    gw = gw_metrics(w)

    thresholds = {
        "mass_mean_rel_pct_max": 15.0,
        "mass_max_rel_pct_max": 35.0,
        "tau_charm_ratio_rel_err_pct_max": 20.0,
        "ckm_mean_rel_pct_max": 15.0,
        "pmns_mean_rel_pct_max": 15.0,
        "gw_sep_min": 0.0020,
        "gw_adv_min": 0.30,
        "gw_auc_min": 0.75,
        "gw_control_gap_max": 0.0025,
    }

    flags = {
        "mass_mean_rel_pct_le_max": bool(mass["mean_rel_err_pct"] <= thresholds["mass_mean_rel_pct_max"]),
        "mass_max_rel_pct_le_max": bool(mass["max_rel_err_pct"] <= thresholds["mass_max_rel_pct_max"]),
        "mass_tau_charm_ratio_err_le_max": bool(mass["tau_charm_ratio_rel_err_pct"] <= thresholds["tau_charm_ratio_rel_err_pct_max"]),
        "ckm_mean_rel_pct_le_max": bool(flavor["ckm_mean_rel_pct"] <= thresholds["ckm_mean_rel_pct_max"]),
        "pmns_mean_rel_pct_le_max": bool(flavor["pmns_mean_rel_pct"] <= thresholds["pmns_mean_rel_pct_max"]),
        "gw_sep_ge_min": bool(gw["sep_median_h1l1_minus_ctrl"] >= thresholds["gw_sep_min"]),
        "gw_adv_ge_min": bool(gw["adv_shared_minus_ctrl_q90"] >= thresholds["gw_adv_min"]),
        "gw_auc_ge_min": bool(gw["auc_h1l1_vs_ctrl"] >= thresholds["gw_auc_min"]),
        "gw_control_gap_le_max": bool(gw["control_median_gap"] <= thresholds["gw_control_gap_max"]),
        "gw_weights_derived_no_scan": True,
        "strict_first_principles_from_kernel_only": False,
    }

    pass_count = int(sum(1 for v in flags.values() if v))
    total_flags = int(len(flags))
    physical_flags = {k: v for k, v in flags.items() if not k.startswith("strict_") and not k.endswith("no_scan")}
    physical_pass = bool(all(physical_flags.values()))

    verdict = (
        "TRIAD_PHYSICAL_THRESHOLDS_PASS_WITH_LOCKED_FLAVOR_AND_DERIVED_GW_WEIGHTS"
        if physical_pass
        else "TRIAD_PHYSICAL_THRESHOLDS_FAIL"
    )
    required_next = (
        "DERIVE_LOCKED_FLAVOR_BASIS_FROM_KERNEL_INVARIANTS_TO_PROMOTE_TO_STRICT_FIRST_PRINCIPLES"
        if physical_pass
        else "REPAIR_REMAINING_PHYSICAL_BLOCKERS"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw2049_spectral_micro_stagec_intersection_gate.json:stagec_pool.selected_kernel",
        "kernel": kernel,
        "context_source": "report_qw2060_locked_shared_flavor_basis_no_rescan_gate.json",
        "locked_flavor_basis": r2060["locked_flavor_basis"],
        "derived_gw_weights": w,
        "metrics": {"mass": mass, "flavor": flavor, "gw": gw},
        "thresholds": thresholds,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "physical_pass": physical_pass,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2062: TRIAD STATUS WITH DERIVED GW WEIGHTS",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {pass_count}/{total_flags}",
        f"- physical_pass: {physical_pass}",
        "",
        "## Metrics",
        (
            f"- mass mean/max/tau-charm rel%: "
            f"{mass['mean_rel_err_pct']:.3f}/{mass['max_rel_err_pct']:.3f}/{mass['tau_charm_ratio_rel_err_pct']:.3f}"
        ),
        f"- flavor CKM/PMNS mean rel%: {flavor['ckm_mean_rel_pct']:.3f}/{flavor['pmns_mean_rel_pct']:.3f}",
        (
            f"- GW auc/adv/sep/gap: "
            f"{gw['auc_h1l1_vs_ctrl']:.4f}/{gw['adv_shared_minus_ctrl_q90']:.4f}/"
            f"{gw['sep_median_h1l1_minus_ctrl']:.6f}/{gw['control_median_gap']:.6f}"
        ),
        "",
        "## Derived GW Weights",
        f"- {w}",
        "",
        "## Flags",
    ]
    for k, v in flags.items():
        lines.append(f"- {k}: {v}")
    lines.extend(
        [
            "",
            "## Required Next Step",
            f"- {required_next}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2062] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2062] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2062] verdict={verdict} physical_pass={physical_pass}")


if __name__ == "__main__":
    main()
