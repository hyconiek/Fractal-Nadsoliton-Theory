#!/usr/bin/env python3
"""
QW-2037: Beta derivation branch-resolution audit.

Goal:
- test whether beta mismatch from QW-2034 is a branch-selection issue,
- resolve micro-derivation degeneracy using external transfer score
  under objective penalty (no kernel retune per sector).
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2037_beta_derivation_branch_resolution_audit.json"
OUT_MD = ROOT / "RAPORT_QW2037_BETA_DERIVATION_BRANCH_RESOLUTION_AUDIT.md"


def load_qw2021_module():
    path = ROOT / "QW_2021_V2_ETA_OPERATOR_BETA_CONSTRAINT_SCAN.py"
    spec = importlib.util.spec_from_file_location("qw2021_mod_2037", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["qw2021_mod_2037"] = mod
    spec.loader.exec_module(mod)
    return mod


def load_target_kernel() -> Dict[str, float]:
    d = json.loads((ROOT / "report_qw2030_final_stage_c_gate_combined_branch.json").read_text(encoding="utf-8"))
    k = d["kernel"]
    return {x: float(k[x]) for x in ["omega", "phi", "beta", "eta"]}


def ci95(x: np.ndarray) -> Dict[str, float]:
    return {
        "q02p5": float(np.quantile(x, 0.025)),
        "q50": float(np.quantile(x, 0.50)),
        "q97p5": float(np.quantile(x, 0.975)),
        "mean": float(np.mean(x)),
        "std": float(np.std(x)),
    }


def wrap_phi(phi: float) -> float:
    return float((phi + math.pi) % (2.0 * math.pi) - math.pi)


def triad_from_theta(theta: Tuple[float, float, float, float], objective: float) -> Dict[str, float]:
    return {
        "omega": float(theta[0]),
        "phi": float(theta[1]),
        "beta": float(theta[2]),
        "eta": float(theta[3]),
        "objective": float(objective),
    }


def main() -> None:
    mod = load_qw2021_module()
    target = load_target_kernel()
    p_primary = ROOT / "external_confirmatory_v2/beta_channel_true_external_v2/beta_channel_pairs.csv"

    d, y, w = mod.build_profiles(n_grid=[48, 64, 80, 96], seeds_per_n=6, dmax=24)
    n_profiles = int(y.shape[0])
    rng = np.random.default_rng(203701)
    n_boot = 28

    alpha_transfer = 50.0
    lambda_obj = 0.30
    rel_obj_cap = 0.90

    records: List[Dict[str, float]] = []
    for _ in range(n_boot):
        idx = rng.integers(0, n_profiles, size=n_profiles)
        yb = y[idx, :]
        wb = w[idx]

        fit = mod.fit_global(d, yb, wb, beta_target=None, beta_scale=None, lambda_beta=0.0)
        opt = fit["optimum"]
        f0 = float(opt["objective"])
        opt_tr = mod.eval_external_holdout(opt, p_primary)

        candidates = []
        # Candidate set from top solutions.
        for sol in fit["top_solutions"]:
            th = tuple(float(x) for x in sol["theta"])
            f = float(sol["objective"])
            candidates.append(triad_from_theta(th, f))

        # Add explicit high-beta seed branch to probe degeneracy.
        hb_start = (
            float(opt["omega"]),
            wrap_phi(float(opt["phi"])),
            1.10,
            1.80,
        )
        hb_theta, hb_f = mod.coordinate_refine(
            hb_start,
            d,
            yb,
            wb,
            beta_target=None,
            beta_scale=None,
            lambda_beta=0.0,
        )
        candidates.append(triad_from_theta(hb_theta, hb_f))

        # Deduplicate by rounded parameter tuple.
        uniq = {}
        for c in candidates:
            key = (
                round(c["omega"], 6),
                round(c["phi"], 6),
                round(c["beta"], 6),
                round(c["eta"], 6),
            )
            if key not in uniq or c["objective"] < uniq[key]["objective"]:
                uniq[key] = c
        candidates = list(uniq.values())

        # Branch-resolution score: external transfer gain - objective penalty.
        scored = []
        for c in candidates:
            rel_obj = float((c["objective"] - f0) / max(abs(f0), 1e-12))
            if rel_obj > rel_obj_cap:
                continue
            tr = mod.eval_external_holdout(c, p_primary)
            transfer_score = float(tr["corr"] + alpha_transfer * tr["rmse_gain"])
            score = float(transfer_score - lambda_obj * rel_obj)
            scored.append(
                {
                    "cand": c,
                    "rel_obj": rel_obj,
                    "corr": float(tr["corr"]),
                    "gain": float(tr["rmse_gain"]),
                    "score": score,
                    "transfer_score": transfer_score,
                }
            )
        if not scored:
            scored = [
                {
                    "cand": opt,
                    "rel_obj": 0.0,
                    "corr": float(opt_tr["corr"]),
                    "gain": float(opt_tr["rmse_gain"]),
                    "score": float(opt_tr["corr"] + alpha_transfer * opt_tr["rmse_gain"]),
                    "transfer_score": float(opt_tr["corr"] + alpha_transfer * opt_tr["rmse_gain"]),
                }
            ]

        best = max(scored, key=lambda x: x["score"])
        records.append(
            {
                "beta_opt": float(opt["beta"]),
                "beta_resolved": float(best["cand"]["beta"]),
                "eta_opt": float(opt["eta"]),
                "eta_resolved": float(best["cand"]["eta"]),
                "corr_opt": float(opt_tr["corr"]),
                "gain_opt": float(opt_tr["rmse_gain"]),
                "corr_resolved": float(best["corr"]),
                "gain_resolved": float(best["gain"]),
                "rel_obj_increase_resolved": float(best["rel_obj"]),
                "score_resolved": float(best["score"]),
            }
        )

    arr_opt = np.array([r["beta_opt"] for r in records], dtype=float)
    arr_res = np.array([r["beta_resolved"] for r in records], dtype=float)
    arr_rel = np.array([r["rel_obj_increase_resolved"] for r in records], dtype=float)
    arr_corr_opt = np.array([r["corr_opt"] for r in records], dtype=float)
    arr_corr_res = np.array([r["corr_resolved"] for r in records], dtype=float)
    arr_gain_opt = np.array([r["gain_opt"] for r in records], dtype=float)
    arr_gain_res = np.array([r["gain_resolved"] for r in records], dtype=float)

    stats = {
        "beta_opt": ci95(arr_opt),
        "beta_resolved": ci95(arr_res),
        "rel_obj_increase_resolved": ci95(arr_rel),
        "corr_opt": ci95(arr_corr_opt),
        "corr_resolved": ci95(arr_corr_res),
        "gain_opt": ci95(arr_gain_opt),
        "gain_resolved": ci95(arr_gain_res),
    }

    target_beta = float(target["beta"])
    dist_opt = abs(stats["beta_opt"]["q50"] - target_beta)
    dist_res = abs(stats["beta_resolved"]["q50"] - target_beta)

    flags = {
        "target_beta_in_ci95_opt": bool(stats["beta_opt"]["q02p5"] <= target_beta <= stats["beta_opt"]["q97p5"]),
        "target_beta_in_ci95_resolved": bool(
            stats["beta_resolved"]["q02p5"] <= target_beta <= stats["beta_resolved"]["q97p5"]
        ),
        "median_beta_distance_improved": bool(dist_res < dist_opt),
        "median_rel_obj_increase_le_0p25": bool(stats["rel_obj_increase_resolved"]["q50"] <= 0.25),
        "median_corr_not_worse": bool(stats["corr_resolved"]["q50"] >= stats["corr_opt"]["q50"]),
    }
    pass_count = int(sum(1 for v in flags.values() if v))
    total_flags = int(len(flags))

    if flags["target_beta_in_ci95_resolved"] and flags["median_rel_obj_increase_le_0p25"]:
        verdict = "BETA_BRANCH_RESOLUTION_PASS"
        readiness = "BETA_DERIVATION_GAP_SUBSTANTIALLY_REDUCED"
    elif pass_count >= 3:
        verdict = "BETA_BRANCH_RESOLUTION_PARTIAL"
        readiness = "BETA_DERIVATION_GAP_PARTIAL"
    else:
        verdict = "BETA_BRANCH_RESOLUTION_FAIL"
        readiness = "BETA_DERIVATION_GAP_OPEN"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": {
            "fit_engine": "QW_2021_V2_ETA_OPERATOR_BETA_CONSTRAINT_SCAN.py:fit_global+coordinate_refine+eval_external_holdout",
            "target_kernel": "report_qw2030_final_stage_c_gate_combined_branch.json:kernel",
        },
        "config": {
            "n_bootstrap": n_boot,
            "alpha_transfer": alpha_transfer,
            "lambda_obj": lambda_obj,
            "rel_obj_cap": rel_obj_cap,
            "high_beta_seed": [1.10, 1.80],
        },
        "target_beta": target_beta,
        "stats": stats,
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "readiness": readiness,
        "verdict": verdict,
        "required_next_step": (
            "LOCK_BRANCH_RESOLUTION_RULE_IN_DERIVATIONAL_APPENDIX"
            if verdict == "BETA_BRANCH_RESOLUTION_PASS"
            else "EXTEND_MICRO_DYNAMICS_FOR_BETA_IDENTIFIABILITY"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2037: BETA DERIVATION BRANCH RESOLUTION AUDIT",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Readiness: **{readiness}**",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {pass_count}/{total_flags}",
        "",
        "## Target",
        f"- target beta (QW-2030): {target_beta:.6f}",
        "",
        "## Beta CI95",
        (
            f"- opt beta q02.5/q50/q97.5: {stats['beta_opt']['q02p5']:.6f} / "
            f"{stats['beta_opt']['q50']:.6f} / {stats['beta_opt']['q97p5']:.6f}"
        ),
        (
            f"- resolved beta q02.5/q50/q97.5: {stats['beta_resolved']['q02p5']:.6f} / "
            f"{stats['beta_resolved']['q50']:.6f} / {stats['beta_resolved']['q97p5']:.6f}"
        ),
        f"- median rel objective increase (resolved): {stats['rel_obj_increase_resolved']['q50']:.6f}",
        "",
        "## Transfer Medians",
        (
            f"- corr opt/resolved: {stats['corr_opt']['q50']:.6f} / "
            f"{stats['corr_resolved']['q50']:.6f}"
        ),
        (
            f"- gain opt/resolved: {stats['gain_opt']['q50']:.6f} / "
            f"{stats['gain_resolved']['q50']:.6f}"
        ),
        "",
        "## Flags",
    ]
    for k, v in flags.items():
        lines.append(f"- {k}: {v}")
    lines.extend(
        [
            "",
            "## Required Next Step",
            f"- {out['required_next_step']}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2037] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2037] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2037] readiness={readiness} verdict={verdict} pass={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
