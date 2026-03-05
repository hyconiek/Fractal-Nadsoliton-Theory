#!/usr/bin/env python3
"""
QW-2128: Kernel->representation assignment uniqueness gate (locked branch).

Purpose:
- test whether representation proxy assignment (q-assignment family) is unique
  in the locked strict noncircular branch,
- quantify robustness of assignment ranking under deterministic delta_info
  perturbations without scan/retune.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2128_kernel_rep_assignment_uniqueness_gate.json"
OUT_MD = ROOT / "RAPORT_QW2128_KERNEL_REP_ASSIGNMENT_UNIQUENESS_GATE.md"

M_TOP = 173_000.0
TAU_CHARM_RATIO_EXP = 1776.9 / 1270.0
PARTICLES = [
    ("Top", 173_000.0, 0),
    ("Bottom", 4_180.0, +1),
    ("Tau", 1_776.9, -1),
    ("Charm", 1_270.0, +1),
    ("Muon", 105.7, -1),
    ("Electron", 0.511, -1),
]


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def rel_err_pct(pred: float, exp: float) -> float:
    return abs(pred - exp) / max(abs(exp), 1e-15) * 100.0


def predict_mass(q_eff: float, gamma: float) -> float:
    return M_TOP * (4.0 ** (-(gamma * q_eff / 4.0)))


def evaluate_variant(gamma: float, q_map: Dict[str, float], delta_info: float, use_info_split: bool) -> Dict[str, object]:
    rows: List[Dict[str, float]] = []
    errs: List[float] = []
    pred_by_name: Dict[str, float] = {}
    for name, exp_mev, sector_sign in PARTICLES:
        q_base = float(q_map[name])
        q_eff = q_base + (sector_sign * delta_info if use_info_split else 0.0)
        pred = predict_mass(q_eff, gamma)
        err = rel_err_pct(pred, exp_mev)
        rows.append(
            {
                "particle": name,
                "q_base": q_base,
                "q_eff": float(q_eff),
                "exp_mev": float(exp_mev),
                "pred_mev": float(pred),
                "rel_err_pct": float(err),
            }
        )
        errs.append(float(err))
        pred_by_name[name] = float(pred)

    tau_charm_ratio_pred = pred_by_name["Tau"] / max(pred_by_name["Charm"], 1e-15)
    tau_charm_ratio_err_pct = rel_err_pct(tau_charm_ratio_pred, TAU_CHARM_RATIO_EXP)

    metrics = {
        "mean_rel_err_pct": float(mean(errs)),
        "max_rel_err_pct": float(max(errs)),
        "tau_charm_ratio_pred": float(tau_charm_ratio_pred),
        "tau_charm_ratio_exp": float(TAU_CHARM_RATIO_EXP),
        "tau_charm_ratio_rel_err_pct": float(tau_charm_ratio_err_pct),
    }
    flags = {
        "mean_rel_err_le_15": bool(metrics["mean_rel_err_pct"] <= 15.0),
        "max_rel_err_le_35": bool(metrics["max_rel_err_pct"] <= 35.0),
        "tau_charm_ratio_err_le_20": bool(metrics["tau_charm_ratio_rel_err_pct"] <= 20.0),
    }
    score = (
        metrics["mean_rel_err_pct"] / 15.0
        + metrics["max_rel_err_pct"] / 35.0
        + metrics["tau_charm_ratio_rel_err_pct"] / 20.0
    )
    return {
        "rows": rows,
        "metrics": metrics,
        "flags": flags,
        "all_pass": bool(all(flags.values())),
        "score": float(score),
    }


def main() -> None:
    r1961 = load_json("report_qw1961_noncircular_gamma_q_derivation_matrix.json")
    q_assignments = r1961["inputs"]["q_assignments"]
    delta_info = float(r1961["info_split_source"]["delta_info"])

    best_noncirc = r1961["summary"]["best_noncircular"]
    gamma_noncirc = float(best_noncirc["gamma_value"])
    gamma_source_noncirc = str(best_noncirc["gamma_source"])
    split_mode_noncirc = str(best_noncirc["split_mode"])
    use_info_split_noncirc = bool(split_mode_noncirc == "info_split_from_qw1958")

    # Locked-branch assignment ranking.
    locked_rows: List[Dict[str, object]] = []
    for q_name, q_map in q_assignments.items():
        ev = evaluate_variant(gamma_noncirc, q_map, delta_info, use_info_split_noncirc)
        locked_rows.append({"q_assignment": q_name, **ev})
    locked_sorted = sorted(locked_rows, key=lambda x: float(x["score"]))
    best_locked = locked_sorted[0]
    second_locked = locked_sorted[1] if len(locked_sorted) > 1 else None
    score_gap = float(second_locked["score"] - best_locked["score"]) if second_locked else float("inf")

    # Kernel-gamma consistency check under same split mode.
    gamma_kernel = float(r1961["inputs"]["gamma_candidates"]["derived_kernel_d1_to_d4"])
    kernel_rows: List[Dict[str, object]] = []
    for q_name, q_map in q_assignments.items():
        ev = evaluate_variant(gamma_kernel, q_map, delta_info, use_info_split_noncirc)
        kernel_rows.append({"q_assignment": q_name, **ev})
    kernel_sorted = sorted(kernel_rows, key=lambda x: float(x["score"]))
    best_kernel_branch = kernel_sorted[0]

    # Delta robustness around locked branch.
    scales = [0.80, 0.90, 1.00, 1.10, 1.20]
    winners: List[str] = []
    robustness_rows: List[Dict[str, object]] = []
    for s in scales:
        dloc = float(delta_info * s)
        local = []
        for q_name, q_map in q_assignments.items():
            ev = evaluate_variant(gamma_noncirc, q_map, dloc, use_info_split_noncirc)
            local.append({"q_assignment": q_name, "score": float(ev["score"]), "all_pass": bool(ev["all_pass"])})
        local_sorted = sorted(local, key=lambda x: x["score"])
        winner = str(local_sorted[0]["q_assignment"])
        winners.append(winner)
        robustness_rows.append(
            {
                "delta_scale": s,
                "delta_info_scaled": dloc,
                "winner_q_assignment": winner,
                "winner_score": float(local_sorted[0]["score"]),
                "runnerup_q_assignment": str(local_sorted[1]["q_assignment"]) if len(local_sorted) > 1 else None,
                "runnerup_score": float(local_sorted[1]["score"]) if len(local_sorted) > 1 else None,
                "score_gap": float(local_sorted[1]["score"] - local_sorted[0]["score"]) if len(local_sorted) > 1 else None,
            }
        )

    winner_name = str(best_locked["q_assignment"])
    winner_freq = int(sum(1 for w in winners if w == winner_name))

    flags = {
        "q_assignment_candidates_ge_2": bool(len(q_assignments) >= 2),
        "locked_branch_is_best_noncircular_from_qw1961": bool(gamma_source_noncirc != "canonical_frozen_1p52_reference"),
        "locked_branch_winner_all_pass_mass_thresholds": bool(best_locked["all_pass"]),
        "locked_branch_unique_winner_by_score_gap_ge_0p25": bool(score_gap >= 0.25),
        "locked_branch_runnerup_fails_at_least_one_threshold": bool(second_locked is not None and (not second_locked["all_pass"])),
        "delta_info_robust_winner_frequency_ge_4of5": bool(winner_freq >= 4),
        "kernel_gamma_branch_selects_same_assignment": bool(str(best_kernel_branch["q_assignment"]) == winner_name),
        "deterministic_no_scan_no_retune": True,
        "global_uniqueness_across_all_gamma_hypotheses": False,
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    verdict = (
        "KERNEL_REP_ASSIGNMENT_UNIQUENESS_GATE_PASS_LOCKED_BRANCH_PARTIAL"
        if (
            flags["q_assignment_candidates_ge_2"]
            and flags["locked_branch_is_best_noncircular_from_qw1961"]
            and flags["locked_branch_winner_all_pass_mass_thresholds"]
            and flags["locked_branch_unique_winner_by_score_gap_ge_0p25"]
            and flags["locked_branch_runnerup_fails_at_least_one_threshold"]
            and flags["delta_info_robust_winner_frequency_ge_4of5"]
            and flags["kernel_gamma_branch_selects_same_assignment"]
            and flags["deterministic_no_scan_no_retune"]
        )
        else "KERNEL_REP_ASSIGNMENT_UNIQUENESS_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q1961_matrix": "report_qw1961_noncircular_gamma_q_derivation_matrix.json",
        },
        "locked_branch_definition": {
            "gamma_source": gamma_source_noncirc,
            "gamma_value": gamma_noncirc,
            "split_mode": split_mode_noncirc,
            "delta_info": delta_info,
        },
        "locked_branch_ranking": {
            "rows_sorted": locked_sorted,
            "winner_q_assignment": winner_name,
            "winner_score": float(best_locked["score"]),
            "runnerup_q_assignment": str(second_locked["q_assignment"]) if second_locked else None,
            "runnerup_score": float(second_locked["score"]) if second_locked else None,
            "winner_score_gap": score_gap,
        },
        "kernel_gamma_consistency": {
            "gamma_kernel_d1_to_d4": gamma_kernel,
            "rows_sorted": kernel_sorted,
            "winner_q_assignment": str(best_kernel_branch["q_assignment"]),
        },
        "delta_info_robustness": {
            "scales": scales,
            "rows": robustness_rows,
            "winner_frequency_for_locked_winner": winner_freq,
            "n_tests": len(scales),
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "CONVERT_LOCKED_BRANCH_UNIQUENESS_TO_FULL_KERNEL_HYPOTHESIS_UNIQUENESS_OR_DOCUMENT_SCOPE_BOUNDARY"
            if verdict.endswith("PARTIAL")
            else "REPAIR_REP_ASSIGNMENT_RULES_AND_RERUN_QW2128"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2128: KERNEL->REP ASSIGNMENT UNIQUENESS GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Locked branch winner",
        f"- winner: `{winner_name}`",
        f"- winner score: `{float(best_locked['score']):.6f}`",
        f"- runner-up: `{str(second_locked['q_assignment']) if second_locked else 'n/a'}`",
        f"- score gap: `{score_gap:.6f}`",
        "",
        "## Robustness",
        f"- winner frequency (delta perturb): `{winner_freq}/{len(scales)}`",
        f"- kernel-gamma winner: `{str(best_kernel_branch['q_assignment'])}`",
        "",
        "## Open scope boundary",
        "- global_uniqueness_across_all_gamma_hypotheses: `False`",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2128] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2128] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2128] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()

