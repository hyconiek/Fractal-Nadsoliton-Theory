#!/usr/bin/env python3
"""
QW-2011: Strict robustness audit for non-circular mass branch (Task 1).

Goal:
- execute task #1 from closure list: non-circular mass route,
- verify deterministic + uncertainty robustness without Q<-mass inversion.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean
from typing import Dict, List, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2011_noncircular_mass_branch_strict_robustness.json"
OUT_MD = ROOT / "RAPORT_QW2011_NONCIRCULAR_MASS_BRANCH_STRICT_ROBUSTNESS.md"

PARTICLES = [
    ("Top", 173000.0, 0),
    ("Bottom", 4180.0, +1),
    ("Tau", 1776.9, -1),
    ("Charm", 1270.0, +1),
    ("Muon", 105.7, -1),
    ("Electron", 0.511, -1),
]

N_GRAV_BASE = 2.26
M_TOP = 173000.0
TAU_CHARM_RATIO_EXP = 1776.9 / 1270.0


def rel_err_pct(pred: float, exp: float) -> float:
    return abs(pred - exp) / max(exp, 1e-15) * 100.0


def predict_mass(q_eff: float, gamma: float) -> float:
    return M_TOP * (4.0 ** (-(gamma * q_eff / 4.0)))


def eval_mass(q_map: Dict[str, float], gamma: float, delta_info: float) -> Dict[str, float]:
    rows: List[Dict[str, float]] = []
    errs = []
    pred = {}
    for name, exp, sign in PARTICLES:
        q = float(q_map[name]) + sign * float(delta_info)
        m = predict_mass(q, gamma)
        e = rel_err_pct(m, exp)
        errs.append(e)
        pred[name] = m
        rows.append(
            {
                "particle": name,
                "q_eff": float(q),
                "pred_mev": float(m),
                "exp_mev": float(exp),
                "rel_err_pct": float(e),
            }
        )

    tau_charm_pred = pred["Tau"] / max(pred["Charm"], 1e-15)
    tau_charm_err = rel_err_pct(tau_charm_pred, TAU_CHARM_RATIO_EXP)

    return {
        "rows": rows,
        "mean_rel_err_pct": float(mean(errs)),
        "max_rel_err_pct": float(max(errs)),
        "tau_charm_ratio_pred": float(tau_charm_pred),
        "tau_charm_ratio_exp": float(TAU_CHARM_RATIO_EXP),
        "tau_charm_ratio_rel_err_pct": float(tau_charm_err),
    }


def main() -> None:
    d1961 = json.loads((ROOT / "report_qw1961_noncircular_gamma_q_derivation_matrix.json").read_text(encoding="utf-8"))
    d1958 = json.loads((ROOT / "report_qw1958_formal_nogo_and_minimal_lagrangian_term.json").read_text(encoding="utf-8"))

    best = d1961["summary"]["best_noncircular"]
    q_name = best["q_assignment"]
    q_map = d1961["inputs"]["q_assignments"][q_name]
    gamma_ref = float(best["gamma_value"])

    coeffs = d1958["minimal_eft_repair_term"]["deterministic_coefficients"]
    lam = float(coeffs["lambda_I"])
    mu = float(coeffs["mu_I"])
    delta_info_ref = (lam - mu) / max(lam + mu, 1e-15)

    thresholds = {
        "mean_rel_err_pct_max": 15.0,
        "max_rel_err_pct_max": 35.0,
        "tau_charm_ratio_rel_err_pct_max": 20.0,
    }

    det = eval_mass(q_map=q_map, gamma=gamma_ref, delta_info=delta_info_ref)
    det_flags = {
        "mean_rel_err_le_max": bool(det["mean_rel_err_pct"] <= thresholds["mean_rel_err_pct_max"]),
        "max_rel_err_le_max": bool(det["max_rel_err_pct"] <= thresholds["max_rel_err_pct_max"]),
        "tau_charm_ratio_err_le_max": bool(det["tau_charm_ratio_rel_err_pct"] <= thresholds["tau_charm_ratio_rel_err_pct_max"]),
    }
    det_pass = bool(all(det_flags.values()))

    # Uncertainty robustness: no fitting, only uncertainty propagation around derivational sources.
    rng = np.random.default_rng(2011)
    n_mc = 10000
    mc_rows = []
    pass_count = 0

    for _ in range(n_mc):
        n_grav = float(rng.normal(N_GRAV_BASE, 0.03))  # conservative uncertainty model
        gamma = float((2.0 * n_grav) / 3.0)

        delta_jitter = float(rng.normal(0.0, abs(delta_info_ref) * 0.10 + 1e-6))
        delta_info = float(delta_info_ref + delta_jitter)

        em = eval_mass(q_map=q_map, gamma=gamma, delta_info=delta_info)
        flags = {
            "mean_rel_err_le_max": bool(em["mean_rel_err_pct"] <= thresholds["mean_rel_err_pct_max"]),
            "max_rel_err_le_max": bool(em["max_rel_err_pct"] <= thresholds["max_rel_err_pct_max"]),
            "tau_charm_ratio_err_le_max": bool(em["tau_charm_ratio_rel_err_pct"] <= thresholds["tau_charm_ratio_rel_err_pct_max"]),
        }
        ok = bool(all(flags.values()))
        if ok:
            pass_count += 1
        mc_rows.append(
            {
                "gamma": gamma,
                "delta_info": delta_info,
                "mean_rel_err_pct": em["mean_rel_err_pct"],
                "max_rel_err_pct": em["max_rel_err_pct"],
                "tau_charm_ratio_rel_err_pct": em["tau_charm_ratio_rel_err_pct"],
                "all_pass": ok,
            }
        )

    pass_rate = pass_count / n_mc
    pass_gate = 0.95

    verdict = (
        "NONCIRCULAR_MASS_BRANCH_STRICT_ROBUST_PASS"
        if det_pass and pass_rate >= pass_gate
        else "NONCIRCULAR_MASS_BRANCH_STRICT_ROBUST_FAIL"
    )
    required_next = (
        "FREEZE_THIS_NONCIRCULAR_MASS_BRANCH_FOR_TRIAD_WITHOUT_MASS_INVERSION"
        if verdict.endswith("PASS")
        else "MASS_BRANCH_IS_NONCIRCULAR_BUT_NOT_ROBUST_ENOUGH_REWORK_DERIVATIONAL_UNCERTAINTY"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input_best_noncircular_from_q1961": {
            "q_assignment": q_name,
            "gamma_source": best["gamma_source"],
            "gamma_value": gamma_ref,
            "split_mode": best["split_mode"],
        },
        "q_assignment_used": q_map,
        "delta_info_reference": float(delta_info_ref),
        "thresholds": thresholds,
        "deterministic": {
            "metrics": det,
            "flags": det_flags,
            "all_pass": det_pass,
        },
        "uncertainty_mc": {
            "n_mc": int(n_mc),
            "n_grav_model": "Normal(2.26, 0.03)",
            "delta_info_jitter_model": "Normal(0, 0.1*|delta_info_ref|)",
            "pass_rate": float(pass_rate),
            "pass_gate": float(pass_gate),
            "samples_head": mc_rows[:40],
        },
        "verdict": verdict,
        "required_next_step": required_next,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2011: NONCIRCULAR MASS BRANCH STRICT ROBUSTNESS",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Deterministic Noncircular Branch",
        f"- q_assignment: {q_name}",
        f"- gamma: {gamma_ref:.6f} ({best['gamma_source']})",
        f"- delta_info_ref: {delta_info_ref:.6f}",
        (
            f"- mean/max/tau-charm err%: {det['mean_rel_err_pct']:.3f} / "
            f"{det['max_rel_err_pct']:.3f} / {det['tau_charm_ratio_rel_err_pct']:.3f}"
        ),
        f"- deterministic pass: {det_pass}",
        "",
        "## Uncertainty Robustness (No Refit)",
        f"- n_mc: {n_mc}",
        f"- pass_rate: {pass_rate:.4f} (gate={pass_gate:.2f})",
        "",
        "## Required Next Step",
        f"- {required_next}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2011] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2011] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2011] verdict={verdict}, pass_rate={pass_rate:.4f}")


if __name__ == "__main__":
    main()
