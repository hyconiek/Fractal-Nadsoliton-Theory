#!/usr/bin/env python3
"""
QW-1960: Audit of mass-formula derivation chain.

Goal:
- verify arithmetic consistency in documented derivation,
- quantify circularity risk (Q inferred from masses, then reused),
- check frozen-kernel compatibility of gamma.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean, pstdev
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1960_mass_formula_derivation_audit.json"
OUT_MD = ROOT / "RAPORT_QW1960_MASS_FORMULA_DERIVATION_AUDIT.md"


PARTICLES: List[Tuple[str, float, float]] = [
    ("Top", 0.0, 173_000.0),
    ("Bottom", 7.0, 4_180.0),
    ("Tau", 9.0, 1_776.9),
    ("Charm", 9.0, 1_270.0),
    ("Muon", 14.0, 105.7),
    ("Electron", 24.0, 0.511),
]

M_TOP = 173_000.0
GAMMA_CANON = 1.52


def rel_err_pct(pred: float, exp: float) -> float:
    return abs(pred - exp) / max(abs(exp), 1e-15) * 100.0


def predict_mass(q: float, gamma: float) -> float:
    return M_TOP * (4.0 ** (-(gamma * q / 4.0)))


def infer_q_from_mass(mass_mev: float, gamma: float) -> float:
    return -4.0 * math.log(max(mass_mev / M_TOP, 1e-300), 4.0) / gamma


def gamma_from_q_mass(q: float, mass_mev: float) -> float:
    return -4.0 * math.log(max(mass_mev / M_TOP, 1e-300), 4.0) / max(q, 1e-15)


def load_kernel_gamma_1to4() -> Dict[str, float]:
    p = ROOT / "report_qw1932_physical_reparameterization_eta_scan.json"
    if not p.exists():
        return {"available": False, "gamma_kernel_1to4": float("nan")}
    d = json.loads(p.read_text(encoding="utf-8"))
    sel = d["selected"]
    omega = float(sel["fit"]["omega"])
    phi = float(sel["fit"]["phi"])
    beta = float(sel["fit"]["beta"])
    eta = float(sel["eta"])

    def kernel_fn(dist: np.ndarray) -> np.ndarray:
        return np.cos(omega * dist + phi) / (1.0 + beta * (dist**eta))

    k1 = abs(float(kernel_fn(np.array([1.0]))[0]))
    k4 = abs(float(kernel_fn(np.array([4.0]))[0]))
    gamma_k = -4.0 * math.log(max(k4 / max(k1, 1e-15), 1e-15), 4.0) / 3.0
    return {
        "available": True,
        "omega": omega,
        "phi": phi,
        "beta": beta,
        "eta": eta,
        "abs_k1": k1,
        "abs_k4": k4,
        "gamma_kernel_1to4": float(gamma_k),
    }


def main() -> None:
    rows_canon = []
    canon_errs = []
    gamma_implied = []
    q_inferred = []
    q_rounded_reconstruct_errs = []
    q_delta_abs_non_top = []

    for name, q_model, m_exp in PARTICLES:
        pred = predict_mass(q_model, GAMMA_CANON)
        err = rel_err_pct(pred, m_exp)
        canon_errs.append(err)
        rows_canon.append(
            {
                "particle": name,
                "q_model": q_model,
                "exp_mev": m_exp,
                "pred_mev_gamma_1p52": pred,
                "rel_err_pct_gamma_1p52": err,
            }
        )
        if q_model > 0:
            gamma_implied.append(gamma_from_q_mass(q_model, m_exp))

        q_calc = infer_q_from_mass(m_exp, GAMMA_CANON)
        q_round = round(q_calc)
        if q_model > 0:
            q_delta_abs_non_top.append(abs(q_calc - q_model))
        q_inferred.append(
            {
                "particle": name,
                "q_inferred_from_mass_gamma_1p52": q_calc,
                "q_rounded": q_round,
                "delta_q": q_calc - q_round,
            }
        )
        pred_from_q_round = predict_mass(float(q_round), GAMMA_CANON)
        q_rounded_reconstruct_errs.append(rel_err_pct(pred_from_q_round, m_exp))

    gamma_mean = mean(gamma_implied)
    gamma_std = pstdev(gamma_implied)

    gamma_text_numerator = 2.26
    gamma_text_denominator = 1.77
    gamma_text_value = gamma_text_numerator / gamma_text_denominator
    gamma_text_delta_vs_152 = gamma_text_value - GAMMA_CANON

    rows_text_gamma = []
    text_gamma_errs = []
    for name, q_model, m_exp in PARTICLES:
        pred_text = predict_mass(q_model, gamma_text_value)
        err_text = rel_err_pct(pred_text, m_exp)
        rows_text_gamma.append(
            {
                "particle": name,
                "q_model": q_model,
                "exp_mev": m_exp,
                "pred_mev_gamma_text_ratio": pred_text,
                "rel_err_pct_gamma_text_ratio": err_text,
            }
        )
        text_gamma_errs.append(err_text)

    tau = next(x for x in PARTICLES if x[0] == "Tau")
    charm = next(x for x in PARTICLES if x[0] == "Charm")
    same_q_pred = predict_mass(9.0, GAMMA_CANON)
    tau_charm_mass_ratio_exp = tau[2] / charm[2]
    tau_charm_degeneracy_gap_pct = rel_err_pct(tau[2], charm[2])

    kernel_info = load_kernel_gamma_1to4()
    kernel_gamma = kernel_info.get("gamma_kernel_1to4", float("nan"))
    kernel_delta = float(kernel_gamma - GAMMA_CANON) if kernel_info["available"] else float("nan")
    mean_abs_q_delta_vs_model = float(mean(q_delta_abs_non_top)) if q_delta_abs_non_top else float("nan")

    flags = {
        "arithmetic_inconsistency_gamma_origin": bool(abs(gamma_text_delta_vs_152) > 0.10),
        "q_mapping_circularity_risk": bool(mean_abs_q_delta_vs_model < 0.40),
        "tau_charm_same_q_non_degeneracy": bool(tau_charm_degeneracy_gap_pct > 20.0),
        "frozen_kernel_gamma_incompatibility": bool(kernel_info["available"] and abs(kernel_delta) > 0.30),
    }

    if sum(bool(v) for v in flags.values()) >= 3:
        verdict = "DERIVATION_CONTAINS_MATERIAL_ERRORS_AND_CIRCULAR_STEPS"
    elif any(flags.values()):
        verdict = "DERIVATION_PARTIALLY_INCONSISTENT_NEEDS_REPAIR"
    else:
        verdict = "DERIVATION_AUDIT_NO_CRITICAL_ISSUES"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "formula": "M(Q)=M_top*4^(-(gamma*Q/4))",
            "canonical_gamma": GAMMA_CANON,
            "m_top_mev": M_TOP,
            "particle_table": [{"name": n, "q_model": q, "mass_mev": m} for n, q, m in PARTICLES],
        },
        "canonical_gamma_check": {
            "rows": rows_canon,
            "mean_rel_err_pct": float(mean(canon_errs)),
            "max_rel_err_pct": float(max(canon_errs)),
            "gamma_implied_from_each_particle_q_gt_0": gamma_implied,
            "gamma_implied_mean": float(gamma_mean),
            "gamma_implied_std": float(gamma_std),
        },
        "gamma_origin_arithmetic_check": {
            "text_claim_expression": "gamma = 2.26/1.77",
            "computed_value": float(gamma_text_value),
            "delta_vs_1p52": float(gamma_text_delta_vs_152),
            "rows_using_gamma_text_ratio": rows_text_gamma,
            "mean_rel_err_pct": float(mean(text_gamma_errs)),
            "max_rel_err_pct": float(max(text_gamma_errs)),
        },
        "circularity_check": {
            "method": "infer Q from experimental masses using gamma=1.52, round Q, reconstruct masses",
            "q_inferred_rows": q_inferred,
            "mean_abs_delta_q_inferred_vs_q_model_non_top": mean_abs_q_delta_vs_model,
            "mean_reconstruct_err_pct_after_rounding": float(mean(q_rounded_reconstruct_errs)),
            "max_reconstruct_err_pct_after_rounding": float(max(q_rounded_reconstruct_errs)),
        },
        "degeneracy_check": {
            "tau_mass_mev": tau[2],
            "charm_mass_mev": charm[2],
            "shared_q_model": 9.0,
            "formula_prediction_for_q9_mev": float(same_q_pred),
            "tau_to_charm_ratio_exp": float(tau_charm_mass_ratio_exp),
            "tau_charm_gap_pct": float(tau_charm_degeneracy_gap_pct),
        },
        "kernel_compatibility_check": {
            **kernel_info,
            "delta_vs_canonical_gamma_1p52": kernel_delta,
        },
        "flags": flags,
        "verdict": verdict,
        "required_next_step": (
            "RE-DERIVE_GAMMA_AND_Q_FROM_INDEPENDENT_TOPOLOGY_WITHOUT_MASS_INVERSION"
            if verdict != "DERIVATION_AUDIT_NO_CRITICAL_ISSUES"
            else "PROCEED_TO_EXTERNAL_CONFIRMATORY_TESTS"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1960: MASS FORMULA DERIVATION AUDIT",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Canonical Formula Check (gamma=1.52)",
        f"- mean rel error: {out['canonical_gamma_check']['mean_rel_err_pct']:.3f}%",
        f"- max rel error: {out['canonical_gamma_check']['max_rel_err_pct']:.3f}%",
        (
            "- implied gamma mean/std from fixed Q table: "
            f"{out['canonical_gamma_check']['gamma_implied_mean']:.4f} / "
            f"{out['canonical_gamma_check']['gamma_implied_std']:.4f}"
        ),
        "",
        "## Arithmetic Check of Documented Gamma Origin",
        "- expression tested: gamma = 2.26 / 1.77",
        f"- computed gamma: {gamma_text_value:.6f}",
        f"- delta vs 1.52: {gamma_text_delta_vs_152:+.6f}",
        (
            "- if this gamma is used in mass law: "
            f"mean/max rel error = {mean(text_gamma_errs):.3f}% / {max(text_gamma_errs):.3f}%"
        ),
        "",
        "## Circularity Check",
        f"- mean |Q_inferred - Q_model| (non-top): {mean_abs_q_delta_vs_model:.3f}",
        (
            "- infer Q from masses (gamma=1.52) -> round Q -> reconstruct masses: "
            f"mean/max rel error = {mean(q_rounded_reconstruct_errs):.3f}% / "
            f"{max(q_rounded_reconstruct_errs):.3f}%"
        ),
        "- interpretation: low error here is not independent confirmation.",
        "",
        "## Structural Degeneracy Check",
        "- tau and charm share Q=9 in current mapping.",
        f"- experimental tau/charm ratio: {tau_charm_mass_ratio_exp:.4f}",
        f"- tau/charm gap: {tau_charm_degeneracy_gap_pct:.3f}%",
        f"- formula prediction at Q=9: {same_q_pred:.3f} MeV (single value for both)",
        "",
        "## Frozen-Kernel Compatibility",
    ]
    if kernel_info["available"]:
        lines.extend(
            [
                f"- gamma from kernel (1->4): {kernel_gamma:.6f}",
                f"- delta vs 1.52: {kernel_delta:+.6f}",
            ]
        )
    else:
        lines.append("- kernel report unavailable, check skipped.")

    lines.extend(
        [
            "",
            "## Flags",
            f"- arithmetic_inconsistency_gamma_origin: {flags['arithmetic_inconsistency_gamma_origin']}",
            f"- q_mapping_circularity_risk: {flags['q_mapping_circularity_risk']}",
            f"- tau_charm_same_q_non_degeneracy: {flags['tau_charm_same_q_non_degeneracy']}",
            f"- frozen_kernel_gamma_incompatibility: {flags['frozen_kernel_gamma_incompatibility']}",
            "",
            "## Required Next Step",
            f"- {out['required_next_step']}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1960] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1960] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1960] verdict={verdict}")


if __name__ == "__main__":
    main()
