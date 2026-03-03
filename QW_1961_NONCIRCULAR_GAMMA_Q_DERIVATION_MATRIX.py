#!/usr/bin/env python3
"""
QW-1961: Non-circular gamma/Q derivation matrix for mass sector.

Protocol:
- no inversion Q <- mass,
- no parameter fitting/optimization,
- only documented formulas + frozen artifacts.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1961_noncircular_gamma_q_derivation_matrix.json"
OUT_MD = ROOT / "RAPORT_QW1961_NONCIRCULAR_GAMMA_Q_DERIVATION_MATRIX.md"


M_TOP = 173_000.0
N_GRAV = 2.26
D_F = 4.0 * math.log(2.0)
TAU_CHARM_RATIO_EXP = 1776.9 / 1270.0


PARTICLES = [
    ("Top", 173_000.0, 0),
    ("Bottom", 4_180.0, +1),
    ("Tau", 1_776.9, -1),
    ("Charm", 1_270.0, +1),
    ("Muon", 105.7, -1),
    ("Electron", 0.511, -1),
]


Q_ASSIGNMENTS: Dict[str, Dict[str, float]] = {
    # Historical ladder from QW-1159/QW-1160 style mapping.
    "legacy_fibonacci": {
        "Top": 0.0,
        "Bottom": 7.0,
        "Tau": 9.0,
        "Charm": 9.0,
        "Muon": 14.0,
        "Electron": 24.0,
    },
    # Alternative non-degenerate charm node suggested in historical notes.
    "charm_split_fibonacci": {
        "Top": 0.0,
        "Bottom": 7.0,
        "Tau": 9.0,
        "Charm": 10.0,
        "Muon": 14.0,
        "Electron": 24.0,
    },
}


def rel_err_pct(pred: float, exp: float) -> float:
    return abs(pred - exp) / max(abs(exp), 1e-15) * 100.0


def kernel_gamma_1to4() -> Dict[str, float]:
    d1932 = json.loads((ROOT / "report_qw1932_physical_reparameterization_eta_scan.json").read_text(encoding="utf-8"))
    sel = d1932["selected"]
    omega = float(sel["fit"]["omega"])
    phi = float(sel["fit"]["phi"])
    beta = float(sel["fit"]["beta"])
    eta = float(sel["eta"])

    def k(dist: np.ndarray) -> np.ndarray:
        return np.cos(omega * dist + phi) / (1.0 + beta * (dist**eta))

    k1 = abs(float(k(np.array([1.0]))[0]))
    k4 = abs(float(k(np.array([4.0]))[0]))
    gamma = -4.0 * math.log(max(k4 / max(k1, 1e-15), 1e-15), 4.0) / 3.0
    return {
        "omega": omega,
        "phi": phi,
        "beta": beta,
        "eta": eta,
        "abs_k1": k1,
        "abs_k4": k4,
        "gamma_kernel_1to4": float(gamma),
    }


def info_split_delta() -> Dict[str, float]:
    d1958 = json.loads((ROOT / "report_qw1958_formal_nogo_and_minimal_lagrangian_term.json").read_text(encoding="utf-8"))
    coeffs = d1958["minimal_eft_repair_term"]["deterministic_coefficients"]
    lam = float(coeffs["lambda_I"])
    mu = float(coeffs["mu_I"])
    delta = (lam - mu) / max(lam + mu, 1e-15)
    return {"lambda_I": lam, "mu_I": mu, "delta_info": float(delta)}


def predict_mass(q_eff: float, gamma: float) -> float:
    return M_TOP * (4.0 ** (-(gamma * q_eff / 4.0)))


def evaluate_variant(
    gamma: float,
    q_map: Dict[str, float],
    use_info_split: bool,
    delta_info: float,
) -> Dict[str, object]:
    rows: List[Dict[str, float]] = []
    errs = []
    pred_by_name: Dict[str, float] = {}

    for name, exp_mev, sector_sign in PARTICLES:
        q_base = float(q_map[name])
        q_eff = q_base + (sector_sign * delta_info if use_info_split else 0.0)
        pred = predict_mass(q_eff, gamma)
        err = rel_err_pct(pred, exp_mev)
        rows.append(
            {
                "particle": name,
                "exp_mev": exp_mev,
                "q_base": q_base,
                "q_eff": float(q_eff),
                "pred_mev": float(pred),
                "rel_err_pct": float(err),
            }
        )
        errs.append(err)
        pred_by_name[name] = pred

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
    kg = kernel_gamma_1to4()
    info = info_split_delta()

    gamma_candidates = {
        # Reference (historical frozen constant, not a derivation baseline).
        "canonical_frozen_1p52_reference": 1.52,
        # Documented formula in tex around line ~2812.
        "derived_force_energy_2n_over_3": (2.0 * N_GRAV) / 3.0,
        # Documented but arithmetically inconsistent path around line ~3322.
        "derived_ratio_n_over_df_minus_1": N_GRAV / (D_F - 1.0),
        # Strictly from frozen kernel report.
        "derived_kernel_d1_to_d4": float(kg["gamma_kernel_1to4"]),
    }

    variants: List[Dict[str, object]] = []
    for q_name, q_map in Q_ASSIGNMENTS.items():
        for g_name, gamma in gamma_candidates.items():
            for split_name, use_split in [
                ("no_info_split", False),
                ("info_split_from_qw1958", True),
            ]:
                eval_out = evaluate_variant(
                    gamma=float(gamma),
                    q_map=q_map,
                    use_info_split=use_split,
                    delta_info=float(info["delta_info"]),
                )
                variants.append(
                    {
                        "q_assignment": q_name,
                        "gamma_source": g_name,
                        "gamma_value": float(gamma),
                        "split_mode": split_name,
                        **eval_out,
                    }
                )

    variants_sorted = sorted(variants, key=lambda x: x["score"])
    best_overall = variants_sorted[0]

    # Strictly non-circular derivation branch: exclude canonical frozen reference.
    noncircular_only = [v for v in variants if v["gamma_source"] != "canonical_frozen_1p52_reference"]
    noncircular_sorted = sorted(noncircular_only, key=lambda x: x["score"])
    best_noncircular = noncircular_sorted[0]

    noncircular_pass_count = sum(1 for v in noncircular_only if v["all_pass"])
    overall_pass_count = sum(1 for v in variants if v["all_pass"])

    verdict = (
        "NONCIRCULAR_DERIVATION_HAS_STRICT_PASS_CANDIDATE"
        if noncircular_pass_count > 0
        else "NONCIRCULAR_DERIVATION_NO_STRICT_PASS_CANDIDATE"
    )
    required_next = (
        "PROMOTE_BEST_NONCIRCULAR_BRANCH_TO_UNIFIED_TRIAD_TEST"
        if noncircular_pass_count > 0
        else "REWORK_MASS_LINK_OR_KERNEL_DERIVATION_BEFORE_TRIAD_LOCK"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "m_top_mev": M_TOP,
            "n_grav": N_GRAV,
            "d_f": D_F,
            "tau_charm_ratio_exp": TAU_CHARM_RATIO_EXP,
            "q_assignments": Q_ASSIGNMENTS,
            "gamma_candidates": gamma_candidates,
        },
        "kernel_gamma_source": kg,
        "info_split_source": info,
        "thresholds": {
            "mean_rel_err_pct_max": 15.0,
            "max_rel_err_pct_max": 35.0,
            "tau_charm_ratio_rel_err_pct_max": 20.0,
        },
        "summary": {
            "n_variants_total": len(variants),
            "overall_pass_count": int(overall_pass_count),
            "noncircular_pass_count": int(noncircular_pass_count),
            "best_overall": best_overall,
            "best_noncircular": best_noncircular,
            "top5": variants_sorted[:5],
        },
        "variants": variants_sorted,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    bo = best_overall
    bn = best_noncircular
    lines = [
        "# RAPORT QW-1961: NONCIRCULAR GAMMA/Q DERIVATION MATRIX",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- variants total: {len(variants)}",
        f"- noncircular pass count: {noncircular_pass_count}",
        "",
        "## Key Inputs",
        f"- D_f = 4 ln 2 = {D_F:.6f}",
        f"- n_grav = {N_GRAV:.4f}",
        f"- gamma_kernel(d1->d4) = {kg['gamma_kernel_1to4']:.6f}",
        f"- delta_info (from QW-1958 lambda/mu) = {info['delta_info']:.6f}",
        "",
        "## Best Overall Variant",
        (
            f"- q/gamma/split: {bo['q_assignment']} / {bo['gamma_source']}="
            f"{bo['gamma_value']:.6f} / {bo['split_mode']}"
        ),
        (
            f"- mean/max err: {bo['metrics']['mean_rel_err_pct']:.3f}% / "
            f"{bo['metrics']['max_rel_err_pct']:.3f}%"
        ),
        (
            f"- tau/charm ratio pred/exp/error: "
            f"{bo['metrics']['tau_charm_ratio_pred']:.4f} / "
            f"{bo['metrics']['tau_charm_ratio_exp']:.4f} / "
            f"{bo['metrics']['tau_charm_ratio_rel_err_pct']:.3f}%"
        ),
        f"- all_pass: {bo['all_pass']}",
        "",
        "## Best Noncircular Variant",
        (
            f"- q/gamma/split: {bn['q_assignment']} / {bn['gamma_source']}="
            f"{bn['gamma_value']:.6f} / {bn['split_mode']}"
        ),
        (
            f"- mean/max err: {bn['metrics']['mean_rel_err_pct']:.3f}% / "
            f"{bn['metrics']['max_rel_err_pct']:.3f}%"
        ),
        (
            f"- tau/charm ratio pred/exp/error: "
            f"{bn['metrics']['tau_charm_ratio_pred']:.4f} / "
            f"{bn['metrics']['tau_charm_ratio_exp']:.4f} / "
            f"{bn['metrics']['tau_charm_ratio_rel_err_pct']:.3f}%"
        ),
        f"- all_pass: {bn['all_pass']}",
        "",
        "## Required Next Step",
        f"- {required_next}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1961] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1961] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1961] verdict={verdict} noncircular_pass={noncircular_pass_count}")


if __name__ == "__main__":
    main()
