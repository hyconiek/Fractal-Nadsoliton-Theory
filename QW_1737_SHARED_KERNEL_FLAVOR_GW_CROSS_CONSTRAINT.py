#!/usr/bin/env python3
"""
QW-1737: Shared kernel flavor+GW cross-constraint and identifiability.

This script is conservative:
1) Uses empirical gates from prior strict reports (1720/1723 flavor, 1725/1726 GW).
2) Combines with parameter distributions inferred in 1734/1735/1736.
3) Reports whether any scientifically credible shared region exists.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1737_shared_kernel_flavor_gw_cross_constraint.json"
OUT_MD = ROOT / "RAPORT_QW1737_SHARED_KERNEL_FLAVOR_GW_CROSS_CONSTRAINT.md"


def load(path: str) -> Dict:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def ckm_matrix_abs(theta12: float, theta23: float, theta13: float, delta: float = math.pi / 2.0) -> np.ndarray:
    s12, c12 = math.sin(theta12), math.cos(theta12)
    s23, c23 = math.sin(theta23), math.cos(theta23)
    s13, c13 = math.sin(theta13), math.cos(theta13)
    e = complex(math.cos(delta), -math.sin(delta))

    # Standard PDG parameterization.
    v = np.array(
        [
            [c12 * c13, s12 * c13, s13 * e],
            [
                -s12 * c23 - c12 * s23 * s13 * np.conjugate(e),
                c12 * c23 - s12 * s23 * s13 * np.conjugate(e),
                s23 * c13,
            ],
            [
                s12 * s23 - c12 * c23 * s13 * np.conjugate(e),
                -c12 * s23 - s12 * c23 * s13 * np.conjugate(e),
                c23 * c13,
            ],
        ],
        dtype=complex,
    )
    return np.abs(v)


def kernel_values(alpha: float, omega: float, phi: float, beta: float, dmax: int = 12) -> np.ndarray:
    d = np.arange(1, dmax + 1, dtype=float)
    return alpha * np.cos(omega * d + phi) / (1.0 + beta * d)


def flavor_proxy_score(omega: float, phi: float, beta: float, alpha: float = 4.0 * math.log(2.0)) -> Dict[str, float]:
    k = np.abs(kernel_values(alpha, omega, phi, beta, dmax=6))
    # Conservative, fixed mapping from kernel magnitudes to mixing angles.
    th12 = math.atan(min(1.0, k[0]))
    th23 = math.atan(min(1.0, k[1]))
    th13 = 0.35 * math.atan(min(1.0, k[2]))

    # PMNS proxy: larger angles from farther-octave couplings.
    tn12 = math.atan(min(1.0, k[2] + 0.15))
    tn23 = math.atan(min(1.0, k[3] + 0.15))
    tn13 = 0.7 * math.atan(min(1.0, k[4] + 0.05))

    ckm_pred = ckm_matrix_abs(th12, th23, th13, delta=math.pi / 2.0)
    pmns_pred = ckm_matrix_abs(tn12, tn23, tn13, delta=math.pi / 2.0)

    ckm_ref = np.array(
        [
            [0.97401, 0.22650, 0.00361],
            [0.22636, 0.97320, 0.04053],
            [0.00854, 0.04040, 0.99914],
        ],
        dtype=float,
    )
    pmns_ref = np.array(
        [
            [0.821, 0.550, 0.150],
            [0.432, 0.577, 0.695],
            [0.376, 0.604, 0.702],
        ],
        dtype=float,
    )

    ckm_rel = np.abs(ckm_pred - ckm_ref) / np.clip(ckm_ref, 1e-8, None)
    pmns_rel = np.abs(pmns_pred - pmns_ref) / np.clip(pmns_ref, 1e-8, None)
    ckm_mean_pct = float(100.0 * np.mean(ckm_rel))
    pmns_mean_pct = float(100.0 * np.mean(pmns_rel))
    return {
        "ckm_mean_rel_pct": ckm_mean_pct,
        "pmns_mean_rel_pct": pmns_mean_pct,
        "pass_proxy": (ckm_mean_pct <= 15.0 and pmns_mean_pct <= 15.0),
    }


def gw_proxy_score(omega: float, phi: float, beta: float, alpha: float = 4.0 * math.log(2.0)) -> Dict[str, float]:
    k = np.abs(kernel_values(alpha, omega, phi, beta, dmax=12))
    near = float(np.mean(k[:3]))
    far = float(np.mean(k[6:12]))
    ratio = far / max(near, 1e-12)

    # Conservative null-consistent criterion from 1725/1726:
    # no robust detector-spanning anomaly => avoid extreme far/near amplification.
    # sensitivity band chosen to test identifiability robustness.
    pass_loose = ratio <= 2.0
    pass_strict = ratio <= 1.5
    return {"far_to_near_ratio": ratio, "pass_loose": pass_loose, "pass_strict": pass_strict}


def sample_joint_parameters(d1734: Dict, d1735: Dict, n_samples: int = 3000, seed: int = 1737) -> np.ndarray:
    rng = np.random.default_rng(seed)
    b_rows = d1734.get("rows_head_120", [])
    o_rows = d1735.get("rows", [])
    if not b_rows or not o_rows:
        raise RuntimeError("Missing rows in 1734/1735 reports.")

    beta_pool = np.array([float(r["beta_hat"]) for r in b_rows], dtype=float)
    omega_pool = np.array([float(r["omega_hat"]) for r in o_rows], dtype=float)
    phi_pool = np.array([float(r["phi_hat"]) for r in o_rows], dtype=float)

    idx_b = rng.integers(0, len(beta_pool), size=n_samples)
    idx_o = rng.integers(0, len(omega_pool), size=n_samples)
    beta = beta_pool[idx_b]
    omega = omega_pool[idx_o]
    phi = phi_pool[idx_o]
    return np.stack([omega, phi, beta], axis=1)


def entropy_of_weights(w: np.ndarray) -> float:
    w = np.asarray(w, dtype=float)
    w = np.clip(w, 1e-15, None)
    w = w / np.sum(w)
    return float(-np.sum(w * np.log(w)))


def main() -> None:
    d1720 = load("report_qw1720_flavor_extended_operator.json")
    d1723 = load("report_qw1723_unified_effective_hamiltonian_integration.json")
    d1725 = load("report_qw1725_gw_strict_cross_hurst_reanalysis.json")
    d1726 = load("report_qw1726_gw_fin_projection_retest.json")
    d1734 = load("report_qw1734_micro_beta_tors_derivation.json")
    d1735 = load("report_qw1735_omega_phi_from_lattice_dynamics.json")
    d1736 = load("report_qw1736_kernel_node_bayesian_model_selection.json")

    # Empirical hard gates from previous strict studies.
    flavor_gate_empirical = (
        d1720["best"]["ckm_error"]["mean_rel_pct"] <= d1720["threshold_mean_rel_pct"]
        and d1720["best"]["pmns_error"]["mean_rel_pct"] <= d1720["threshold_mean_rel_pct"]
        and d1723["best"]["flavor"]["ckm_error"]["mean_rel_pct"] <= d1723["thresholds"]["flavor_ckm_mean_rel_pct_max"]
        and d1723["best"]["flavor"]["pmns_error"]["mean_rel_pct"] <= d1723["thresholds"]["flavor_pmns_mean_rel_pct_max"]
    )
    gw_gate_empirical = (
        d1725["null_phase_randomized"]["p_lower"] <= d1725["thresholds"]["null_p_max"]
        and d1725["stability"]["length_spread"] <= d1725["thresholds"]["length_spread_max"]
        and d1725["lag_profile"]["corr_at_plus_10ms"] >= d1725["thresholds"]["lag_pos_corr_at_10ms_min"]
        and d1726["best_row_closest_to_031"]["shared_background"]["prob_near_031_pm_002"] > 0.05
    )

    samples = sample_joint_parameters(d1734, d1735, n_samples=3000, seed=1737)
    node_model = d1736["best_model"]

    rows: List[Dict[str, float]] = []
    hard_shared_count = 0
    loose_shared_count = 0

    flavor_weights = []
    gw_weights = []
    joint_weights = []

    for omega, phi, beta in samples:
        fs = flavor_proxy_score(float(omega), float(phi), float(beta))
        gs = gw_proxy_score(float(omega), float(phi), float(beta))

        # Node-model consistency soft factor.
        if node_model == "M_C_4over3_plus_4n":
            node_factor = math.exp(-((omega - math.pi / 4.0) / 0.18) ** 2)
        elif node_model == "M_A_2_plus_3n":
            node_factor = math.exp(-((omega - (math.pi / 3.0)) / 0.18) ** 2)
        else:
            node_factor = math.exp(-((omega - (math.pi / 6.0)) / 0.18) ** 2)

        wf = math.exp(-0.015 * fs["ckm_mean_rel_pct"] - 0.010 * fs["pmns_mean_rel_pct"])
        wg = math.exp(-0.9 * max(0.0, gs["far_to_near_ratio"] - 1.5))
        wj = wf * wg * node_factor

        flavor_weights.append(wf)
        gw_weights.append(wg)
        joint_weights.append(wj)

        pass_loose = fs["pass_proxy"] and gs["pass_loose"] and node_factor > 0.35
        pass_strict = fs["pass_proxy"] and gs["pass_strict"] and node_factor > 0.55
        if pass_loose:
            loose_shared_count += 1
        if pass_strict:
            hard_shared_count += 1

        rows.append(
            {
                "omega": float(omega),
                "phi": float(phi),
                "beta": float(beta),
                "flavor_ckm_mean_rel_pct": fs["ckm_mean_rel_pct"],
                "flavor_pmns_mean_rel_pct": fs["pmns_mean_rel_pct"],
                "gw_far_to_near_ratio": gs["far_to_near_ratio"],
                "node_factor": node_factor,
                "pass_loose": pass_loose,
                "pass_strict": pass_strict,
                "joint_weight": wj,
            }
        )

    n = len(rows)
    frac_loose = loose_shared_count / n
    frac_strict = hard_shared_count / n

    wf = np.array(flavor_weights, dtype=float)
    wg = np.array(gw_weights, dtype=float)
    wj = np.array(joint_weights, dtype=float)

    h_uniform = math.log(n)
    h_joint = entropy_of_weights(wj)
    identifiability_reduction = max(0.0, (h_uniform - h_joint) / max(h_uniform, 1e-12))

    # Final readiness is strict: empirical gates + sampled shared region must both exist.
    if flavor_gate_empirical and gw_gate_empirical and frac_strict >= 0.05:
        verdict = "SHARED_KERNEL_REGION_SUPPORTED"
    elif frac_loose > 0.02 and identifiability_reduction > 0.10:
        verdict = "SHARED_KERNEL_REGION_WEAK_PROXY_ONLY"
    else:
        verdict = "SHARED_KERNEL_REGION_NOT_SUPPORTED"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "empirical_hard_gates": {
            "flavor_gate_from_1720_1723": flavor_gate_empirical,
            "gw_gate_from_1725_1726": gw_gate_empirical,
            "note": "Both must pass for strict empirical support.",
        },
        "sampling": {
            "n_samples": n,
            "node_model_from_1736": node_model,
            "shared_fraction_loose": frac_loose,
            "shared_fraction_strict": frac_strict,
        },
        "identifiability": {
            "entropy_uniform": h_uniform,
            "entropy_joint_weighted": h_joint,
            "relative_entropy_reduction": identifiability_reduction,
            "interpretation": "0 means non-identifiable; larger means stronger concentration.",
        },
        "weight_summary": {
            "flavor_weight_mean": float(np.mean(wf)),
            "gw_weight_mean": float(np.mean(wg)),
            "joint_weight_mean": float(np.mean(wj)),
        },
        "verdict": verdict,
        "rows_head_200": rows[:200],
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1737: SHARED KERNEL FLAVOR+GW CROSS-CONSTRAINT",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Empirical flavor gate (1720/1723): {flavor_gate_empirical}",
        f"- Empirical GW gate (1725/1726): {gw_gate_empirical}",
        f"- Shared fraction (loose): {frac_loose:.4f}",
        f"- Shared fraction (strict): {frac_strict:.4f}",
        f"- Identifiability reduction: {identifiability_reduction:.4f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Note",
        "- If empirical hard gates fail, proxy overlap alone is not treated as closure evidence.",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1737] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1737] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
