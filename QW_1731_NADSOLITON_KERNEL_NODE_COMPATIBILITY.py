#!/usr/bin/env python3
"""
QW-1731: Compatibility test between Nadsoliton characteristics and kernel nodes.

Purpose:
1) Check if characteristic priors (omega=pi/4, phi=pi/6) reproduce claimed node sets.
2) Find best-fit (omega, phi) for claimed node sets and compare distance from priors.
3) Quantify whether current node narrative is mathematically compatible.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1731_nadsoliton_kernel_node_compatibility.json"
OUT_MD = ROOT / "RAPORT_QW1731_NADSOLITON_KERNEL_NODE_COMPATIBILITY.md"


def kernel(d: np.ndarray, alpha: float, omega: float, phi: float, beta: float) -> np.ndarray:
    return alpha * np.cos(omega * d + phi) / (1.0 + beta * d)


def node_loss(omega: float, phi: float, nodes: np.ndarray) -> float:
    vals = np.cos(omega * nodes + phi)
    return float(np.mean(vals * vals))


def nearest_distance_to_sequence(values: np.ndarray, seq_start: float, seq_step: float, n_terms: int = 8) -> float:
    seq = np.array([seq_start + i * seq_step for i in range(n_terms)], dtype=float)
    dists = []
    for v in values:
        dists.append(float(np.min(np.abs(seq - v))))
    return float(np.mean(dists))


def grid_optimize_nodes(nodes: np.ndarray) -> Dict[str, float]:
    omega_grid = np.linspace(0.35, 1.25, 1401)
    phi_grid = np.linspace(-math.pi, math.pi, 2001)
    best = {
        "omega": float("nan"),
        "phi": float("nan"),
        "loss": float("inf"),
    }
    for omega in omega_grid:
        vals = omega * nodes
        for phi in phi_grid:
            loss = float(np.mean(np.cos(vals + phi) ** 2))
            if loss < best["loss"]:
                best = {"omega": float(omega), "phi": float(phi), "loss": loss}
    return best


def wrap_angle_diff(a: float, b: float) -> float:
    x = (a - b + math.pi) % (2.0 * math.pi) - math.pi
    return abs(x)


def main() -> None:
    alpha_prior = 4.0 * math.log(2.0)
    omega_prior = math.pi / 4.0
    phi_prior = math.pi / 6.0
    beta_prior = 0.01

    nodes_claim_a = np.array([2.0, 5.0, 8.0, 11.0], dtype=float)
    nodes_claim_b = np.array([2.0, 8.0, 14.0], dtype=float)

    prior_loss_a = node_loss(omega_prior, phi_prior, nodes_claim_a)
    prior_loss_b = node_loss(omega_prior, phi_prior, nodes_claim_b)

    opt_a = grid_optimize_nodes(nodes_claim_a)
    opt_b = grid_optimize_nodes(nodes_claim_b)

    delta_omega_a = abs(opt_a["omega"] - omega_prior) / omega_prior
    delta_phi_a = wrap_angle_diff(opt_a["phi"], phi_prior) / math.pi
    delta_omega_b = abs(opt_b["omega"] - omega_prior) / omega_prior
    delta_phi_b = wrap_angle_diff(opt_b["phi"], phi_prior) / math.pi

    d = np.arange(1.0, 13.0, dtype=float)
    k_prior = kernel(d, alpha_prior, omega_prior, phi_prior, beta_prior)
    hierarchy_ratio = float(abs(k_prior[6]) / max(abs(k_prior[0]), 1e-12))  # d=7 vs d=1

    # Analytic zero sequence from priors:
    # d_n = (pi/2 - phi + n*pi)/omega
    first_zero = (math.pi / 2.0 - phi_prior) / omega_prior
    zero_spacing = math.pi / omega_prior

    dist_a_to_prior_seq = nearest_distance_to_sequence(
        nodes_claim_a, seq_start=first_zero, seq_step=zero_spacing, n_terms=10
    )
    dist_b_to_prior_seq = nearest_distance_to_sequence(
        nodes_claim_b, seq_start=first_zero, seq_step=zero_spacing, n_terms=10
    )

    incompatibility_flags: List[str] = []
    if prior_loss_a > 0.2 and opt_a["loss"] < 1e-3:
        incompatibility_flags.append("NODE_SET_A_INCOMPATIBLE_WITH_CHARACTERISTIC_PRIORS")
    if prior_loss_b > 0.2 and opt_b["loss"] < 1e-3:
        incompatibility_flags.append("NODE_SET_B_INCOMPATIBLE_WITH_CHARACTERISTIC_PRIORS")
    if dist_a_to_prior_seq > 0.5:
        incompatibility_flags.append("NODE_SET_A_FAR_FROM_ANALYTIC_ZERO_SEQUENCE")
    if dist_b_to_prior_seq > 0.5:
        incompatibility_flags.append("NODE_SET_B_FAR_FROM_ANALYTIC_ZERO_SEQUENCE")

    if len(incompatibility_flags) >= 3:
        verdict = "NODE_NARRATIVE_STRONGLY_INCONSISTENT"
    elif incompatibility_flags:
        verdict = "NODE_NARRATIVE_PARTIALLY_INCONSISTENT"
    else:
        verdict = "NODE_NARRATIVE_COMPATIBLE"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "characteristic_priors": {
            "alpha_geo": alpha_prior,
            "omega": omega_prior,
            "phi": phi_prior,
            "beta_tors": beta_prior,
        },
        "analytic_zero_sequence_from_priors": {
            "first_zero": first_zero,
            "zero_spacing": zero_spacing,
            "formula": "d_n = (pi/2 - phi + n*pi)/omega",
        },
        "claim_set_A_2_5_8_11": {
            "prior_loss": prior_loss_a,
            "best_fit": opt_a,
            "relative_delta_from_priors": {
                "omega_rel": delta_omega_a,
                "phi_wrapped_over_pi": delta_phi_a,
            },
            "mean_distance_to_prior_zero_sequence": dist_a_to_prior_seq,
        },
        "claim_set_B_2_8_14": {
            "prior_loss": prior_loss_b,
            "best_fit": opt_b,
            "relative_delta_from_priors": {
                "omega_rel": delta_omega_b,
                "phi_wrapped_over_pi": delta_phi_b,
            },
            "mean_distance_to_prior_zero_sequence": dist_b_to_prior_seq,
        },
        "hierarchy_check_from_priors": {
            "abs_K_d7_over_abs_K_d1": hierarchy_ratio,
            "note": "d=7 and d=1 evaluated with prior kernel parameters.",
        },
        "incompatibility_flags": incompatibility_flags,
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1731: NODE COMPATIBILITY AUDIT",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Priors from Nadsoliton Characteristics",
        f"- alpha_geo = {alpha_prior:.12f}",
        f"- omega = {omega_prior:.12f} (pi/4)",
        f"- phi = {phi_prior:.12f} (pi/6)",
        f"- beta_tors = {beta_prior:.5f}",
        "",
        "## Analytic Zero Sequence from Priors",
        f"- first_zero = {first_zero:.4f}",
        f"- zero_spacing = {zero_spacing:.4f}",
        "",
        "## Claim Set A: nodes [2,5,8,11]",
        f"- prior_loss = {prior_loss_a:.6f}",
        f"- best_fit omega = {opt_a['omega']:.6f}, phi = {opt_a['phi']:.6f}, loss = {opt_a['loss']:.6e}",
        f"- relative shift from priors: domega={delta_omega_a:.3f}, dphi/pi={delta_phi_a:.3f}",
        f"- mean distance to prior zero sequence = {dist_a_to_prior_seq:.4f}",
        "",
        "## Claim Set B: nodes [2,8,14]",
        f"- prior_loss = {prior_loss_b:.6f}",
        f"- best_fit omega = {opt_b['omega']:.6f}, phi = {opt_b['phi']:.6f}, loss = {opt_b['loss']:.6e}",
        f"- relative shift from priors: domega={delta_omega_b:.3f}, dphi/pi={delta_phi_b:.3f}",
        f"- mean distance to prior zero sequence = {dist_b_to_prior_seq:.4f}",
        "",
        "## Hierarchy Check",
        f"- abs(K(7))/abs(K(1)) with priors = {hierarchy_ratio:.3f}",
        "",
        "## Flags",
    ]
    if incompatibility_flags:
        for flag in incompatibility_flags:
            lines.append(f"- {flag}")
    else:
        lines.append("- none")
    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1731] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1731] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
