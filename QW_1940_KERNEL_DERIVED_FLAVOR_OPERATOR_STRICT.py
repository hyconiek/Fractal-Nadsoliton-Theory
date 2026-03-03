#!/usr/bin/env python3
"""
QW-1940: Strict kernel-derived flavor operator (no CKM/PMNS separate tuning).

Rules:
- one frozen kernel from QW-1932,
- one deterministic mapping from kernel moments -> operator parameters,
- same mapping used for CKM and PMNS,
- no per-sector calibration loops.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict

import numpy as np
from scipy.linalg import polar


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1940_kernel_derived_flavor_operator_strict.json"
OUT_MD = ROOT / "RAPORT_QW1940_KERNEL_DERIVED_FLAVOR_OPERATOR_STRICT.md"


CKM_REF = np.array(
    [
        [0.97401, 0.22650, 0.00361],
        [0.22636, 0.97320, 0.04053],
        [0.00854, 0.03978, 0.99917],
    ],
    dtype=float,
)
PMNS_REF = np.array(
    [
        [0.821, 0.550, 0.150],
        [0.432, 0.582, 0.693],
        [0.378, 0.598, 0.707],
    ],
    dtype=float,
)


def kernel_fn(d: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d + phi) / (1.0 + beta * (d**eta))


def cyclic_distance_matrix(q_left: np.ndarray, q_right: np.ndarray, modulus: int = 24) -> np.ndarray:
    dq = np.abs(q_left[:, None] - q_right[None, :]) % float(modulus)
    return np.minimum(dq, float(modulus) - dq)


def derive_shared_params(omega: float, phi: float, beta: float, eta: float) -> Dict[str, float]:
    d = np.arange(1.0, 13.0, dtype=float)
    k = np.abs(kernel_fn(d, omega, phi, beta, eta))
    w = k / max(np.sum(k), 1e-15)

    mean_d = float(np.sum(w * d))
    var_d = float(np.sum(w * (d - mean_d) ** 2))
    decay_ratio = float(k[3] / max(k[0], 1e-15))

    # Deterministic mapping only (no fit): one shared parameter vector.
    p_amp = float(np.clip(1.0 + 0.60 * np.tanh((mean_d - 2.0) / 2.0), 0.30, 2.20))
    r_dist = float(np.clip(np.tanh((var_d - 1.0) / 2.5), -1.20, 1.80))
    phase_scale = float(np.clip(1.0 + 0.70 * np.tanh(abs(phi)) + 0.25 * np.tanh(1.0 - decay_ratio), 0.00, 3.00))

    return {
        "p_amp": p_amp,
        "r_dist": r_dist,
        "phase_scale": phase_scale,
        "moments": {
            "mean_d": mean_d,
            "var_d": var_d,
            "decay_ratio_k4_over_k1": decay_ratio,
        },
    }


def flavor_prediction_abs(
    q_left: np.ndarray,
    q_right: np.ndarray,
    p_amp: float,
    r_dist: float,
    phase_scale: float,
    omega: float,
    phi: float,
    beta: float,
    eta: float,
) -> np.ndarray:
    n = len(q_left)
    d = 1.0 + cyclic_distance_matrix(q_left, q_right, modulus=24)
    k = kernel_fn(d, omega=omega, phi=phi, beta=beta, eta=eta)

    amp = np.sign(k) * ((np.abs(k) ** p_amp) * (d**r_dist))
    idx = np.arange(n, dtype=float)
    gen_gap = np.abs(idx[:, None] - idx[None, :])
    sign = np.where((idx[:, None] - idx[None, :]) < 0.0, 1.0, -1.0)
    phase = np.exp(1j * (phi + phase_scale * omega * sign * gen_gap))

    m = amp * phase
    u = polar(m)[0]
    return np.abs(u)


def matrix_metrics(pred: np.ndarray, ref: np.ndarray) -> Dict[str, float]:
    rel = np.abs(pred - ref) / np.clip(ref, 1e-12, None)
    return {
        "mean_rel_pct": float(100.0 * np.mean(rel)),
        "max_rel_pct": float(100.0 * np.max(rel)),
    }


def main() -> None:
    d1932 = json.loads((ROOT / "report_qw1932_physical_reparameterization_eta_scan.json").read_text(encoding="utf-8"))
    sel = d1932["selected"]
    omega = float(sel["fit"]["omega"])
    phi = float(sel["fit"]["phi"])
    beta = float(sel["fit"]["beta"])
    eta = float(sel["eta"])

    shared = derive_shared_params(omega, phi, beta, eta)
    p_amp = float(shared["p_amp"])
    r_dist = float(shared["r_dist"])
    phase_scale = float(shared["phase_scale"])

    q_up = np.array([0.0, 9.0, 14.0], dtype=float)
    q_down = np.array([7.0, 9.0, 14.0], dtype=float)
    q_nu = np.array([0.0, 1.0, 2.0], dtype=float)
    q_lep = np.array([24.0, 14.0, 9.0], dtype=float)

    ckm_pred = flavor_prediction_abs(q_up, q_down, p_amp, r_dist, phase_scale, omega, phi, beta, eta)
    pmns_pred = flavor_prediction_abs(q_nu, q_lep, p_amp, r_dist, phase_scale, omega, phi, beta, eta)

    m_ckm = matrix_metrics(ckm_pred, CKM_REF)
    m_pmns = matrix_metrics(pmns_pred, PMNS_REF)

    thresholds = {
        "ckm_mean_rel_pct_max": 15.0,
        "pmns_mean_rel_pct_max": 15.0,
    }
    flags = {
        "ckm_mean_rel_pct_le_max": bool(m_ckm["mean_rel_pct"] <= thresholds["ckm_mean_rel_pct_max"]),
        "pmns_mean_rel_pct_le_max": bool(m_pmns["mean_rel_pct"] <= thresholds["pmns_mean_rel_pct_max"]),
    }
    all_pass = bool(all(flags.values()))

    verdict = "KERNEL_DERIVED_FLAVOR_OPERATOR_PASS" if all_pass else "KERNEL_DERIVED_FLAVOR_OPERATOR_FAIL"
    required_next = (
        "FREEZE_FLAVOR_OPERATOR_FOR_TRIPLE_SECTOR_GATE"
        if all_pass
        else "REDESIGN_FLAVOR_MECHANISM_WITHOUT_PER_SECTOR_TUNING"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw1932_physical_reparameterization_eta_scan.json:selected",
        "kernel": {"omega": omega, "phi": phi, "beta": beta, "eta": eta},
        "shared_params_from_kernel": shared,
        "predictions": {
            "ckm_pred_abs": ckm_pred.tolist(),
            "pmns_pred_abs": pmns_pred.tolist(),
        },
        "metrics": {
            "ckm": m_ckm,
            "pmns": m_pmns,
        },
        "thresholds": thresholds,
        "flags": flags,
        "all_pass": all_pass,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1940: KERNEL-DERIVED FLAVOR OPERATOR (STRICT)",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Kernel: omega={omega:.6f}, phi={phi:.6f}, beta={beta:.6f}, eta={eta:.2f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Shared Params (Derived, No Fit)",
        f"- p_amp: {p_amp:.6f}",
        f"- r_dist: {r_dist:.6f}",
        f"- phase_scale: {phase_scale:.6f}",
        "",
        "## Flavor Errors",
        f"- CKM mean/max rel%: {m_ckm['mean_rel_pct']:.3f}/{m_ckm['max_rel_pct']:.3f}",
        f"- PMNS mean/max rel%: {m_pmns['mean_rel_pct']:.3f}/{m_pmns['max_rel_pct']:.3f}",
        "",
        "## Flags",
        f"- ckm_mean_rel_pct_le_max: {flags['ckm_mean_rel_pct_le_max']}",
        f"- pmns_mean_rel_pct_le_max: {flags['pmns_mean_rel_pct_le_max']}",
        f"- all_pass: {all_pass}",
        "",
        "## Required Next Step",
        f"- {required_next}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1940] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1940] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1940] verdict={verdict}")


if __name__ == "__main__":
    main()

