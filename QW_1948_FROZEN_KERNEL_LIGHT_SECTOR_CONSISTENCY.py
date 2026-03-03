#!/usr/bin/env python3
"""
QW-1948: Frozen-kernel light-sector consistency test.

Goal:
- verify whether the same frozen kernel supports a photon-like branch
  with near-linear low-k dispersion and stable effective propagation speed,
- quantify polarization splitting from signed kernel structure.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1948_frozen_kernel_light_sector_consistency.json"
OUT_MD = ROOT / "RAPORT_QW1948_FROZEN_KERNEL_LIGHT_SECTOR_CONSISTENCY.md"


THRESHOLDS = {
    "lowk_linearity_r2_min": 0.995,
    "lowk_group_cv_max": 0.12,
    "polarization_speed_split_max": 0.02,
}


def kernel_fn(d: np.ndarray, omega: float, phi: float, beta: float, eta: float) -> np.ndarray:
    return np.cos(omega * d + phi) / (1.0 + beta * (d**eta))


def linear_fit_through_origin(x: np.ndarray, y: np.ndarray) -> Dict[str, float]:
    num = float(np.sum(x * y))
    den = float(np.sum(x * x))
    slope = num / max(den, 1e-15)
    y_hat = slope * x

    y_bar = float(np.mean(y))
    ss_res = float(np.sum((y - y_hat) ** 2))
    ss_tot = float(np.sum((y - y_bar) ** 2))
    r2 = float(1.0 - ss_res / max(ss_tot, 1e-15))

    rel_max = float(np.max(np.abs(y - y_hat) / np.clip(np.abs(y), 1e-15, None)))
    return {
        "slope": slope,
        "r2": r2,
        "max_rel_residual": rel_max,
    }


def dispersion_from_couplings(k: np.ndarray, d: np.ndarray, j: np.ndarray) -> np.ndarray:
    out = np.zeros_like(k, dtype=float)
    for i, kv in enumerate(k):
        val = 2.0 * np.sum(j * (1.0 - np.cos(kv * d)))
        out[i] = math.sqrt(max(float(val), 1e-15))
    return out


def channel_metrics(k: np.ndarray, omega_k: np.ndarray, lowk_mask: np.ndarray) -> Dict[str, float]:
    fit = linear_fit_through_origin(k[lowk_mask], omega_k[lowk_mask])
    vg = np.gradient(omega_k, k)
    vg_low = vg[lowk_mask]

    vg_mean = float(np.mean(vg_low))
    vg_std = float(np.std(vg_low))
    vg_cv = float(vg_std / max(abs(vg_mean), 1e-15))

    return {
        "c_eff_slope_lowk": float(fit["slope"]),
        "lowk_r2": float(fit["r2"]),
        "lowk_max_rel_residual": float(fit["max_rel_residual"]),
        "group_v_mean_lowk": vg_mean,
        "group_v_std_lowk": vg_std,
        "group_v_cv_lowk": vg_cv,
        "group_v_min_lowk": float(np.min(vg_low)),
        "group_v_max_lowk": float(np.max(vg_low)),
    }


def build_channels(kernel: Dict[str, float]) -> Dict[str, object]:
    d = np.arange(1.0, 25.0, dtype=float)
    k_signed = kernel_fn(d, kernel["omega"], kernel["phi"], kernel["beta"], kernel["eta"])
    k_abs = np.abs(k_signed)
    sign = np.sign(k_signed)
    sign[sign == 0.0] = 1.0

    j_base = k_abs / max(np.sum(k_abs), 1e-15)

    flips = float(np.mean((sign[1:] != sign[:-1]).astype(float)))
    parity_imb = float(
        abs(np.sum(k_abs[(d.astype(int) % 2) == 0]) - np.sum(k_abs[(d.astype(int) % 2) == 1]))
        / max(np.sum(k_abs), 1e-15)
    )
    eps_pol = float(np.clip(0.08 * flips + 0.05 * parity_imb + 0.03 * abs(math.sin(kernel["phi"])), 0.0, 0.15))

    j_plus = np.clip(j_base * (1.0 + eps_pol * sign), 1e-15, None)
    j_minus = np.clip(j_base * (1.0 - eps_pol * sign), 1e-15, None)
    j_plus /= max(np.sum(j_plus), 1e-15)
    j_minus /= max(np.sum(j_minus), 1e-15)

    return {
        "d_1_to_24": d,
        "k_signed_1_to_24": k_signed,
        "k_abs_1_to_24": k_abs,
        "j_base": j_base,
        "j_plus": j_plus,
        "j_minus": j_minus,
        "eps_pol": eps_pol,
        "flip_rate": flips,
        "parity_imbalance": parity_imb,
    }


def main() -> None:
    d1932 = json.loads((ROOT / "report_qw1932_physical_reparameterization_eta_scan.json").read_text(encoding="utf-8"))
    sel = d1932["selected"]
    kernel = {
        "omega": float(sel["fit"]["omega"]),
        "phi": float(sel["fit"]["phi"]),
        "beta": float(sel["fit"]["beta"]),
        "eta": float(sel["eta"]),
    }

    light = build_channels(kernel)
    d = np.array(light["d_1_to_24"], dtype=float)
    j_plus = np.array(light["j_plus"], dtype=float)
    j_minus = np.array(light["j_minus"], dtype=float)
    j_base = np.array(light["j_base"], dtype=float)

    k = np.linspace(0.01, 1.20, 600, dtype=float)
    lowk_mask = (k >= 0.02) & (k <= 0.25)

    w_base = dispersion_from_couplings(k, d, j_base)
    w_plus = dispersion_from_couplings(k, d, j_plus)
    w_minus = dispersion_from_couplings(k, d, j_minus)

    m_base = channel_metrics(k, w_base, lowk_mask)
    m_plus = channel_metrics(k, w_plus, lowk_mask)
    m_minus = channel_metrics(k, w_minus, lowk_mask)

    c_plus = float(m_plus["c_eff_slope_lowk"])
    c_minus = float(m_minus["c_eff_slope_lowk"])
    c_split = float(abs(c_plus - c_minus) / max((abs(c_plus) + abs(c_minus)) / 2.0, 1e-15))

    lowk_r2_min = float(min(m_plus["lowk_r2"], m_minus["lowk_r2"]))
    lowk_group_cv_max = float(max(m_plus["group_v_cv_lowk"], m_minus["group_v_cv_lowk"]))

    flags = {
        "lowk_linearity_r2_ge_min": bool(lowk_r2_min >= THRESHOLDS["lowk_linearity_r2_min"]),
        "lowk_group_cv_le_max": bool(lowk_group_cv_max <= THRESHOLDS["lowk_group_cv_max"]),
        "polarization_speed_split_le_max": bool(c_split <= THRESHOLDS["polarization_speed_split_max"]),
    }
    light_pass = bool(all(flags.values()))

    mass_status = None
    p1947 = ROOT / "report_qw1947_coupling_completeness_and_mass_operator_frontier.json"
    if p1947.exists():
        d1947 = json.loads(p1947.read_text(encoding="utf-8"))
        mass_status = bool(d1947.get("summary", {}).get("strict_pass_exists", False))

    if light_pass and mass_status is False:
        verdict = "LIGHT_SECTOR_PASS_MASS_SECTOR_STILL_BLOCKED"
        required_next = "KEEP_LIGHT_KERNEL_CONSTRAINTS_AND_REPAIR_MASS_LINK_ONLY"
    elif light_pass:
        verdict = "LIGHT_SECTOR_PASS"
        required_next = "INTEGRATE_LIGHT_CONSTRAINTS_IN_TRIPLE_SECTOR_GATE"
    else:
        verdict = "LIGHT_SECTOR_FAIL"
        required_next = "REWORK_KERNEL_OR_EFFECTIVE_PROPAGATION_OPERATOR_FOR_LIGHT"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_source": "report_qw1932_physical_reparameterization_eta_scan.json:selected",
        "kernel": kernel,
        "thresholds": THRESHOLDS,
        "light_channel_construction": {
            "eps_pol": float(light["eps_pol"]),
            "flip_rate": float(light["flip_rate"]),
            "parity_imbalance": float(light["parity_imbalance"]),
        },
        "metrics": {
            "base": m_base,
            "plus": m_plus,
            "minus": m_minus,
            "lowk_r2_min": lowk_r2_min,
            "lowk_group_cv_max": lowk_group_cv_max,
            "polarization_speed_split": c_split,
        },
        "flags": flags,
        "light_pass": light_pass,
        "mass_strict_pass_from_qw1947": mass_status,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1948: FROZEN KERNEL LIGHT-SECTOR CONSISTENCY",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Core Light Metrics",
        f"- low-k linearity R2 (min over +/-): {lowk_r2_min:.6f}",
        f"- low-k group velocity CV (max over +/-): {lowk_group_cv_max:.6f}",
        f"- polarization speed split: {c_split:.6f}",
        "",
        "## Flags",
        f"- lowk_linearity_r2_ge_min: {flags['lowk_linearity_r2_ge_min']}",
        f"- lowk_group_cv_le_max: {flags['lowk_group_cv_le_max']}",
        f"- polarization_speed_split_le_max: {flags['polarization_speed_split_le_max']}",
        f"- light_pass: {light_pass}",
        "",
        "## Link to QW-1947",
        f"- mass_strict_pass_from_qw1947: {mass_status}",
        "",
        "## Required Next Step",
        f"- {required_next}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1948] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1948] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1948] verdict={verdict}")


if __name__ == "__main__":
    main()

