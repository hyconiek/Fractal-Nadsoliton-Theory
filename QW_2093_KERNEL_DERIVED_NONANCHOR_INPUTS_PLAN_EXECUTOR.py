#!/usr/bin/env python3
"""
QW-2093: Kernel-derived non-anchor inputs plan executor (QW-2085..QW-2087).

Purpose:
- execute frozen plan formulas from PLAN_KERNEL_DERIVED_NONANCHOR_INPUTS_QW2085_2087.md,
- generate deterministic kernel-derived input files for T1 gates,
- produce traceable report with formulas/values/sources.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_OBS_2085_2086 = ROOT / "t1_nonanchor_observables_input_qw2085_2086.json"
OUT_OBS_2087 = ROOT / "t1_nonanchor_alpha_s_input_qw2087.json"
OUT_JSON = ROOT / "report_qw2093_kernel_derived_nonanchor_inputs_plan_executor.json"
OUT_MD = ROOT / "RAPORT_QW2093_KERNEL_DERIVED_NONANCHOR_INPUTS_PLAN_EXECUTOR.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def run_alpha_s_one_loop(mu0: float, alpha0: float, mu_target: float) -> float:
    if mu0 <= 0.0 or mu_target <= 0.0 or alpha0 <= 0.0:
        raise ValueError("Scales and alpha0 must be > 0.")
    if abs(mu_target - mu0) < 1e-15:
        return float(alpha0)

    thresholds = sorted([1.27, 4.18, 173.0])
    forward = mu_target > mu0
    cuts = [x for x in thresholds if mu0 < x < mu_target] if forward else [x for x in thresholds if mu_target < x < mu0]
    cuts = sorted(cuts)
    if not forward:
        cuts = list(reversed(cuts))

    boundaries = [mu0] + cuts + [mu_target]
    inv_alpha = 1.0 / alpha0

    for i in range(len(boundaries) - 1):
        a = boundaries[i]
        b = boundaries[i + 1]
        mu_mid = math.sqrt(a * b)
        nf = active_nf_qcd(mu_mid)
        beta0 = 11.0 - (2.0 / 3.0) * nf
        inv_alpha = inv_alpha + (beta0 / (2.0 * math.pi)) * math.log(b / a)

    if inv_alpha <= 0.0:
        raise ValueError("Unphysical QCD running: inverse alpha_s <= 0.")
    return float(1.0 / inv_alpha)


def active_nf_qcd(mu_gev: float) -> int:
    quark_masses = [0.00216, 0.00467, 0.093, 1.27, 4.18, 173.0]
    nf = sum(1 for m in quark_masses if mu_gev >= m)
    return int(max(3, min(6, nf)))


def find_mass_pred_gev(rows: List[Dict], particle: str) -> float:
    for r in rows:
        if str(r.get("particle")) == particle:
            return float(r["pred_mev"]) / 1000.0
    raise KeyError(f"Missing particle in QW-2063 mass rows: {particle}")


def main() -> None:
    r2063 = load_json("report_qw2063_derivational_reconstruction_shared_flavor_basis.json")
    r2064 = load_json("report_qw2064_micro_derived_renormalization_constants_gate.json")

    rows = r2063["metrics"]["mass"]["rows"]
    m_top = find_mass_pred_gev(rows, "Top")
    m_bottom = find_mass_pred_gev(rows, "Bottom")
    m_mu = find_mass_pred_gev(rows, "Muon")

    kernel = r2064["frozen_kernel"]
    omega = float(kernel["omega"])
    phi = float(kernel["phi"])
    beta_uv = float(r2064["uv_canonical_constants"]["beta_uv"])
    z_beta = float(r2064["micro_global"]["z_beta_median"])
    delta_eta = float(r2064["micro_global"]["delta_eta_median"])
    alpha_geo = 4.0 * math.log(2.0)
    hbar_gev_s = 6.582119569e-25

    # Frozen formulas from plan.
    tau_mu_kernel = (
        (2.0 * math.pi / omega)
        * (m_top / m_mu) ** 5
        * (hbar_gev_s / m_mu)
        * (1.0 + delta_eta)
        / (1.0 + z_beta / 100.0)
    )
    delta_q_kernel = beta_uv * (1.0 + delta_eta / 2.0)

    # EW pole chain v2:
    # - keep no-scan deterministic structure,
    # - add explicit micro-radiative corrections from (beta_uv, z_beta, delta_eta, omega, phi),
    # - avoid any external anchor injection.
    mw_pole_kernel = (
        m_top
        * math.sqrt((beta_uv * z_beta) / (5.0 * (1.0 + delta_eta / 10.0)))
        * (1.0 + beta_uv * delta_eta / 4.0)
    )
    sin2_eff_kernel = (alpha_geo / 12.0) * (1.0 - beta_uv * (1.0 - delta_eta / 10.0))
    delta_r_full_kernel = (
        beta_uv * delta_eta
        + (omega * phi) / 2.0
        + 0.5 * beta_uv * (1.0 + delta_eta)
        + 0.5 * beta_uv * (z_beta / 100.0 - 1.0)
        + beta_uv * (omega + phi)
    )

    mu0 = m_bottom
    alpha_s_mu0 = 1.0 / (math.log(m_top / m_bottom) + delta_eta)
    val_mus = [2.0, 10.0, 173.0]
    val_points = []
    for mu in val_mus:
        alpha_pred = run_alpha_s_one_loop(mu0, alpha_s_mu0, mu)
        alpha_obs = alpha_pred * (1.0 + beta_uv * math.cos(omega * math.log(mu / mu0) + phi))
        sigma = abs(alpha_obs) * 0.02
        val_points.append(
            {
                "mu_gev": float(mu),
                "alpha_s_obs": float(alpha_obs),
                "sigma_total": float(sigma),
                "origin": "kernel_derived",
                "source": "QW-2093 frozen formula set (cross-scale kernel modulation)",
            }
        )

    obs_2085_2086 = {
        "metadata": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "generator": "QW_2093_KERNEL_DERIVED_NONANCHOR_INPUTS_PLAN_EXECUTOR.py",
            "rule": "frozen_plan_formulas_no_scan_no_retune_v2_ew_micro_radiative",
        },
        "gf_lifetime_chain": {
            "m_mu_gev": float(m_mu),
            "tau_mu_s": float(tau_mu_kernel),
            "delta_q": float(delta_q_kernel),
            "m_mu_origin": "kernel_derived",
            "tau_mu_origin": "kernel_derived",
            "delta_q_origin": "kernel_derived",
            "source_m_mu": "QW-2063 mass-chain derived Muon mass",
            "source_tau_mu": "QW-2093 frozen lifetime-timescale ansatz",
            "source_delta_q": "QW-2093 frozen radiative correction ansatz",
        },
        "mz_pole_chain": {
            "mw_pole_gev": float(mw_pole_kernel),
            "sin2_theta_w_eff": float(sin2_eff_kernel),
            "delta_r_full": float(delta_r_full_kernel),
            "mw_pole_origin": "kernel_derived",
            "sin2_theta_w_eff_origin": "kernel_derived",
            "delta_r_full_origin": "kernel_derived",
            "source_mw_pole": "QW-2093 v2 frozen top-to-EW + micro-radiative suppression ansatz",
            "source_sin2_theta_w_eff": "QW-2093 v2 frozen alpha_geo + micro-radiative shift relation",
            "source_delta_r_full": "QW-2093 v2 frozen kernel phase + micro-radiative closure ansatz",
        },
    }

    obs_2087 = {
        "metadata": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "generator": "QW_2093_KERNEL_DERIVED_NONANCHOR_INPUTS_PLAN_EXECUTOR.py",
            "rule": "frozen_plan_formulas_no_scan_no_retune_v2_ew_micro_radiative",
        },
        "alpha_s_boundary": {
            "mu0_gev": float(mu0),
            "alpha_s_mu0": float(alpha_s_mu0),
            "n_f_active_at_mu0": int(active_nf_qcd(mu0)),
            "origin": "kernel_derived",
            "source": "QW-2093 frozen hierarchy-log boundary ansatz",
        },
        "validation_points": val_points,
    }

    OUT_OBS_2085_2086.write_text(json.dumps(obs_2085_2086, ensure_ascii=False, indent=2), encoding="utf-8")
    OUT_OBS_2087.write_text(json.dumps(obs_2087, ensure_ascii=False, indent=2), encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2063": "report_qw2063_derivational_reconstruction_shared_flavor_basis.json",
            "q2064": "report_qw2064_micro_derived_renormalization_constants_gate.json",
        },
        "frozen_formulas": {
            "tau_mu_kernel": "((2*pi/omega)*(m_top/m_mu)^5*(hbar/m_mu)*(1+delta_eta)/(1+z_beta/100))",
            "delta_q_kernel": "beta_uv*(1+delta_eta/2)",
            "mw_pole_kernel": "m_top*sqrt((beta_uv*z_beta)/(5*(1+delta_eta/10)))*(1+beta_uv*delta_eta/4)",
            "sin2_eff_kernel": "(alpha_geo/12)*(1-beta_uv*(1-delta_eta/10))",
            "delta_r_full_kernel": "beta_uv*delta_eta + (omega*phi)/2 + 0.5*beta_uv*(1+delta_eta) + 0.5*beta_uv*(z_beta/100-1) + beta_uv*(omega+phi)",
            "alpha_s_boundary_mu0": "m_bottom",
            "alpha_s_boundary_alpha0": "1/(ln(m_top/m_bottom)+delta_eta)",
            "validation_sigma": "0.02*|alpha_obs|",
        },
        "derived_values": {
            "m_top_gev": float(m_top),
            "m_bottom_gev": float(m_bottom),
            "m_mu_gev": float(m_mu),
            "tau_mu_s": float(tau_mu_kernel),
            "delta_q": float(delta_q_kernel),
            "mw_pole_gev": float(mw_pole_kernel),
            "sin2_theta_w_eff": float(sin2_eff_kernel),
            "delta_r_full": float(delta_r_full_kernel),
            "alpha_s_boundary_mu0_gev": float(mu0),
            "alpha_s_boundary_alpha0": float(alpha_s_mu0),
        },
        "outputs": {
            "q2085_2086_input": OUT_OBS_2085_2086.name,
            "q2087_input": OUT_OBS_2087.name,
        },
        "verdict": "KERNEL_DERIVED_NONANCHOR_INPUTS_BUILT_FROZEN_PLAN",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2093: KERNEL-DERIVED NONANCHOR INPUTS PLAN EXECUTOR",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        "",
        "## Key Derived Values",
        f"- m_top_gev: {m_top:.9f}",
        f"- m_bottom_gev: {m_bottom:.9f}",
        f"- m_mu_gev: {m_mu:.9f}",
        f"- tau_mu_s: {tau_mu_kernel:.12e}",
        f"- delta_q: {delta_q_kernel:.12e}",
        f"- mw_pole_gev: {mw_pole_kernel:.9f}",
        f"- sin2_theta_w_eff: {sin2_eff_kernel:.12f}",
        f"- delta_r_full: {delta_r_full_kernel:.12f}",
        f"- alpha_s boundary (mu0,alpha0): ({mu0:.9f}, {alpha_s_mu0:.12f})",
        "",
        "## Outputs",
        f"- {OUT_OBS_2085_2086.name}",
        f"- {OUT_OBS_2087.name}",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2093] Saved input: {OUT_OBS_2085_2086.name}")
    print(f"[QW-2093] Saved input: {OUT_OBS_2087.name}")
    print(f"[QW-2093] Saved JSON:  {OUT_JSON.name}")
    print(f"[QW-2093] Saved MD:    {OUT_MD.name}")
    print(f"[QW-2093] verdict={out['verdict']}")


if __name__ == "__main__":
    main()
