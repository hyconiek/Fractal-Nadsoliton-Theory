#!/usr/bin/env python3
"""
QW-2095: Kernel-derived T2 non-anchor inputs plan executor (QW-2088..QW-2089).
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_Q2088 = ROOT / "t2_nonanchor_light_quark_input_qw2088.json"
OUT_Q2089 = ROOT / "t2_nonanchor_higgs_input_qw2089.json"
OUT_JSON = ROOT / "report_qw2095_kernel_derived_t2_nonanchor_inputs_plan_executor.json"
OUT_MD = ROOT / "RAPORT_QW2095_KERNEL_DERIVED_T2_NONANCHOR_INPUTS_PLAN_EXECUTOR.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def find_mass_pred_gev(rows: List[Dict], particle: str) -> float:
    for r in rows:
        if str(r.get("particle")) == particle:
            return float(r["pred_mev"]) / 1000.0
    raise KeyError(f"Missing particle in QW-2063 mass rows: {particle}")


def lambda_run_proxy(lambda0: float, mu0: float, mu: float, beta_uv: float) -> float:
    if lambda0 <= 0.0 or mu0 <= 0.0 or mu <= 0.0:
        raise ValueError("lambda0/mu0/mu must be positive.")
    return float(lambda0 / (1.0 + beta_uv * math.log(mu / mu0)))


def main() -> None:
    r2063 = load_json("report_qw2063_derivational_reconstruction_shared_flavor_basis.json")
    r2064 = load_json("report_qw2064_micro_derived_renormalization_constants_gate.json")

    rows = r2063["metrics"]["mass"]["rows"]
    m_top = find_mass_pred_gev(rows, "Top")

    kernel = r2064["frozen_kernel"]
    omega = float(kernel["omega"])
    phi = float(kernel["phi"])
    beta_uv = float(r2064["uv_canonical_constants"]["beta_uv"])
    z_beta = float(r2064["micro_global"]["z_beta_median"])
    delta_eta = float(r2064["micro_global"]["delta_eta_median"])

    alpha_geo = 4.0 * math.log(2.0)

    # Frozen T2 formulas from plan.
    m_up = m_top * (beta_uv**2.5) * (1.0 + delta_eta / 10.0)
    m_down = m_up * (1.0 + delta_eta)
    m_strange = m_down * (beta_uv ** (-0.6)) * (1.0 + omega)

    lambda_eff = (alpha_geo / (8.0 * math.pi)) * (1.0 + delta_eta / 10.0 + beta_uv + omega / 2.0)
    v_eff = (m_top / (omega * math.sqrt(z_beta / 100.0))) * (beta_uv**0.25) * (1.0 - delta_eta / 10.0)
    m_h = math.sqrt(max(2.0 * lambda_eff, 1e-15)) * v_eff

    q2088_in = {
        "metadata": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "generator": "QW_2095_KERNEL_DERIVED_T2_NONANCHOR_INPUTS_PLAN_EXECUTOR.py",
            "rule": "frozen_plan_formulas_no_scan_no_retune",
        },
        "light_quark_chain": {
            "m_up_gev": float(m_up),
            "m_down_gev": float(m_down),
            "m_strange_gev": float(m_strange),
            "m_up_origin": "kernel_derived",
            "m_down_origin": "kernel_derived",
            "m_strange_origin": "kernel_derived",
            "source_m_up": "QW-2095 frozen light-quark hierarchy ansatz",
            "source_m_down": "QW-2095 frozen light-quark hierarchy ansatz",
            "source_m_strange": "QW-2095 frozen light-quark hierarchy ansatz",
        },
        "hadronic_crosschecks": {
            "ratio_md_mu_obs": 2.16,
            "ratio_md_mu_sigma": 0.4,
            "ratio_ms_md_obs": 19.9,
            "ratio_ms_md_sigma": 3.0,
            "origin": "external_hadronic",
            "source": "PDG/FLAG light-quark mass ratio summaries (independent hadronic constraints)",
        },
    }

    q2089_in = {
        "metadata": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "generator": "QW_2095_KERNEL_DERIVED_T2_NONANCHOR_INPUTS_PLAN_EXECUTOR.py",
            "rule": "frozen_plan_formulas_no_scan_no_retune",
        },
        "higgs_chain": {
            "lambda_eff": float(lambda_eff),
            "v_eff_gev": float(v_eff),
            "m_h_gev": float(m_h),
            "lambda_eff_origin": "kernel_derived",
            "v_eff_origin": "kernel_derived",
            "m_h_origin": "kernel_derived",
            "source_lambda_eff": "QW-2095 frozen Higgs self-coupling ansatz",
            "source_v_eff": "QW-2095 frozen electroweak-scale ansatz",
            "source_m_h": "QW-2095 frozen Higgs mass closure ansatz",
        },
        "lambda_validation_points": [
            {
                "mu_gev": 125.0,
                "lambda_obs": 0.129,
                "sigma_total": 0.02,
                "origin": "external_higgs",
                "source": "Independent Higgs self-coupling summary",
            },
            {
                "mu_gev": 500.0,
                "lambda_obs": 0.120,
                "sigma_total": 0.03,
                "origin": "external_higgs",
                "source": "Independent Higgs running summary",
            },
        ],
    }

    OUT_Q2088.write_text(json.dumps(q2088_in, ensure_ascii=False, indent=2), encoding="utf-8")
    OUT_Q2089.write_text(json.dumps(q2089_in, ensure_ascii=False, indent=2), encoding="utf-8")

    ratio_md_mu = m_down / max(m_up, 1e-15)
    ratio_ms_md = m_strange / max(m_down, 1e-15)
    lambda_125 = lambda_run_proxy(lambda_eff, 125.0, 125.0, beta_uv)
    lambda_500 = lambda_run_proxy(lambda_eff, 125.0, 500.0, beta_uv)

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2063": "report_qw2063_derivational_reconstruction_shared_flavor_basis.json",
            "q2064": "report_qw2064_micro_derived_renormalization_constants_gate.json",
        },
        "frozen_formulas": {
            "m_up": "m_top * beta_uv^(2.5) * (1 + delta_eta/10)",
            "m_down": "m_up * (1 + delta_eta)",
            "m_strange": "m_down * beta_uv^(-0.6) * (1 + omega)",
            "lambda_eff": "(alpha_geo/(8*pi)) * (1 + delta_eta/10 + beta_uv + omega/2)",
            "v_eff": "(m_top/(omega*sqrt(z_beta/100))) * beta_uv^(0.25) * (1 - delta_eta/10)",
            "m_h": "sqrt(2*lambda_eff) * v_eff",
        },
        "derived_values": {
            "m_up_gev": float(m_up),
            "m_down_gev": float(m_down),
            "m_strange_gev": float(m_strange),
            "ratio_md_mu_pred": float(ratio_md_mu),
            "ratio_ms_md_pred": float(ratio_ms_md),
            "lambda_eff": float(lambda_eff),
            "v_eff_gev": float(v_eff),
            "m_h_gev": float(m_h),
            "lambda_125_pred": float(lambda_125),
            "lambda_500_pred": float(lambda_500),
            "omega": float(omega),
            "phi": float(phi),
            "beta_uv": float(beta_uv),
            "z_beta": float(z_beta),
            "delta_eta": float(delta_eta),
        },
        "outputs": {
            "q2088_input": OUT_Q2088.name,
            "q2089_input": OUT_Q2089.name,
        },
        "verdict": "KERNEL_DERIVED_T2_NONANCHOR_INPUTS_BUILT_FROZEN_PLAN",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2095: KERNEL-DERIVED T2 NONANCHOR INPUTS PLAN EXECUTOR",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        "",
        "## Derived Values",
        f"- m_up_gev: {m_up:.9e}",
        f"- m_down_gev: {m_down:.9e}",
        f"- m_strange_gev: {m_strange:.9e}",
        f"- ratio_md_mu_pred: {ratio_md_mu:.6f}",
        f"- ratio_ms_md_pred: {ratio_ms_md:.6f}",
        f"- lambda_eff: {lambda_eff:.9f}",
        f"- v_eff_gev: {v_eff:.9f}",
        f"- m_h_gev: {m_h:.9f}",
        "",
        "## Outputs",
        f"- {OUT_Q2088.name}",
        f"- {OUT_Q2089.name}",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2095] Saved input: {OUT_Q2088.name}")
    print(f"[QW-2095] Saved input: {OUT_Q2089.name}")
    print(f"[QW-2095] Saved JSON:  {OUT_JSON.name}")
    print(f"[QW-2095] Saved MD:    {OUT_MD.name}")
    print(f"[QW-2095] verdict={out['verdict']}")


if __name__ == "__main__":
    main()

