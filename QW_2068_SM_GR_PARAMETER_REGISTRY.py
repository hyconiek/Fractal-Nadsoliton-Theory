#!/usr/bin/env python3
"""
QW-2068: SM+GR parameter registry for strict closure package.

Purpose:
- define an explicit target set for "full SM+GR precision derivation package",
- attach source tags and nominal tolerances,
- avoid ambiguity about what "all known values" means in this project.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2068_sm_gr_parameter_registry.json"
OUT_MD = ROOT / "RAPORT_QW2068_SM_GR_PARAMETER_REGISTRY.md"


def main() -> None:
    # Reference values are standard literature-level anchors used for benchmarking.
    # They are not claimed as new derivations here.
    registry = {
        "gauge_and_electroweak": [
            {"id": "alpha_em_inv_mz", "value": 137.035999084, "unit": "dimensionless", "source": "reference_anchor", "tolerance_rel_pct": 1.0},
            {"id": "sin2_theta_w_mz", "value": 0.23122, "unit": "dimensionless", "source": "reference_anchor", "tolerance_rel_pct": 1.0},
            {"id": "alpha_s_mz", "value": 0.1179, "unit": "dimensionless", "source": "reference_anchor", "tolerance_rel_pct": 5.0},
            {"id": "g_f", "value": 1.1663787e-5, "unit": "GeV^-2", "source": "reference_anchor", "tolerance_rel_pct": 5.0},
            {"id": "m_w", "value": 80.377, "unit": "GeV", "source": "reference_anchor", "tolerance_rel_pct": 2.0},
            {"id": "m_z", "value": 91.1876, "unit": "GeV", "source": "reference_anchor", "tolerance_rel_pct": 2.0},
            {"id": "m_h", "value": 125.25, "unit": "GeV", "source": "reference_anchor", "tolerance_rel_pct": 2.0},
            {"id": "v_higgs", "value": 246.22, "unit": "GeV", "source": "reference_anchor", "tolerance_rel_pct": 2.0},
        ],
        "fermion_masses": [
            {"id": "m_top", "value": 173.0, "unit": "GeV", "source": "project_mass_chain_ref", "tolerance_rel_pct": 35.0},
            {"id": "m_bottom", "value": 4.18, "unit": "GeV", "source": "project_mass_chain_ref", "tolerance_rel_pct": 35.0},
            {"id": "m_charm", "value": 1.27, "unit": "GeV", "source": "project_mass_chain_ref", "tolerance_rel_pct": 35.0},
            {"id": "m_tau", "value": 1.7769, "unit": "GeV", "source": "project_mass_chain_ref", "tolerance_rel_pct": 35.0},
            {"id": "m_muon", "value": 0.1057, "unit": "GeV", "source": "project_mass_chain_ref", "tolerance_rel_pct": 35.0},
            {"id": "m_electron", "value": 0.000511, "unit": "GeV", "source": "project_mass_chain_ref", "tolerance_rel_pct": 35.0},
            {"id": "m_up", "value": 0.00216, "unit": "GeV", "source": "reference_anchor", "tolerance_rel_pct": 50.0},
            {"id": "m_down", "value": 0.00467, "unit": "GeV", "source": "reference_anchor", "tolerance_rel_pct": 50.0},
            {"id": "m_strange", "value": 0.093, "unit": "GeV", "source": "reference_anchor", "tolerance_rel_pct": 50.0},
            {"id": "m_nu1", "value": None, "unit": "eV", "source": "unknown_absolute_scale", "tolerance_rel_pct": None},
            {"id": "m_nu2", "value": None, "unit": "eV", "source": "unknown_absolute_scale", "tolerance_rel_pct": None},
            {"id": "m_nu3", "value": None, "unit": "eV", "source": "unknown_absolute_scale", "tolerance_rel_pct": None},
        ],
        "flavor": [
            {"id": "ckm_matrix_abs", "value": [[0.97401, 0.22650, 0.00361], [0.22636, 0.97320, 0.04053], [0.00854, 0.03978, 0.99917]], "unit": "dimensionless", "source": "project_gate_ref", "tolerance_rel_pct": 15.0},
            {"id": "pmns_matrix_abs", "value": [[0.821, 0.550, 0.150], [0.432, 0.582, 0.693], [0.378, 0.598, 0.707]], "unit": "dimensionless", "source": "project_gate_ref", "tolerance_rel_pct": 15.0},
            {"id": "delta_cp_ckm", "value": 1.20, "unit": "rad", "source": "reference_anchor", "tolerance_rel_pct": 20.0},
            {"id": "delta_cp_pmns", "value": None, "unit": "rad", "source": "currently_weakly_constrained", "tolerance_rel_pct": None},
        ],
        "gr_and_cosmology": [
            {"id": "G_newton", "value": 6.67430e-11, "unit": "m^3 kg^-1 s^-2", "source": "reference_anchor", "tolerance_rel_pct": 5.0},
            {"id": "c_light", "value": 299792458.0, "unit": "m s^-1", "source": "exact_si", "tolerance_rel_pct": 0.0},
            {"id": "hbar", "value": 1.054571817e-34, "unit": "J s", "source": "reference_anchor", "tolerance_rel_pct": 5.0},
            {"id": "lambda_cosmological", "value": 1.1056e-52, "unit": "m^-2", "source": "reference_anchor", "tolerance_rel_pct": 20.0},
            {"id": "h0", "value": 67.4, "unit": "km s^-1 Mpc^-1", "source": "reference_anchor", "tolerance_rel_pct": 10.0},
        ],
        "project_specific_bridge_quantities": [
            {"id": "gravity_hierarchy_beta20", "value": 1.0e-40, "unit": "dimensionless", "source": "project_formula_anchor", "tolerance_rel_pct": 50.0},
            {"id": "z_beta", "value": 100.0, "unit": "dimensionless", "source": "project_refrozen_kernel_vs_uv", "tolerance_rel_pct": 100.0},
            {"id": "delta_eta", "value": 0.8, "unit": "dimensionless", "source": "project_refrozen_kernel_vs_uv", "tolerance_rel_pct": 50.0},
        ],
    }

    n_total = sum(len(v) for v in registry.values())

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "registry_version": "1.0",
        "purpose": "full_sm_gr_precision_derivation_target_set",
        "groups": registry,
        "n_total_parameters": int(n_total),
        "notes": [
            "Registry defines closure target set, not a claim that all values are already derived.",
            "Some entries intentionally have value=None where external physics lacks fixed absolute value in this context.",
        ],
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2068: SM+GR PARAMETER REGISTRY",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Registry version: {out['registry_version']}",
        f"- Total target parameters: {n_total}",
        "",
        "## Groups",
    ]
    for g, items in registry.items():
        lines.append(f"- {g}: {len(items)}")

    lines.extend([
        "",
        "## Interpretation",
        "- This registry is the explicit closure target map for the full SM+GR package.",
        "- It separates project-derived quantities from reference-anchor quantities.",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2068] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2068] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2068] n_total_parameters={n_total}")


if __name__ == "__main__":
    main()
