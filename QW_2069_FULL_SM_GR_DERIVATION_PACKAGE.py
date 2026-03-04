#!/usr/bin/env python3
"""
QW-2069: Full SM+GR derivation package (strict internal snapshot).

This script assembles all currently available first-principles/internal outputs
into one explicit parameter-by-parameter package and quantifies coverage.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2069_full_sm_gr_derivation_package.json"
OUT_MD = ROOT / "RAPORT_QW2069_FULL_SM_GR_DERIVATION_PACKAGE.md"


def rel_err_pct(pred: float, ref: float) -> float:
    return abs(pred - ref) / max(abs(ref), 1e-15) * 100.0


def mean_matrix_rel_pct(pred: np.ndarray, ref: np.ndarray) -> float:
    rel = np.abs(pred - ref) / np.clip(ref, 1e-12, None)
    return float(100.0 * np.mean(rel))


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def load_optional_json(name: str) -> Dict | None:
    path = ROOT / name
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    reg = load_json("report_qw2068_sm_gr_parameter_registry.json")
    r2063 = load_json("report_qw2063_derivational_reconstruction_shared_flavor_basis.json")
    r2064 = load_json("report_qw2064_micro_derived_renormalization_constants_gate.json")
    r2067 = load_json("report_qw2067_strict_first_principles_internal_closure_strengthened_gate.json")
    r2074 = load_optional_json("report_qw2074_strict_nofit_missing_parameter_derivations.json")
    r2075 = load_optional_json("report_qw2075_strict_cp_phase_derivation_gate.json")

    # Canonical bridge constants (project-level formula anchors, not strict independent derivation claims).
    alpha_geo = 4.0 * math.log(2.0)
    beta_uv = 0.01
    alpha_em_inv_pred = (alpha_geo / (2.0 * beta_uv)) * (1.0 - beta_uv)
    sin2_theta_w_pred = alpha_geo / 12.0
    gravity_hierarchy_pred = beta_uv**20

    # Pull strict internal triad outputs.
    mass_rows = r2063["metrics"]["mass"]["rows"]
    mass_pred = {r["particle"]: float(r["pred_mev"]) / 1000.0 for r in mass_rows}  # GeV
    mass_ref = {r["particle"]: float(r["exp_mev"]) / 1000.0 for r in mass_rows}

    ckm_pred = np.array(r2063["metrics"]["flavor"]["ckm_pred_abs"], dtype=float)
    pmns_pred = np.array(r2063["metrics"]["flavor"]["pmns_pred_abs"], dtype=float)
    ckm_ref = np.array(reg["groups"]["flavor"][0]["value"], dtype=float)
    pmns_ref = np.array(reg["groups"]["flavor"][1]["value"], dtype=float)

    # Reference map from registry.
    ref_map: Dict[str, Dict] = {}
    for group, items in reg["groups"].items():
        for item in items:
            ref_map[item["id"]] = {**item, "group": group}

    entries: List[Dict] = []

    def add_entry(
        param_id: str,
        pred_value,
        method: str,
        status: str,
        strict_level: str,
        notes: str,
    ) -> None:
        ref = ref_map.get(param_id)
        rel = None
        within = None
        tol = None
        if ref is not None:
            tol = ref.get("tolerance_rel_pct")
            ref_val = ref.get("value")
            if isinstance(ref_val, (int, float)) and ref_val is not None and isinstance(pred_value, (int, float)):
                rel = rel_err_pct(float(pred_value), float(ref_val))
                if isinstance(tol, (int, float)):
                    within = bool(rel <= float(tol))
        entries.append(
            {
                "id": param_id,
                "group": ref["group"] if ref else None,
                "reference_value": ref.get("value") if ref else None,
                "predicted_value": pred_value,
                "unit": ref.get("unit") if ref else None,
                "tolerance_rel_pct": tol,
                "rel_err_pct": rel,
                "within_tolerance": within,
                "method": method,
                "status": status,
                "strict_level": strict_level,
                "notes": notes,
            }
        )

    # Gauge / EW anchors from canonical formulas.
    add_entry(
        "alpha_em_inv_mz",
        alpha_em_inv_pred,
        "canonical_formula_alpha_geo_beta_uv",
        "model_level_estimate",
        "model_formula",
        "Formula anchor used historically; not promoted as independent strict first-principles proof.",
    )
    add_entry(
        "sin2_theta_w_mz",
        sin2_theta_w_pred,
        "canonical_formula_alpha_geo_over_12",
        "model_level_estimate",
        "model_formula",
        "Formula anchor used as internal benchmark.",
    )

    # Masses (strict internal from QW-2063/1961 chain).
    map_mass_ids = {
        "Top": "m_top",
        "Bottom": "m_bottom",
        "Charm": "m_charm",
        "Tau": "m_tau",
        "Muon": "m_muon",
        "Electron": "m_electron",
    }
    for part, pid in map_mass_ids.items():
        add_entry(
            pid,
            float(mass_pred[part]),
            "qw2063_mass_chain_from_qw1961",
            "derived",
            "strict_internal_gate",
            f"Reference exp={mass_ref[part]:.6g} GeV in mass-chain report.",
        )

    # Flavor matrices.
    ckm_mean = mean_matrix_rel_pct(ckm_pred, ckm_ref)
    pmns_mean = mean_matrix_rel_pct(pmns_pred, pmns_ref)
    add_entry(
        "ckm_matrix_abs",
        ckm_pred.tolist(),
        "qw2063_deterministic_flavor_basis",
        "derived",
        "strict_internal_gate",
        f"Mean relative error={ckm_mean:.3f}%.",
    )
    add_entry(
        "pmns_matrix_abs",
        pmns_pred.tolist(),
        "qw2063_deterministic_flavor_basis",
        "derived",
        "strict_internal_gate",
        f"Mean relative error={pmns_mean:.3f}%.",
    )

    # Project bridge constants from QW-2064.
    add_entry(
        "z_beta",
        float(r2064["micro_global"]["z_beta_median"]),
        "qw2064_micro_derived_renormalization_constants_gate",
        "derived_with_warning_reduced",
        "strict_internal_gate",
        "Supported by micro derivation; dispersion tightened in QW-2066.",
    )
    add_entry(
        "delta_eta",
        float(r2064["micro_global"]["delta_eta_median"]),
        "qw2064_micro_derived_renormalization_constants_gate",
        "derived_with_warning_reduced",
        "strict_internal_gate",
        "Supported by micro derivation; compatibility-filtered tightening available.",
    )
    add_entry(
        "gravity_hierarchy_beta20",
        gravity_hierarchy_pred,
        "canonical_beta_uv_power_20",
        "model_level_estimate",
        "model_formula",
        "Project bridge quantity, not a standalone external strict proof.",
    )

    # Explicitly mark currently missing direct derivations.
    missing_ids = [
        "alpha_s_mz",
        "g_f",
        "m_w",
        "m_z",
        "m_h",
        "v_higgs",
        "m_up",
        "m_down",
        "m_strange",
        "m_nu1",
        "m_nu2",
        "m_nu3",
        "delta_cp_ckm",
        "delta_cp_pmns",
        "G_newton",
        "c_light",
        "hbar",
        "lambda_cosmological",
        "h0",
    ]
    for pid in missing_ids:
        add_entry(
            pid,
            None,
            "not_available_in_current_first_principles_chain",
            "missing",
            "not_derived",
            "No direct strict derivation in current package chain.",
        )

    # Optional no-fit round updates (strictly labeled as anchor-dependent/definition constants).
    if r2074 is not None:
        upd_map = {u["id"]: u for u in r2074.get("updates", []) if isinstance(u, dict) and "id" in u}
        for e in entries:
            if e.get("status") != "missing":
                continue
            u = upd_map.get(e["id"])
            if u is None:
                continue
            e["predicted_value"] = u.get("predicted_value")
            e["method"] = u.get("method", e["method"])
            e["status"] = u.get("status", e["status"])
            e["strict_level"] = u.get("strict_level", e["strict_level"])
            e["notes"] = u.get("notes", e["notes"])

            ref_val = e.get("reference_value")
            pred_val = e.get("predicted_value")
            tol = e.get("tolerance_rel_pct")
            if isinstance(ref_val, (int, float)) and isinstance(pred_val, (int, float)):
                rel = rel_err_pct(float(pred_val), float(ref_val))
                e["rel_err_pct"] = rel
                if isinstance(tol, (int, float)):
                    e["within_tolerance"] = bool(rel <= float(tol))

    # Optional strict CP-phase updates (QW-2075).
    if r2075 is not None:
        upd_map = {u["id"]: u for u in r2075.get("updates", []) if isinstance(u, dict) and "id" in u}
        for e in entries:
            if e.get("status") != "missing":
                continue
            u = upd_map.get(e["id"])
            if u is None:
                continue
            e["predicted_value"] = u.get("predicted_value")
            e["method"] = u.get("method", e["method"])
            e["status"] = u.get("status", e["status"])
            e["strict_level"] = u.get("strict_level", e["strict_level"])
            e["notes"] = u.get("notes", e["notes"])

            ref_val = e.get("reference_value")
            pred_val = e.get("predicted_value")
            tol = e.get("tolerance_rel_pct")
            if isinstance(ref_val, (int, float)) and isinstance(pred_val, (int, float)):
                rel = rel_err_pct(float(pred_val), float(ref_val))
                e["rel_err_pct"] = rel
                if isinstance(tol, (int, float)):
                    e["within_tolerance"] = bool(rel <= float(tol))

    # Coverage summary.
    n_total = int(reg["n_total_parameters"])
    n_derived_strict = sum(1 for e in entries if e["strict_level"] == "strict_internal_gate" and e["status"].startswith("derived"))
    n_model_formula = sum(1 for e in entries if e["strict_level"] == "model_formula")
    n_anchor_dependent = sum(1 for e in entries if e["strict_level"] == "physical_relation_anchor_dependent")
    n_definition_constants = sum(1 for e in entries if e["strict_level"] == "si_definition")
    n_missing = sum(1 for e in entries if e["status"] == "missing")

    # Tolerance stats only where numeric comparison exists.
    numeric_comp = [e for e in entries if isinstance(e.get("rel_err_pct"), float)]
    n_numeric = len(numeric_comp)
    n_numeric_within = sum(1 for e in numeric_comp if bool(e.get("within_tolerance")))

    strict_internal_strengthened_pass = bool(r2067.get("strengthened_pass", False))
    full_smgr_package_pass = bool(strict_internal_strengthened_pass and n_missing == 0)

    if full_smgr_package_pass:
        verdict = "FULL_SM_GR_DERIVATION_PACKAGE_PASS"
    elif strict_internal_strengthened_pass:
        verdict = "FULL_SM_GR_DERIVATION_PACKAGE_PARTIAL_STRONG_INTERNAL"
    else:
        verdict = "FULL_SM_GR_DERIVATION_PACKAGE_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "registry": "report_qw2068_sm_gr_parameter_registry.json",
            "strict_internal_closure": "report_qw2067_strict_first_principles_internal_closure_strengthened_gate.json",
            "triad_derivation": "report_qw2063_derivational_reconstruction_shared_flavor_basis.json",
            "micro_constants": "report_qw2064_micro_derived_renormalization_constants_gate.json",
            "strict_nofit_round_optional": (
                "report_qw2074_strict_nofit_missing_parameter_derivations.json"
                if r2074 is not None
                else None
            ),
            "strict_cp_phase_round_optional": (
                "report_qw2075_strict_cp_phase_derivation_gate.json"
                if r2075 is not None
                else None
            ),
        },
        "entries": entries,
        "coverage": {
            "n_total_registry": n_total,
            "n_entries_evaluated": len(entries),
            "n_derived_strict_internal": n_derived_strict,
            "n_model_formula_only": n_model_formula,
            "n_anchor_dependent_nofit": n_anchor_dependent,
            "n_definition_constants": n_definition_constants,
            "n_missing": n_missing,
            "strict_internal_coverage_fraction": float(n_derived_strict / max(n_total, 1)),
            "numeric_comparisons": {
                "n": n_numeric,
                "n_within_tolerance": n_numeric_within,
                "within_tolerance_fraction": float(n_numeric_within / max(n_numeric, 1)),
            },
        },
        "strict_internal_strengthened_pass": strict_internal_strengthened_pass,
        "verdict": verdict,
        "required_next_step": (
            "IMPLEMENT_MISSING_SM_GR_DERIVATIONS_AND_RADIATIVE_PRECISION_CLOSURE"
            if verdict != "FULL_SM_GR_DERIVATION_PACKAGE_PASS"
            else "PROCEED_TO_EXTERNAL_MULTITEAM_CONFIRMATION"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2069: FULL SM+GR DERIVATION PACKAGE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- strict_internal_strengthened_pass: {strict_internal_strengthened_pass}",
        "",
        "## Coverage",
        f"- total registry parameters: {n_total}",
        f"- strict internal derived: {n_derived_strict}",
        f"- model-formula only: {n_model_formula}",
        f"- no-fit anchor-dependent derived: {n_anchor_dependent}",
        f"- SI definition constants mapped: {n_definition_constants}",
        f"- missing direct derivation: {n_missing}",
        f"- strict internal coverage fraction: {out['coverage']['strict_internal_coverage_fraction']:.3f}",
        "",
        "## Numeric Comparison Stats",
        f"- compared: {n_numeric}",
        f"- within tolerance: {n_numeric_within}",
        f"- within tolerance fraction: {out['coverage']['numeric_comparisons']['within_tolerance_fraction']:.3f}",
        "",
        "## Key Derived Items",
        f"- CKM mean rel%: {ckm_mean:.3f}",
        f"- PMNS mean rel%: {pmns_mean:.3f}",
        f"- z_beta (micro median): {float(r2064['micro_global']['z_beta_median']):.6f}",
        f"- delta_eta (micro median): {float(r2064['micro_global']['delta_eta_median']):.6f}",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2069] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2069] Saved MD:   {OUT_MD.name}")
    print(
        f"[QW-2069] verdict={verdict} strict_derived={n_derived_strict}/{n_total} "
        f"missing={n_missing}"
    )


if __name__ == "__main__":
    main()
