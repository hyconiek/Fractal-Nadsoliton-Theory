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
    denom = abs(ref) if ref != 0.0 else 1e-300
    return abs(pred - ref) / denom * 100.0


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


def is_closure_resolved_entry(e: Dict) -> bool:
    lvl = str(e.get("strict_level", ""))
    st = str(e.get("status", ""))

    # SI definitions are treated as explicitly resolved-by-definition in package logic.
    if lvl == "si_definition" and st == "definition_constant":
        return True

    # Strict internal derivation must not be a target-miss and must satisfy tolerance when numeric.
    if lvl == "strict_internal_gate" and st.startswith("derived") and st != "derived_strict_target_miss":
        within = e.get("within_tolerance")
        if isinstance(within, bool):
            return within
        return True

    return False


def main() -> None:
    reg = load_json("report_qw2068_sm_gr_parameter_registry.json")
    r2063 = load_json("report_qw2063_derivational_reconstruction_shared_flavor_basis.json")
    r2064 = load_json("report_qw2064_micro_derived_renormalization_constants_gate.json")
    r2067 = load_json("report_qw2067_strict_first_principles_internal_closure_strengthened_gate.json")
    r2074 = load_optional_json("report_qw2074_strict_nofit_missing_parameter_derivations.json")
    r2075 = load_optional_json("report_qw2075_strict_cp_phase_derivation_gate.json")
    r2083 = load_optional_json("report_qw2083_missing14_epistemic_status_gate.json")
    r2084 = load_optional_json("report_qw2084_t1_nonanchor_strict_gate.json")
    r2085 = load_optional_json("report_qw2085_gf_nonanchor_lifetime_gate.json")
    r2086 = load_optional_json("report_qw2086_mz_nonanchor_ew_pole_gate.json")
    r2087 = load_optional_json("report_qw2087_alpha_s_nonanchor_boundary_gate.json")
    r2088 = load_optional_json("report_qw2088_light_quark_mass_nonanchor_gate.json")
    r2089 = load_optional_json("report_qw2089_higgs_selfcoupling_strict_gate.json")
    r2090 = load_optional_json("report_qw2090_h0_lambda_decoupling_gate.json")
    r2091 = load_optional_json("report_qw2091_neutrino_absolute_scale_gate.json")
    r2092 = load_optional_json("report_qw2092_gnewton_si_bridge_gate.json")
    r2096 = load_optional_json("report_qw2096_t2_nonanchor_strict_gate.json")
    r2097 = load_optional_json("report_qw2097_ckm_cp_target_refinement_gate.json")
    r2098 = load_optional_json("report_qw2098_ew_secondary_nonanchor_closure_gate.json")
    r2115 = load_optional_json("report_qw2115_gravity_hierarchy_strict_bridge_gate.json")

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

    # Optional missing-14 epistemic status gate updates (QW-2083).
    if r2083 is not None:
        upd_map = {u["id"]: u for u in r2083.get("updates", []) if isinstance(u, dict) and "id" in u}
        for e in entries:
            # QW-2083 may update entries even if they were previously touched by QW-2074/2075.
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

    # Optional T1 strict non-anchor audit updates (QW-2084).
    if r2084 is not None:
        upd_map = {u["id"]: u for u in r2084.get("updates", []) if isinstance(u, dict) and "id" in u}
        for e in entries:
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

    # Optional single-update gates (QW-2085/QW-2086).
    for rsingle in [r2085, r2086, r2087, r2092, r2097, r2115]:
        if rsingle is None:
            continue
        u = rsingle.get("update")
        if not (isinstance(u, dict) and "id" in u):
            continue
        for e in entries:
            if e.get("id") != u.get("id"):
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

    # Optional multi-update T2 gates (QW-2088/QW-2089/QW-2096).
    for rmulti in [r2088, r2089, r2090, r2091, r2096, r2098]:
        if rmulti is None:
            continue
        updates = rmulti.get("updates")
        if not isinstance(updates, list):
            upd_single = rmulti.get("update")
            updates = [upd_single] if isinstance(upd_single, dict) else []
        upd_map = {u["id"]: u for u in updates if isinstance(u, dict) and "id" in u}
        if not upd_map:
            continue
        for e in entries:
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
    n_coupled_anchor_dependent = sum(1 for e in entries if e["strict_level"] == "coupled_anchor_dependent")
    n_model_assumption = sum(1 for e in entries if e["strict_level"] == "model_assumption_anchor")
    n_definition_constants = sum(1 for e in entries if e["strict_level"] == "si_definition")
    n_missing = sum(1 for e in entries if e["status"] == "missing")
    strict_unresolved_ids = sorted([e["id"] for e in entries if not is_closure_resolved_entry(e)])
    n_strict_unresolved = len(strict_unresolved_ids)

    # Tolerance stats only where numeric comparison exists.
    numeric_comp = [e for e in entries if isinstance(e.get("rel_err_pct"), float)]
    n_numeric = len(numeric_comp)
    n_numeric_within = sum(1 for e in numeric_comp if bool(e.get("within_tolerance")))

    strict_internal_strengthened_pass = bool(r2067.get("strengthened_pass", False))
    full_smgr_package_pass = bool(strict_internal_strengthened_pass and n_missing == 0 and n_strict_unresolved == 0)

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
            "missing14_epistemic_status_gate_optional": (
                "report_qw2083_missing14_epistemic_status_gate.json"
                if r2083 is not None
                else None
            ),
            "t1_nonanchor_strict_gate_optional": (
                "report_qw2084_t1_nonanchor_strict_gate.json"
                if r2084 is not None
                else None
            ),
            "gf_nonanchor_lifetime_gate_optional": (
                "report_qw2085_gf_nonanchor_lifetime_gate.json"
                if r2085 is not None
                else None
            ),
            "mz_nonanchor_ew_pole_gate_optional": (
                "report_qw2086_mz_nonanchor_ew_pole_gate.json"
                if r2086 is not None
                else None
            ),
            "alpha_s_nonanchor_boundary_gate_optional": (
                "report_qw2087_alpha_s_nonanchor_boundary_gate.json"
                if r2087 is not None
                else None
            ),
            "light_quark_mass_nonanchor_gate_optional": (
                "report_qw2088_light_quark_mass_nonanchor_gate.json"
                if r2088 is not None
                else None
            ),
            "higgs_selfcoupling_strict_gate_optional": (
                "report_qw2089_higgs_selfcoupling_strict_gate.json"
                if r2089 is not None
                else None
            ),
            "h0_lambda_decoupling_gate_optional": (
                "report_qw2090_h0_lambda_decoupling_gate.json"
                if r2090 is not None
                else None
            ),
            "neutrino_absolute_scale_gate_optional": (
                "report_qw2091_neutrino_absolute_scale_gate.json"
                if r2091 is not None
                else None
            ),
            "gnewton_si_bridge_gate_optional": (
                "report_qw2092_gnewton_si_bridge_gate.json"
                if r2092 is not None
                else None
            ),
            "t2_nonanchor_strict_gate_optional": (
                "report_qw2096_t2_nonanchor_strict_gate.json"
                if r2096 is not None
                else None
            ),
            "ckm_cp_target_refinement_gate_optional": (
                "report_qw2097_ckm_cp_target_refinement_gate.json"
                if r2097 is not None
                else None
            ),
            "ew_secondary_nonanchor_closure_gate_optional": (
                "report_qw2098_ew_secondary_nonanchor_closure_gate.json"
                if r2098 is not None
                else None
            ),
            "gravity_hierarchy_strict_bridge_gate_optional": (
                "report_qw2115_gravity_hierarchy_strict_bridge_gate.json"
                if r2115 is not None
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
            "n_coupled_anchor_dependent": n_coupled_anchor_dependent,
            "n_model_assumption_nonclosing": n_model_assumption,
            "n_definition_constants": n_definition_constants,
            "n_missing": n_missing,
            "n_strict_unresolved": n_strict_unresolved,
            "strict_unresolved_ids": strict_unresolved_ids,
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
        f"- coupled anchor-dependent derived: {n_coupled_anchor_dependent}",
        f"- model-assumption nonclosing derived: {n_model_assumption}",
        f"- SI definition constants mapped: {n_definition_constants}",
        f"- missing direct derivation: {n_missing}",
        f"- strict unresolved (closure criterion): {n_strict_unresolved}",
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
