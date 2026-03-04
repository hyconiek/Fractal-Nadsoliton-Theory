#!/usr/bin/env python3
"""
QW-2081: Missing-14 strict-rigor frontier audit.

Goal:
- evaluate what can/cannot be derived for the current 14 missing SM+GR parameters
  under strict constraints (deterministic, no retune, no scan),
- separate strict candidates from anchor-dependent baselines and truly underdetermined
  quantities.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2081_missing14_strict_rigor_frontier.json"
OUT_MD = ROOT / "RAPORT_QW2081_MISSING14_STRICT_RIGOR_FRONTIER.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def find_entry(entries: List[Dict], pid: str) -> Dict:
    for e in entries:
        if e.get("id") == pid:
            return e
    raise KeyError(f"Missing entry id in QW-2069: {pid}")


def is_closure_resolved_entry(e: Dict) -> bool:
    lvl = str(e.get("strict_level", ""))
    st = str(e.get("status", ""))
    if lvl == "si_definition" and st == "definition_constant":
        return True
    if lvl == "strict_internal_gate" and st.startswith("derived") and st != "derived_strict_target_miss":
        within = e.get("within_tolerance")
        if isinstance(within, bool):
            return within
        return True
    return False


def min_rel_err_branch(ckm_phase: Dict) -> tuple[float, float]:
    ref = float(ckm_phase["registry_reference_rad"])
    d1 = float(ckm_phase["delta_primary_rad"])
    d2 = float(ckm_phase["delta_secondary_rad"])
    e1 = abs(d1 - ref) / max(abs(ref), 1e-15) * 100.0
    e2 = abs(d2 - ref) / max(abs(ref), 1e-15) * 100.0
    return (d1, e1) if e1 <= e2 else (d2, e2)


def main() -> None:
    r2069 = load_json("report_qw2069_full_sm_gr_derivation_package.json")
    r2070 = load_json("report_qw2070_full_radiative_program_baseline.json")
    r2072 = load_json("report_qw2072_ew_yukawa_flavor_radiative_baselines.json")
    r2073 = load_json("report_qw2073_radiative_channels_closure_upgrade.json")
    r2075 = load_json("report_qw2075_strict_cp_phase_derivation_gate.json")

    # Keep the tracked historical 14 list explicit and auditable.
    expected_missing_14 = [
        "G_newton",
        "alpha_s_mz",
        "delta_cp_ckm",
        "g_f",
        "h0",
        "lambda_cosmological",
        "m_down",
        "m_h",
        "m_nu1",
        "m_nu2",
        "m_nu3",
        "m_strange",
        "m_up",
        "m_z",
    ]

    entries = r2069["entries"]
    unresolved_ids = sorted([e["id"] for e in entries if not is_closure_resolved_entry(e)])
    tracked_set = set(expected_missing_14)
    unresolved_in_scope = sorted([x for x in unresolved_ids if x in tracked_set])
    resolved_in_scope = sorted([x for x in expected_missing_14 if x not in set(unresolved_in_scope)])

    ew = r2072["electroweak_baseline"]
    yuk = r2072["yukawa_baseline"]
    unresolved_set = set(unresolved_in_scope)

    # Deterministic baseline-derived quantities (explicitly anchor-dependent).
    v = float(ew["v_from_gf"])
    g_f_from_v = float(1.0 / (math.sqrt(2.0) * v * v))
    alpha_s_from_qcd_baseline = None
    mz_ref = float(ew["mz_ref"])
    for row in r2070["qcd_one_loop_baseline"]["samples"]:
        if abs(float(row["mu_gev"]) - mz_ref) < 1e-9:
            alpha_s_from_qcd_baseline = float(row["alpha_s"])
            break
    if alpha_s_from_qcd_baseline is None:
        raise RuntimeError("Could not find alpha_s sample at mu=MZ in QW-2070 output.")

    m_z_from_ew = float(ew["m_w_delta_r0_from_alpha_gf_mz"]) / math.sqrt(max(1.0 - float(ew["sin2_theta_w_input"]), 1e-15))
    y_top = float(yuk["y_top"])
    g = float(ew["g"])
    gp = float(ew["g_prime"])
    # Conservative one-loop inspired proxy with fixed mu=1 TeV (deterministic, no scan).
    mu_loop = 1000.0
    m_top_gev = 173.0
    lambda_tree = (g * g + gp * gp) / 8.0
    delta_lambda_top = (3.0 * y_top**4 / (8.0 * math.pi * math.pi)) * math.log(mu_loop / m_top_gev)
    lambda_eff = max(lambda_tree + delta_lambda_top, 1e-9)
    m_h_proxy = math.sqrt(2.0 * lambda_eff) * v

    m_up_from_y = float(yuk["y_up"]) * v / math.sqrt(2.0)
    m_down_from_y = float(yuk["y_down"]) * v / math.sqrt(2.0)
    m_strange_from_y = float(yuk["y_strange"]) * v / math.sqrt(2.0)

    # CKM CP strict candidate from QW-2075 deterministic phase extraction.
    ckm_best_branch, ckm_rel_err = min_rel_err_branch(r2075["ckm_phase"])
    ckm_tol = float(r2075["ckm_phase"]["registry_tolerance_rel_pct"])

    # Coupled cosmology relation: lambda <-> H0 given Omega_Lambda; still needs one anchor.
    omega_lambda_flat = float(r2073["diagnostics"]["omega_lambda_flat"])
    c_light = 299792458.0
    mpc_m = 3.085677581e22
    lambda_ref = float(find_entry(entries, "lambda_cosmological")["reference_value"])
    h0_from_lambda = math.sqrt(lambda_ref * c_light * c_light / (3.0 * omega_lambda_flat)) * (mpc_m / 1000.0)
    h0_entry = find_entry(entries, "h0")
    lam_entry = find_entry(entries, "lambda_cosmological")

    # Anchor-dependent SI bridge for Newton constant from dimensionless GR sample.
    # NOTE: this remains non-closing because the running baseline itself is anchored.
    gr_samples = r2073["diagnostics"]["gr_samples"]
    g_at_1gev = None
    for s in gr_samples:
        if abs(float(s["mu_gev"]) - 1.0) < 1e-12:
            g_at_1gev = float(s["g_dimensionless"])
            break
    if g_at_1gev is None:
        raise RuntimeError("Could not find GR sample at mu=1 GeV in QW-2073 diagnostics.")
    g_newton_nat_gev_m2 = g_at_1gev
    hbar_si = 1.054571817e-34  # J*s
    c_si = 299792458.0  # m/s
    gev_to_j = 1.602176634e-10  # J
    g_newton_si_from_anchor_bridge = (
        g_newton_nat_gev_m2 * hbar_si * (c_si**5) / (gev_to_j**2)
    )

    # Optional model-assumption neutrino hierarchy (not strict closure).
    dm21 = 7.42e-5  # eV^2 (anchor assumption)
    dm31 = 2.517e-3  # eV^2 (anchor assumption, normal hierarchy)
    m_nu1_assumed = 0.0
    m_nu2_assumed = math.sqrt(dm21)
    m_nu3_assumed = math.sqrt(dm31)

    frontier: List[Dict] = []

    if "delta_cp_ckm" in unresolved_set:
        frontier.append(
            {
                "id": "delta_cp_ckm",
                "classification": (
                    "strict_candidate_target_miss"
                    if ckm_rel_err > ckm_tol
                    else "strict_candidate_within_tolerance"
                ),
                "proposed_value": float(ckm_best_branch),
                "strict_level_proposed": "strict_internal_gate",
                "status_proposed": "derived_strict_target_miss" if ckm_rel_err > ckm_tol else "derived",
                "rel_err_pct_vs_registry": float(ckm_rel_err),
                "tolerance_rel_pct": float(ckm_tol),
                "method": "qw2075_deterministic_complex_ckm_jarlskog_phase_best_branch",
                "notes": "Deterministic no-scan CKM CP phase is available but currently outside registry tolerance.",
            }
        )

    if "alpha_s_mz" in unresolved_set:
        frontier.append(
            {
                "id": "alpha_s_mz",
                "classification": "anchor_dependent_baseline_only",
                "proposed_value": float(alpha_s_from_qcd_baseline),
                "strict_level_proposed": "model_radiative_precision_baseline",
                "status_proposed": "derived_anchor_dependent_baseline",
                "method": "qw2070_qcd_one_loop_sample_at_mz",
                "notes": "Computed from baseline seeded at registry anchor alpha_s(MZ); not independent first-principles closure.",
            }
        )

    if "g_f" in unresolved_set:
        frontier.append(
            {
                "id": "g_f",
                "classification": "anchor_dependent_baseline_only",
                "proposed_value": float(g_f_from_v),
                "strict_level_proposed": "physical_relation_anchor_dependent",
                "status_proposed": "derived_anchor_dependent",
                "method": "G_F=1/(sqrt(2)*v^2) with v from QW-2072 EW baseline",
                "notes": "Identity-level relation under anchored v; no independent strict derivation.",
            }
        )

    if "m_z" in unresolved_set:
        frontier.append(
            {
                "id": "m_z",
                "classification": "anchor_dependent_baseline_only",
                "proposed_value": float(m_z_from_ew),
                "strict_level_proposed": "model_radiative_precision_baseline",
                "status_proposed": "derived_anchor_dependent_baseline",
                "method": "MZ=MW(delta_r=0)/sqrt(1-sin^2(theta_W)) from QW-2072",
                "notes": "Deterministic EW baseline relation, still anchor-dependent.",
            }
        )

    if "m_h" in unresolved_set:
        frontier.append(
            {
                "id": "m_h",
                "classification": "model_assumption_required",
                "proposed_value": float(m_h_proxy),
                "strict_level_proposed": "model_assumption_anchor",
                "status_proposed": "derived_model_assumption_nonclosing",
                "method": "lambda_eff=lambda_tree+top_loop(mu=1TeV), m_h=sqrt(2*lambda_eff)*v",
                "notes": "Requires explicit model assumptions beyond current strict closure chain.",
            }
        )

    if "m_up" in unresolved_set:
        frontier.append(
            {
                "id": "m_up",
                "classification": "anchor_dependent_baseline_only",
                "proposed_value": float(m_up_from_y),
                "strict_level_proposed": "model_radiative_precision_baseline",
                "status_proposed": "derived_anchor_dependent_baseline",
                "method": "m=y*v/sqrt(2) using QW-2072 Yukawa baseline",
                "notes": "Back-computed from Yukawa baseline anchored to reference masses.",
            }
        )

    if "m_down" in unresolved_set:
        frontier.append(
            {
                "id": "m_down",
                "classification": "anchor_dependent_baseline_only",
                "proposed_value": float(m_down_from_y),
                "strict_level_proposed": "model_radiative_precision_baseline",
                "status_proposed": "derived_anchor_dependent_baseline",
                "method": "m=y*v/sqrt(2) using QW-2072 Yukawa baseline",
                "notes": "Back-computed from Yukawa baseline anchored to reference masses.",
            }
        )

    if "m_strange" in unresolved_set:
        frontier.append(
            {
                "id": "m_strange",
                "classification": "anchor_dependent_baseline_only",
                "proposed_value": float(m_strange_from_y),
                "strict_level_proposed": "model_radiative_precision_baseline",
                "status_proposed": "derived_anchor_dependent_baseline",
                "method": "m=y*v/sqrt(2) using QW-2072 Yukawa baseline",
                "notes": "Back-computed from Yukawa baseline anchored to reference masses.",
            }
        )

    if "h0" in unresolved_set:
        h0_status = str(h0_entry.get("status", ""))
        h0_is_target_miss = h0_status == "derived_strict_target_miss"
        h0_proposed = h0_entry.get("predicted_value", h0_from_lambda) if h0_is_target_miss else h0_from_lambda
        frontier.append(
            {
                "id": "h0",
                "classification": "strict_candidate_target_miss" if h0_is_target_miss else "coupled_underdetermined_pair",
                "proposed_value": float(h0_proposed),
                "strict_level_proposed": (
                    "strict_internal_gate" if h0_is_target_miss else "physical_relation_anchor_dependent"
                ),
                "status_proposed": "derived_strict_target_miss" if h0_is_target_miss else "derived_coupled_anchor_dependent",
                "method": (
                    str(h0_entry.get("method", "qw2090_weighted_hz_decoupling_two_parameter_fit_nonanchor"))
                    if h0_is_target_miss
                    else "H0 from lambda and Omega_Lambda flat relation"
                ),
                "notes": (
                    "Strict decoupling candidate exists but misses registry tolerance."
                    if h0_is_target_miss
                    else "Requires lambda anchor; not standalone strict first-principles derivation."
                ),
            }
        )

    if "lambda_cosmological" in unresolved_set:
        lam_status = str(lam_entry.get("status", ""))
        lam_is_target_miss = lam_status == "derived_strict_target_miss"
        lam_proposed = lam_entry.get("predicted_value", lambda_ref) if lam_is_target_miss else lambda_ref
        frontier.append(
            {
                "id": "lambda_cosmological",
                "classification": (
                    "strict_candidate_target_miss" if lam_is_target_miss else "coupled_underdetermined_pair"
                ),
                "proposed_value": float(lam_proposed),
                "strict_level_proposed": (
                    "strict_internal_gate" if lam_is_target_miss else "physical_relation_anchor_dependent"
                ),
                "status_proposed": (
                    "derived_strict_target_miss" if lam_is_target_miss else "derived_coupled_anchor_dependent"
                ),
                "method": (
                    str(lam_entry.get("method", "qw2090_weighted_hz_decoupling_two_parameter_fit_nonanchor"))
                    if lam_is_target_miss
                    else "Lambda-H0-Omega relation (paired with H0)"
                ),
                "notes": (
                    "Strict decoupling candidate exists but misses registry tolerance."
                    if lam_is_target_miss
                    else "Absolute closure remains coupled to external anchor choice."
                ),
            }
        )

    if "G_newton" in unresolved_set:
        frontier.append(
            {
                "id": "G_newton",
                "classification": "anchor_dependent_baseline_only",
                "proposed_value": float(g_newton_si_from_anchor_bridge),
                "strict_level_proposed": "physical_relation_anchor_dependent",
                "status_proposed": "derived_anchor_dependent_baseline",
                "method": "G_nat(mu=1GeV)->SI bridge using hbar*c^5/GeV^2 from QW-2073 gr_samples",
                "notes": (
                    "Numerically mapped via anchored GR running baseline; explicit non-closing "
                    "anchor-dependent status (not independent strict first-principles closure)."
                ),
            }
        )

    if "m_nu1" in unresolved_set:
        frontier.append(
            {
                "id": "m_nu1",
                "classification": "model_assumption_required",
                "proposed_value": float(m_nu1_assumed),
                "strict_level_proposed": "model_assumption_anchor",
                "status_proposed": "derived_model_assumption_nonclosing",
                "method": "normal_hierarchy_assumption_m1_zero",
                "notes": "Absolute neutrino scale is underdetermined in current strict chain.",
            }
        )

    if "m_nu2" in unresolved_set:
        frontier.append(
            {
                "id": "m_nu2",
                "classification": "model_assumption_required",
                "proposed_value": float(m_nu2_assumed),
                "strict_level_proposed": "model_assumption_anchor",
                "status_proposed": "derived_model_assumption_nonclosing",
                "method": "m2=sqrt(dm21) with anchored dm21",
                "notes": "Depends on external oscillation anchor assumptions.",
            }
        )

    if "m_nu3" in unresolved_set:
        frontier.append(
            {
                "id": "m_nu3",
                "classification": "model_assumption_required",
                "proposed_value": float(m_nu3_assumed),
                "strict_level_proposed": "model_assumption_anchor",
                "status_proposed": "derived_model_assumption_nonclosing",
                "method": "m3=sqrt(dm31) with anchored dm31",
                "notes": "Depends on external oscillation anchor assumptions.",
            }
        )

    by_class = {}
    for r in frontier:
        by_class[r["classification"]] = by_class.get(r["classification"], 0) + 1

    unresolved_strict = sorted(
        [
            r["id"]
            for r in frontier
            if r["classification"] in {
                "strict_candidate_target_miss",
                "underdetermined_without_new_observable",
                "model_assumption_required",
                "coupled_underdetermined_pair",
                "anchor_dependent_baseline_only",
            }
        ]
    )

    strict_closure_ready = len(unresolved_in_scope) == 0
    verdict = (
        "MISSING14_STRICT_RIGOR_FRONTIER_PASS_ALL_CLOSED"
        if strict_closure_ready
        else "MISSING14_STRICT_RIGOR_FRONTIER_PARTIAL_ONLY"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "qw2069": "report_qw2069_full_sm_gr_derivation_package.json",
            "qw2070": "report_qw2070_full_radiative_program_baseline.json",
            "qw2072": "report_qw2072_ew_yukawa_flavor_radiative_baselines.json",
            "qw2073": "report_qw2073_radiative_channels_closure_upgrade.json",
            "qw2075": "report_qw2075_strict_cp_phase_derivation_gate.json",
        },
        "expected_missing_14": expected_missing_14,
        "resolved_in_tracked_scope": resolved_in_scope,
        "unresolved_in_tracked_scope": unresolved_in_scope,
        "frontier": frontier,
        "classification_counts": by_class,
        "strict_closure_ready_for_all_14": strict_closure_ready,
        "strict_unresolved_ids": unresolved_strict,
        "verdict": verdict,
        "required_next_step": (
            "EXECUTE_INDEPENDENT_MULTITEAM_CONFIRMATORY_PACKAGE"
            if strict_closure_ready
            else "ADD_NEW_FIRST_PRINCIPLES_OBSERVABLES_FOR_G_EW_HIGGS_LIGHTQUARK_NEUTRINO_COSMO_ABSOLUTE_SCALE"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2081: MISSING-14 STRICT RIGOR FRONTIER",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- strict_closure_ready_for_all_14: `{out['strict_closure_ready_for_all_14']}`",
        "",
        "## Classification Counts",
    ]
    for k in sorted(by_class.keys()):
        lines.append(f"- {k}: {by_class[k]}")

    lines.extend(
        [
            "",
            "## Strict Candidate Snapshot",
            f"- delta_cp_ckm best-branch rel_err_pct: {ckm_rel_err:.6f}",
            f"- delta_cp_ckm tolerance_rel_pct: {ckm_tol:.6f}",
            "",
            "## Strict-Unresolved IDs",
        ]
    )
    for pid in unresolved_strict:
        lines.append(f"- {pid}")

    lines.extend(["", "## Artifact", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2081] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2081] Saved MD:   {OUT_MD.name}")
    print(
        f"[QW-2081] verdict={out['verdict']} "
        f"strict_unresolved={len(unresolved_strict)}/14"
    )


if __name__ == "__main__":
    main()
