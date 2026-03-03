#!/usr/bin/env python3
"""
QW-1914: TOE potential and mathematical-coherence scorecard.

Aggregates existing audit/gate artifacts into a compact status view:
- empirical status,
- robustness against split/multisplit selection,
- derivational closure depth,
- claim hygiene and reproducibility diagnostics,
- qualitative comparison matrix vs major theory programs.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1914_toe_potential_scorecard.json"
OUT_MD = ROOT / "RAPORT_QW1914_TOE_POTENTIAL_SCORECARD.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def clip01(x: float) -> float:
    return float(max(0.0, min(1.0, x)))


def main() -> None:
    d1902 = load("report_qw1902_empirical_closure_gate.json")
    d1912 = load("report_qw1912_external_pta_split_validation.json")
    d1913 = load("report_qw1913_external_pta_multisplit_transfer_stress.json")
    d1745 = load("report_qw1745_micromodel_rigor_gate.json")
    d1890 = load("report_qw1890_toe_closure_decision_gate.json")
    d1702 = load("report_qw1702_reproducibility_audit.json")
    d1703 = load("report_qw1703_claims_vs_computation_audit.json")
    d1907 = load("report_qw1907_pre1700_tuning_boundary_audit.json")

    # Pillar 1: empirical closure
    readiness_1902 = str(d1902.get("readiness", ""))
    metric_score = float(d1902.get("metric_score", 0.0))
    externality_ok = bool(d1902.get("inputs", {}).get("externality_ok", False))
    empirical_score = 0.70 * clip01(metric_score) + 0.30 * float(
        readiness_1902 == "EMPIRICAL_CLOSURE_PASS" and externality_ok
    )

    # Pillar 2: robustness (split + multisplit transfer)
    holdout_single_ok = bool((d1912.get("holdout_validation") or {}).get("all_pass", False))
    ms = d1913.get("summary", {})
    holdout_pass_rate = float(ms.get("holdout_pass_rate_given_selected", 0.0))
    selected_alphas = ms.get("selected_alphas", []) or []
    alpha_stability = 1.0 if len(selected_alphas) >= 2 and float(np.std(selected_alphas)) < 1e-9 else (
        0.8 if len(selected_alphas) >= 2 and float(np.std(selected_alphas)) <= 0.5 else 0.5
    )
    robustness_score = (
        0.40 * float(holdout_single_ok)
        + 0.45 * clip01(holdout_pass_rate)
        + 0.15 * float(alpha_stability)
    )

    # Pillar 3: derivational closure depth
    deriv_1745 = float(d1745.get("global_score", 0.0))
    deriv_1890 = float(d1890.get("global_score", 0.0))
    derivational_score = clip01(0.5 * deriv_1745 + 0.5 * deriv_1890)

    # Pillar 4: claim hygiene
    mass_err = []
    for row in d1703.get("mass_formula_check", []):
        try:
            mass_err.append(abs(float(row.get("error_pct", 0.0))))
        except Exception:
            pass
    mean_mass_err = float(np.mean(mass_err)) if mass_err else 100.0
    claim_obs_score = clip01(1.0 - mean_mass_err / 25.0)

    dd = d1703.get("declaration_density", {})
    exact_mentions = float(dd.get("exact_or_0_00_mentions_in_core_docs", 0.0))
    fit_mentions = float(dd.get("fitting_or_calibration_mentions_in_core_docs", 0.0))
    mention_ratio = exact_mentions / max(1.0, fit_mentions)
    claim_text_score = clip01(1.0 - max(0.0, mention_ratio - 1.0) / 3.0)
    claim_hygiene_score = 0.70 * claim_obs_score + 0.30 * claim_text_score

    # Pillar 5: reproducibility
    tasks = d1702.get("tasks", [])
    if tasks:
        semantic_det_rate = float(
            np.mean([1.0 if bool(t.get("deterministic_outputs_semantic", False)) else 0.0 for t in tasks])
        )
    else:
        semantic_det_rate = 0.0
    reproducibility_score = semantic_det_rate

    # Composite potential
    overall_score = (
        0.32 * empirical_score
        + 0.23 * robustness_score
        + 0.22 * derivational_score
        + 0.13 * claim_hygiene_score
        + 0.10 * reproducibility_score
    )

    if overall_score >= 0.80:
        potential = "HIGH"
    elif overall_score >= 0.65:
        potential = "MODERATE_HIGH"
    elif overall_score >= 0.50:
        potential = "MODERATE"
    else:
        potential = "LOW"

    if empirical_score >= 0.9 and derivational_score < 0.6:
        closure_profile = "EMPIRICAL_STRONG_DERIVATIONAL_OPEN"
    elif empirical_score >= 0.9 and derivational_score >= 0.6:
        closure_profile = "EMPIRICAL_AND_DERIVATIONAL_STRONG"
    elif empirical_score < 0.9 and derivational_score >= 0.6:
        closure_profile = "DERIVATIONAL_AHEAD_OF_EMPIRICAL"
    else:
        closure_profile = "BOTH_OPEN_OR_PARTIAL"

    comparison = [
        {
            "program": "SM+GR (effective benchmark)",
            "empirical_fit_today": "very_high",
            "unification_depth": "low",
            "free_parameter_load": "high",
            "near_term_falsifiability": "very_high",
        },
        {
            "program": "String/M-theory (broad class)",
            "empirical_fit_today": "low_direct",
            "unification_depth": "high_formal",
            "free_parameter_load": "high_landscape",
            "near_term_falsifiability": "low",
        },
        {
            "program": "Loop Quantum Gravity",
            "empirical_fit_today": "low_direct",
            "unification_depth": "medium_gravity_focused",
            "free_parameter_load": "medium",
            "near_term_falsifiability": "low_to_medium",
        },
        {
            "program": "Asymptotic Safety (gravity-first)",
            "empirical_fit_today": "low_direct",
            "unification_depth": "medium",
            "free_parameter_load": "medium",
            "near_term_falsifiability": "medium",
        },
        {
            "program": "Nadsoliton TOE (current repo state)",
            "empirical_fit_today": "medium_to_high_internal_external_archive",
            "unification_depth": "medium_partial_derivational_gap",
            "free_parameter_load": "medium_alpha_bridge_active",
            "near_term_falsifiability": "high_with_new_blind_external_run",
        },
    ]

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "q1902_readiness": readiness_1902,
            "q1912_verdict": d1912.get("verdict"),
            "q1913_verdict": d1913.get("verdict"),
            "q1745_readiness": d1745.get("readiness"),
            "q1890_readiness": d1890.get("readiness"),
            "q1907_overall": d1907.get("verdict", {}).get("overall"),
        },
        "scores": {
            "empirical_closure": empirical_score,
            "robustness_transfer": robustness_score,
            "derivational_closure": derivational_score,
            "claim_hygiene": claim_hygiene_score,
            "reproducibility_semantic": reproducibility_score,
            "overall_potential": overall_score,
        },
        "diagnostics": {
            "metric_score_1902": metric_score,
            "externality_ok_1902": externality_ok,
            "single_holdout_pass_1912": holdout_single_ok,
            "multisplit_holdout_pass_rate_1913": holdout_pass_rate,
            "selected_alphas_1913": selected_alphas,
            "mean_mass_error_pct_q1703": mean_mass_err,
            "semantic_determinism_rate_q1702": semantic_det_rate,
        },
        "classification": {
            "toe_potential": potential,
            "closure_profile": closure_profile,
        },
        "qualitative_comparison_matrix": comparison,
        "verdict": "TOE_POTENTIAL_SCORECARD_COMPLETE",
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1914: TOE POTENTIAL SCORECARD",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- toe_potential: **{potential}**",
        f"- closure_profile: **{closure_profile}**",
        "",
        "## Scores (0..1)",
        f"- empirical_closure: {empirical_score:.3f}",
        f"- robustness_transfer: {robustness_score:.3f}",
        f"- derivational_closure: {derivational_score:.3f}",
        f"- claim_hygiene: {claim_hygiene_score:.3f}",
        f"- reproducibility_semantic: {reproducibility_score:.3f}",
        f"- overall_potential: {overall_score:.3f}",
        "",
        "## Key Diagnostics",
        f"- QW-1902 readiness: {readiness_1902}",
        f"- QW-1912 holdout pass: {holdout_single_ok}",
        f"- QW-1913 holdout pass rate: {holdout_pass_rate:.3f}",
        f"- selected_alphas (QW-1913): {selected_alphas}",
        f"- mean_mass_error_pct (QW-1703): {mean_mass_err:.3f}",
        "",
        "## Qualitative Position",
        "- Nadsoliton branch currently: empirically strong in current frozen pipeline,",
        "- but derivational closure remains partial and should be strengthened before final TOE claim.",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1914] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1914] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1914] toe_potential={potential} overall={overall_score:.3f}")


if __name__ == "__main__":
    main()

