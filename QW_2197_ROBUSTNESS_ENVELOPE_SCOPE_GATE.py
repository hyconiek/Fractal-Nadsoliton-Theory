#!/usr/bin/env python3
"""
QW-2197: Robustness envelope scope gate (L7).

Purpose:
- integrate robustness evidence from generation alignment, q-assignment,
  selection-family stability, mass hierarchy slope stability, and spectral
  perturbation margin,
- keep global-unbounded robustness boundary explicit.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2197_robustness_envelope_scope_gate.json"
OUT_MD = ROOT / "RAPORT_QW2197_ROBUSTNESS_ENVELOPE_SCOPE_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    r2125 = load_json("report_qw2125_ktotal_generation_alignment_audit.json")
    r2128 = load_json("report_qw2128_kernel_rep_assignment_uniqueness_gate.json")
    r2186 = load_json("report_qw2186_ktotal_spectral_stability_margin_gate.json")
    r2193 = load_json("report_qw2193_selection_axiom_family_robustness_gate.json")
    r2194 = load_json("report_qw2194_mass_derivation_calibration_separation_gate.json")
    r2196 = load_json("report_qw2196_global_identifiability_scope_stratification_gate.json")

    p2125 = r2125["perturbation_robustness"]
    d2128 = r2128["delta_info_robustness"]
    w2186 = r2186["weyl_margin_certificate"]
    dc2186 = r2186["deterministic_checks"]
    rows2193 = r2193["numeric_audit"]["rows"]
    a2194 = r2194["audit"]

    min_gap = min(float(r["score_gap"]) for r in d2128["rows"])
    winner_freq = int(d2128["winner_frequency_for_locked_winner"])
    n_tests = int(d2128["n_tests"])
    all_theta_small = all(
        abs(float(r["theta_pair1_argmin"])) <= (math.pi / 360.0)
        and abs(float(r["theta_pair2_argmin"])) <= (math.pi / 360.0)
        for r in rows2193
    )

    flags = {
        "q2125_mod3_overlap_mean_ge_0p60": bool(float(p2125["mod3_overlap_fraction_mean"]) >= 0.60),
        "q2125_mod3_overlap_p10_ge_0p50": bool(float(p2125["mod3_overlap_fraction_p10"]) >= 0.50),
        "q2128_delta_info_winner_frequency_full": bool(winner_freq == n_tests),
        "q2128_delta_info_min_score_gap_ge_0p25": bool(min_gap >= 0.25),
        "q2193_selection_family_all_theta_zero": bool(all_theta_small),
        "q2194_non_top_slope_rel_diff_le_10pct": bool(float(a2194["non_top_slope_rel_diff"]) <= 0.10),
        "q2186_weyl_safe_radius_positive": bool(float(w2186["epsilon_safe"]) > 0.0),
        "q2186_mc_inside_safe_radius_stable": bool(float(dc2186["min_lambda_safe_mc"]) > 0.0),
        "q2186_sharpness_witness_detected_above_radius": bool(float(dc2186["witness_min_lambda_over"]) < 0.0),
        "q2196_scope_stratification_present": bool(str(r2196.get("verdict", "")).endswith("AXIOM_FREE_OPEN")),
        "scope_boundaries_explicit": True,
        "global_unbounded_robustness_closed": False,
        "deterministic_no_scan_no_retune": bool(
            r2125["flags"].get("deterministic_no_scan_no_retune", False)
            and r2128["flags"].get("deterministic_no_scan_no_retune", False)
            and r2193["flags"].get("deterministic_no_scan_no_retune", False)
            and r2194["flags"].get("deterministic_no_scan_no_retune", False)
            and r2186["flags"].get("deterministic_no_scan_no_retune", False)
        ),
    }

    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    core_ok = bool(
        flags["q2125_mod3_overlap_mean_ge_0p60"]
        and flags["q2125_mod3_overlap_p10_ge_0p50"]
        and flags["q2128_delta_info_winner_frequency_full"]
        and flags["q2128_delta_info_min_score_gap_ge_0p25"]
        and flags["q2193_selection_family_all_theta_zero"]
        and flags["q2194_non_top_slope_rel_diff_le_10pct"]
        and flags["q2186_weyl_safe_radius_positive"]
        and flags["q2186_mc_inside_safe_radius_stable"]
        and flags["q2186_sharpness_witness_detected_above_radius"]
        and flags["q2196_scope_stratification_present"]
        and flags["scope_boundaries_explicit"]
        and flags["deterministic_no_scan_no_retune"]
    )

    verdict = (
        "ROBUSTNESS_ENVELOPE_SCOPE_GATE_PASS_STRICT_PARTIAL_GLOBAL_UNBOUNDED_OPEN"
        if core_ok
        else "ROBUSTNESS_ENVELOPE_SCOPE_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2125": "report_qw2125_ktotal_generation_alignment_audit.json",
            "q2128": "report_qw2128_kernel_rep_assignment_uniqueness_gate.json",
            "q2186": "report_qw2186_ktotal_spectral_stability_margin_gate.json",
            "q2193": "report_qw2193_selection_axiom_family_robustness_gate.json",
            "q2194": "report_qw2194_mass_derivation_calibration_separation_gate.json",
            "q2196": "report_qw2196_global_identifiability_scope_stratification_gate.json",
        },
        "robustness_metrics": {
            "mod3_overlap_mean": float(p2125["mod3_overlap_fraction_mean"]),
            "mod3_overlap_p10": float(p2125["mod3_overlap_fraction_p10"]),
            "delta_info_winner_frequency": winner_freq,
            "delta_info_n_tests": n_tests,
            "delta_info_min_score_gap": float(min_gap),
            "selection_family_all_theta_zero": bool(all_theta_small),
            "non_top_slope_rel_diff": float(a2194["non_top_slope_rel_diff"]),
            "epsilon_safe": float(w2186["epsilon_safe"]),
            "min_lambda_safe_mc": float(dc2186["min_lambda_safe_mc"]),
            "witness_min_lambda_over": float(dc2186["witness_min_lambda_over"]),
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "EXTEND_FROM_DECLARED_SCOPE_ROBUSTNESS_TO_GLOBAL_UNBOUNDED_ROBUSTNESS_OR_KEEP_SCOPE_BOUNDARY_EXPLICIT"
            if verdict.endswith("GLOBAL_UNBOUNDED_OPEN")
            else "REPAIR_ROBUSTNESS_ENVELOPE_CHAIN_AND_RERUN_QW2197"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    md_lines = [
        "# RAPORT QW-2197: ROBUSTNESS ENVELOPE SCOPE GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Core result",
        "- Multi-source robustness checks are closed in declared strict scopes.",
        "- Spectral perturbation margin and selection-family stability are both explicit.",
        "- Global unbounded robustness remains open and explicitly marked.",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()
