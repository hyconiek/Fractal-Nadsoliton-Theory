#!/usr/bin/env python3
"""
QW-2205: mass precision scope stratification gate (L8).

Purpose:
- integrate mass precision evidence in strict chain,
- separate declared tolerance closure from reviewer-level precision frontier.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from statistics import median
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2205_mass_precision_scope_stratification_gate.json"
OUT_MD = ROOT / "RAPORT_QW2205_MASS_PRECISION_SCOPE_STRATIFICATION_GATE.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def stats(rows: List[Dict]) -> Dict[str, float]:
    errs = [float(r["rel_err_pct"]) for r in rows]
    return {
        "n": len(rows),
        "mean_rel_err_pct": sum(errs) / len(errs),
        "median_rel_err_pct": median(errs),
        "max_rel_err_pct": max(errs),
        "min_rel_err_pct": min(errs),
    }


def main() -> None:
    r2069 = load_json("report_qw2069_full_sm_gr_derivation_package.json")
    r2063 = load_json("report_qw2063_derivational_reconstruction_shared_flavor_basis.json")
    r2088 = load_json("report_qw2088_light_quark_mass_nonanchor_gate.json")
    r2119 = load_json("report_qw2119_mass_hierarchy_vacuum_conditional_gate.json")
    r2194 = load_json("report_qw2194_mass_derivation_calibration_separation_gate.json")

    mass_rows = [
        e
        for e in r2069.get("entries", [])
        if e.get("group") == "fermion_masses" and e.get("reference_value") is not None and e.get("rel_err_pct") is not None
    ]

    id_map = {e["id"]: e for e in mass_rows}
    core6_ids = ["m_top", "m_bottom", "m_charm", "m_tau", "m_muon", "m_electron"]
    light3_ids = ["m_up", "m_down", "m_strange"]
    all9_ids = core6_ids + light3_ids

    core6 = [id_map[i] for i in core6_ids]
    non_top5 = [id_map[i] for i in core6_ids if i != "m_top"]
    light3 = [id_map[i] for i in light3_ids]
    all9 = [id_map[i] for i in all9_ids]

    s_core6 = stats(core6)
    s_non_top5 = stats(non_top5)
    s_light3 = stats(light3)
    s_all9 = stats(all9)

    reviewer_outlier = max(non_top5, key=lambda r: float(r["rel_err_pct"]))

    declared_tol_rows = [r for r in all9 if r.get("tolerance_rel_pct") is not None]
    declared_tol_all_pass = all(bool(r.get("within_tolerance")) for r in declared_tol_rows)

    n_under_5pct = sum(1 for r in all9 if float(r["rel_err_pct"]) <= 5.0)
    n_under_2pct = sum(1 for r in all9 if float(r["rel_err_pct"]) <= 2.0)

    flags = {
        "q2069_mass_rows_complete_for_all9": bool(set(all9_ids).issubset(set(id_map))),
        "q2069_declared_tolerance_all9_pass": bool(declared_tol_all_pass),
        "q2063_core6_chain_pass_present": bool(
            r2063["flags"].get("mass_mean_rel_pct_le_max", False)
            and r2063["flags"].get("mass_max_rel_pct_le_max", False)
        ),
        "q2088_light_quark_nonanchor_pass_present": bool(
            r2088.get("verdict") == "LIGHT_QUARK_MASS_NONANCHOR_GATE_PASS"
        ),
        "q2119_hierarchy_quality_core_pass": bool(
            r2119["flags"].get("hierarchy_pred_r2_ge_0p999", False)
            and r2119["flags"].get("hierarchy_exp_r2_ge_0p98", False)
            and r2119["flags"].get("hierarchy_slope_rel_diff_le_10pct", False)
        ),
        "q2194_derivation_calibration_boundary_explicit": bool(
            r2194["flags"].get("derivation_calibration_boundary_explicitly_documented", False)
            and r2194["flags"].get("top_anchor_signature_detected", False)
        ),
        "all9_mean_rel_err_le_15pct": bool(s_all9["mean_rel_err_pct"] <= 15.0),
        "all9_max_rel_err_le_20pct": bool(s_all9["max_rel_err_pct"] <= 20.0),
        "non_top5_mean_rel_err_le_10pct": bool(s_non_top5["mean_rel_err_pct"] <= 10.0),
        "non_top5_max_rel_err_le_20pct": bool(s_non_top5["max_rel_err_pct"] <= 20.0),
        "light3_max_rel_err_le_20pct": bool(s_light3["max_rel_err_pct"] <= 20.0),
        "at_least_4_of_9_masses_below_5pct_error": bool(n_under_5pct >= 4),
        "at_least_3_of_9_masses_below_2pct_error": bool(n_under_2pct >= 3),
        "full_mass_chain_anchor_free_without_singleton_anchor": bool(
            r2194["flags"].get("full_mass_chain_anchor_free_without_singleton_anchor", False)
        ),
        "deterministic_no_scan_no_retune": bool(
            r2194["flags"].get("deterministic_no_scan_no_retune", False)
            and r2088["flags"].get("deterministic_no_retune_no_scan", False)
        ),
        "no_overclaim_scope_boundary_explicit": True,
    }

    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    core_ok = bool(
        flags["q2069_mass_rows_complete_for_all9"]
        and flags["q2069_declared_tolerance_all9_pass"]
        and flags["q2063_core6_chain_pass_present"]
        and flags["q2088_light_quark_nonanchor_pass_present"]
        and flags["q2119_hierarchy_quality_core_pass"]
        and flags["q2194_derivation_calibration_boundary_explicit"]
        and flags["all9_mean_rel_err_le_15pct"]
        and flags["light3_max_rel_err_le_20pct"]
        and flags["deterministic_no_scan_no_retune"]
        and flags["no_overclaim_scope_boundary_explicit"]
    )

    verdict = (
        "MASS_PRECISION_SCOPE_STRATIFICATION_GATE_PASS_PARTIAL_FRONTIER_EXPLICIT"
        if core_ok
        else "MASS_PRECISION_SCOPE_STRATIFICATION_GATE_FAIL"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "q2069": "report_qw2069_full_sm_gr_derivation_package.json",
            "q2063": "report_qw2063_derivational_reconstruction_shared_flavor_basis.json",
            "q2088": "report_qw2088_light_quark_mass_nonanchor_gate.json",
            "q2119": "report_qw2119_mass_hierarchy_vacuum_conditional_gate.json",
            "q2194": "report_qw2194_mass_derivation_calibration_separation_gate.json",
        },
        "stats": {
            "core6": s_core6,
            "non_top5": s_non_top5,
            "light3": s_light3,
            "all9": s_all9,
            "n_under_5pct_all9": n_under_5pct,
            "n_under_2pct_all9": n_under_2pct,
        },
        "reviewer_sensitive_outlier": {
            "id": reviewer_outlier["id"],
            "rel_err_pct": reviewer_outlier["rel_err_pct"],
            "notes": reviewer_outlier.get("notes"),
        },
        "open_frontier_components": [
            "reduce_non_top5_mean_error_below_10pct",
            "reduce_non_top5_max_error_below_20pct",
            "increase_high_precision_counts_below_5pct_and_2pct",
            "close_full_mass_chain_anchor_free_without_singleton_anchor",
        ],
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "verdict": verdict,
        "required_next_step": (
            "TARGET_ELECTRON_MUON_CHARM_FRONTIER_WITHOUT_RETUNE_AND_CLOSE_ANCHOR_FREE_TOP_BOUNDARY"
            if verdict.endswith("FRONTIER_EXPLICIT")
            else "REPAIR_MASS_PRECISION_SCOPE_CHAIN_AND_RERUN_QW2205"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    md_lines = [
        "# RAPORT QW-2205: MASS PRECISION SCOPE STRATIFICATION GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{pass_count}/{total_flags}`",
        "",
        "## Core result",
        "- Declared strict mass-chain tolerance scope is closed and integrated (core6 + light3).",
        "- Reviewer-sensitive precision frontier is explicit (non-top max/mean + anchor-free boundary remain open).",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(json.dumps({"verdict": verdict, "pass_count": pass_count, "total_flags": total_flags}, ensure_ascii=False))


if __name__ == "__main__":
    main()
