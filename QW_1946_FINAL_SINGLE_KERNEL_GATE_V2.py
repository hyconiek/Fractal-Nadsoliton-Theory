#!/usr/bin/env python3
"""
QW-1946: Final single-kernel gate v2 after Q-audit + kernel-invariant operator.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1946_final_single_kernel_gate_v2.json"
OUT_MD = ROOT / "RAPORT_QW1946_FINAL_SINGLE_KERNEL_GATE_V2.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def rel_miss_hi(value: float, limit: float) -> float:
    return float(max(0.0, value - limit) / max(limit, 1e-12))


def rel_miss_lo(value: float, limit: float) -> float:
    return float(max(0.0, limit - value) / max(limit, 1e-12))


def main() -> None:
    d1934 = load("report_qw1934_strict_closure_gate_solo.json")
    d1941 = load("report_qw1941_triple_sector_shared_complexity_aic_bic.json")
    d1943 = load("report_qw1943_topological_q_assignment_audit.json")
    d1944 = load("report_qw1944_kernel_invariant_flavor_operator_and_gw.json")
    d1945 = load("report_qw1945_physical_gamma_extraction_hard_mass.json")

    stage_b_closed = d1934.get("readiness") == "TOE_STAGE_B_SOLO_CLOSED"
    mass_pass = bool(d1945.get("primary_pass", False))
    strict_flavor_gw_pass = bool(d1944["baseline"]["all_pass"])
    best_flavor_gw_pass = bool(d1944["best_scenario"]["all_pass"])
    has_joint_q_assignment_pass = int(d1943["counts"]["joint_pass_count"]) > 0

    complexity_guard = bool(d1941["model_comparison"]["delta_bic_m1_minus_m0"] <= -6.0)

    strict_all = bool(stage_b_closed and mass_pass and strict_flavor_gw_pass)
    extended_all = bool(
        stage_b_closed
        and mass_pass
        and best_flavor_gw_pass
        and has_joint_q_assignment_pass
        and complexity_guard
    )

    if strict_all:
        readiness = "TOE_STAGE_C_SINGLE_KERNEL_CLOSED_STRICT_V2"
        verdict = "FINAL_SINGLE_KERNEL_GATE_V2_PASS_STRICT"
        required_next = "LOCK_STRICT_V2_AND_EXECUTE_INDEPENDENT_CONFIRMATORY_PACKAGE"
    elif extended_all:
        readiness = "TOE_STAGE_C_SINGLE_KERNEL_CLOSED_EXTENDED_V2"
        verdict = "FINAL_SINGLE_KERNEL_GATE_V2_PASS_EXTENDED"
        required_next = "PUBLISH_EXTENDED_V2_WITH_EXPLICIT_COMPLEXITY_PENALTY"
    elif stage_b_closed:
        readiness = "TOE_STAGE_C_BLOCKED_V2"
        verdict = "FINAL_SINGLE_KERNEL_GATE_V2_FAIL"
        required_next = "REPAIR_KERNEL_TO_MASS_FLAVOR_LINK_OR_REVISE_CORE_HYPOTHESIS"
    else:
        readiness = "TOE_STAGE_B_BLOCKED"
        verdict = "FINAL_SINGLE_KERNEL_GATE_V2_BLOCKED_BY_STAGE_B"
        required_next = "RECOVER_STAGE_B_FIRST"

    th_m = d1945["thresholds"]
    primary = next(v for v in d1945["variants"] if v["label"] == d1945["primary_variant"])
    th_fg = d1944["thresholds"]
    base = d1944["baseline"]

    blockers: List[Dict[str, float]] = [
        {
            "metric": "mass_mean_rel_pct",
            "relative_miss": rel_miss_hi(primary["mass"]["mean_rel_err_pct"], th_m["mean_rel_err_pct_max"]),
        },
        {
            "metric": "mass_max_rel_pct",
            "relative_miss": rel_miss_hi(primary["mass"]["max_rel_err_pct"], th_m["max_rel_err_pct_max"]),
        },
        {
            "metric": "ckm_mean_rel_pct",
            "relative_miss": rel_miss_hi(base["flavor"]["ckm_mean_rel_pct"], th_fg["ckm_mean_rel_pct_max"]),
        },
        {
            "metric": "pmns_mean_rel_pct",
            "relative_miss": rel_miss_hi(base["flavor"]["pmns_mean_rel_pct"], th_fg["pmns_mean_rel_pct_max"]),
        },
        {
            "metric": "gw_auc",
            "relative_miss": rel_miss_lo(base["gw"]["auc_h1l1_vs_ctrl"], th_fg["gw_auc_min"]),
        },
        {
            "metric": "gw_adv",
            "relative_miss": rel_miss_lo(base["gw"]["adv_shared_minus_ctrl_q90"], th_fg["gw_adv_min"]),
        },
        {
            "metric": "gw_sep",
            "relative_miss": rel_miss_lo(base["gw"]["sep_median_h1l1_minus_ctrl"], th_fg["gw_sep_min"]),
        },
        {
            "metric": "gw_control_gap",
            "relative_miss": rel_miss_hi(base["gw"]["control_median_gap"], th_fg["gw_control_gap_max"]),
        },
        {
            "metric": "q_assignment_joint_feasibility",
            "relative_miss": float(0.0 if has_joint_q_assignment_pass else 1.0),
        },
    ]
    blockers = sorted(blockers, key=lambda x: x["relative_miss"], reverse=True)

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "q1934_readiness": d1934.get("readiness"),
            "q1943_verdict": d1943.get("verdict"),
            "q1944_verdict": d1944.get("verdict"),
            "q1945_verdict": d1945.get("verdict"),
            "q1941_delta_bic_m1_minus_m0": d1941["model_comparison"]["delta_bic_m1_minus_m0"],
        },
        "flags": {
            "stage_b_closed": stage_b_closed,
            "mass_pass_q1945_primary": mass_pass,
            "strict_flavor_gw_pass_q1944_baseline": strict_flavor_gw_pass,
            "best_flavor_gw_pass_q1944": best_flavor_gw_pass,
            "q_assignment_joint_pass_exists_q1943": has_joint_q_assignment_pass,
            "complexity_guard_q1941": complexity_guard,
        },
        "blocker_ranking": blockers,
        "readiness": readiness,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1946: FINAL SINGLE-KERNEL GATE V2",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Readiness: **{readiness}**",
        "",
        "## Flags",
        f"- stage_b_closed: {stage_b_closed}",
        f"- mass_pass_q1945_primary: {mass_pass}",
        f"- strict_flavor_gw_pass_q1944_baseline: {strict_flavor_gw_pass}",
        f"- best_flavor_gw_pass_q1944: {best_flavor_gw_pass}",
        f"- q_assignment_joint_pass_exists_q1943: {has_joint_q_assignment_pass}",
        f"- complexity_guard_q1941: {complexity_guard}",
        "",
        "## Top Blockers",
    ]
    for row in blockers[:6]:
        lines.append(f"- {row['metric']}: {row['relative_miss']:.4f}")
    lines.extend(
        [
            "",
            "## Required Next Step",
            f"- {required_next}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1946] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1946] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1946] readiness={readiness} verdict={verdict}")


if __name__ == "__main__":
    main()

