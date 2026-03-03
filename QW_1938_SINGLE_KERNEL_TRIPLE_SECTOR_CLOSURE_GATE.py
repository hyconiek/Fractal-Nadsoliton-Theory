#!/usr/bin/env python3
"""
QW-1938: Single-kernel triple-sector closure gate.

Integrates:
- prior solo closure stage status (QW-1934),
- new unified mass+flavor+GW no-sector-retune run (QW-1937),
and declares whether single-frozen-kernel closure is achieved.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1938_single_kernel_triple_sector_closure_gate.json"
OUT_MD = ROOT / "RAPORT_QW1938_SINGLE_KERNEL_TRIPLE_SECTOR_CLOSURE_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def blocker_ranking(best: Dict, thresholds: Dict[str, float]) -> List[Tuple[str, float]]:
    m = best["mass"]
    f = best["flavor"]
    g = best["gw"]["summary"]

    misses = {
        "mass_mean_rel_pct": max(0.0, m["mean_rel_err_pct"] - thresholds["mass_mean_rel_pct_max"])
        / max(thresholds["mass_mean_rel_pct_max"], 1e-12),
        "mass_max_rel_pct": max(0.0, m["max_rel_err_pct"] - thresholds["mass_max_rel_pct_max"])
        / max(thresholds["mass_max_rel_pct_max"], 1e-12),
        "ckm_mean_rel_pct": max(0.0, f["ckm_mean_rel_pct"] - thresholds["ckm_mean_rel_pct_max"])
        / max(thresholds["ckm_mean_rel_pct_max"], 1e-12),
        "pmns_mean_rel_pct": max(0.0, f["pmns_mean_rel_pct"] - thresholds["pmns_mean_rel_pct_max"])
        / max(thresholds["pmns_mean_rel_pct_max"], 1e-12),
        "gw_auc": max(0.0, thresholds["gw_auc_min"] - g["auc_h1l1_vs_ctrl"]) / max(thresholds["gw_auc_min"], 1e-12),
        "gw_adv": max(0.0, thresholds["gw_adv_min"] - g["adv_shared_minus_ctrl_q90"])
        / max(thresholds["gw_adv_min"], 1e-12),
        "gw_sep": max(0.0, thresholds["gw_sep_min"] - g["sep_median_h1l1_minus_ctrl"])
        / max(thresholds["gw_sep_min"], 1e-12),
        "gw_control_gap": max(0.0, g["control_median_gap"] - thresholds["gw_control_gap_max"])
        / max(thresholds["gw_control_gap_max"], 1e-12),
    }
    return sorted(misses.items(), key=lambda kv: kv[1], reverse=True)


def main() -> None:
    d1934 = load("report_qw1934_strict_closure_gate_solo.json")
    d1937 = load("report_qw1937_unified_frozen_kernel_mass_flavor_gw_no_sector_retune.json")

    solo_ready = d1934.get("readiness") == "TOE_STAGE_B_SOLO_CLOSED"
    derived_all_pass = bool(d1937["derived_candidate"]["all_pass"])
    global_all_pass = bool(d1937["joint_single_fit_best"]["all_pass"])
    thresholds = d1937["thresholds"]
    best = d1937["joint_single_fit_best"]

    if solo_ready and derived_all_pass:
        readiness = "TOE_STAGE_C_SINGLE_KERNEL_CLOSED_DERIVED"
        verdict = "SINGLE_KERNEL_TRIPLE_SECTOR_CLOSURE_PASS_DERIVED"
        required_next = "FREEZE_PREDICTIONS_AND_EXECUTE_INDEPENDENT_CONFIRMATORY_PACKAGE"
    elif solo_ready and global_all_pass:
        readiness = "TOE_STAGE_C_SINGLE_KERNEL_CLOSED_GLOBAL_SHARED_FIT"
        verdict = "SINGLE_KERNEL_TRIPLE_SECTOR_CLOSURE_PASS_GLOBAL_SHARED_FIT"
        required_next = "PUBLISH_FROZEN_GLOBAL_VECTOR_AND_RUN_TRUE_EXTERNAL_REPLICATION"
    elif solo_ready:
        readiness = "TOE_STAGE_C_BLOCKED_SINGLE_KERNEL_NOT_PASSING_TRIPLE_SECTOR"
        verdict = "SINGLE_KERNEL_TRIPLE_SECTOR_CLOSURE_FAIL"
        required_next = "REPAIR_SHARED_FLAVOR_MECHANISM_WITHOUT_SECTOR_RETUNING"
    else:
        readiness = "TOE_STAGE_B_NOT_CLOSED_CANNOT_EVALUATE_STAGE_C"
        verdict = "SINGLE_KERNEL_TRIPLE_SECTOR_CLOSURE_BLOCKED_BY_STAGE_B"
        required_next = "RECOVER_STAGE_B_SOLO_CLOSURE_FIRST"

    ranking = blocker_ranking(best=best, thresholds=thresholds)

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "q1934_readiness": d1934.get("readiness"),
            "q1934_verdict": d1934.get("verdict"),
            "q1937_verdict": d1937.get("verdict"),
            "q1937_feasible_count": d1937["protocol"]["feasible_count_all_pass"],
        },
        "flags": {
            "stage_b_solo_closed": bool(solo_ready),
            "q1937_derived_all_pass": bool(derived_all_pass),
            "q1937_global_shared_all_pass": bool(global_all_pass),
        },
        "blocker_ranking": [{"metric": k, "relative_miss": float(v)} for k, v in ranking],
        "readiness": readiness,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1938: SINGLE KERNEL TRIPLE-SECTOR CLOSURE GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Readiness: **{readiness}**",
        "",
        "## Inputs",
        f"- QW-1934 readiness/verdict: {d1934.get('readiness')} / {d1934.get('verdict')}",
        f"- QW-1937 verdict: {d1937.get('verdict')}",
        f"- QW-1937 feasible all-pass count: {d1937['protocol']['feasible_count_all_pass']}",
        "",
        "## Flags",
        f"- stage_b_solo_closed: {solo_ready}",
        f"- q1937_derived_all_pass: {derived_all_pass}",
        f"- q1937_global_shared_all_pass: {global_all_pass}",
        "",
        "## Top Blockers (relative miss)",
    ]
    for row in out["blocker_ranking"][:5]:
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

    print(f"[QW-1938] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1938] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1938] readiness={readiness} verdict={verdict}")


if __name__ == "__main__":
    main()

