#!/usr/bin/env python3
"""
QW-1942: Final single-kernel ToE gate (strict branch + shared-lambda branch).

Strict branch:
- stage-B solo closure already pass,
- hard mass formula pass (QW-1939),
- strict kernel-derived flavor pass (QW-1940),
- triple-sector M0 strict-derived pass (QW-1941).

Extended branch:
- one shared lambda branch pass (QW-1941),
- complexity justified by BIC.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1942_final_single_kernel_toe_gate.json"
OUT_MD = ROOT / "RAPORT_QW1942_FINAL_SINGLE_KERNEL_TOE_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def rel_miss_hi(value: float, limit: float) -> float:
    return float(max(0.0, value - limit) / max(limit, 1e-12))


def rel_miss_lo(value: float, limit: float) -> float:
    return float(max(0.0, limit - value) / max(limit, 1e-12))


def compute_blockers_from_m0(d1941: Dict) -> List[Dict[str, float]]:
    th = d1941["thresholds"]
    m0 = d1941["models"]["M0_strict_derived"]
    mass = m0["mass"]
    fl = m0["flavor"]
    gw = m0["gw"]

    rows = [
        {"metric": "mass_mean_rel_pct", "relative_miss": rel_miss_hi(mass["mean_rel_err_pct"], th["mass_mean_rel_pct_max"])},
        {"metric": "mass_max_rel_pct", "relative_miss": rel_miss_hi(mass["max_rel_err_pct"], th["mass_max_rel_pct_max"])},
        {"metric": "ckm_mean_rel_pct", "relative_miss": rel_miss_hi(fl["ckm_mean_rel_pct"], th["ckm_mean_rel_pct_max"])},
        {"metric": "pmns_mean_rel_pct", "relative_miss": rel_miss_hi(fl["pmns_mean_rel_pct"], th["pmns_mean_rel_pct_max"])},
        {"metric": "gw_auc", "relative_miss": rel_miss_lo(gw["auc_h1l1_vs_ctrl"], th["gw_auc_min"])},
        {"metric": "gw_adv", "relative_miss": rel_miss_lo(gw["adv_shared_minus_ctrl_q90"], th["gw_adv_min"])},
        {"metric": "gw_sep", "relative_miss": rel_miss_lo(gw["sep_median_h1l1_minus_ctrl"], th["gw_sep_min"])},
        {"metric": "gw_control_gap", "relative_miss": rel_miss_hi(gw["control_median_gap"], th["gw_control_gap_max"])},
    ]
    return sorted(rows, key=lambda x: x["relative_miss"], reverse=True)


def main() -> None:
    d1934 = load("report_qw1934_strict_closure_gate_solo.json")
    d1939 = load("report_qw1939_hard_mass_formula_baseline.json")
    d1940 = load("report_qw1940_kernel_derived_flavor_operator_strict.json")
    d1941 = load("report_qw1941_triple_sector_shared_complexity_aic_bic.json")

    stage_b_closed = d1934.get("readiness") == "TOE_STAGE_B_SOLO_CLOSED"
    hard_mass_pass = bool(d1939.get("hard_pass", False))
    flavor_derived_pass = bool(d1940.get("all_pass", False))

    m0 = d1941["models"]["M0_strict_derived"]
    m1 = d1941["models"]["M1_one_shared_lambda_best"]
    cmp = d1941["model_comparison"]

    m0_pass = bool(m0["all_pass"])
    m1_pass = bool(m1["all_pass"])
    m1_bic_justified = bool(cmp["delta_bic_m1_minus_m0"] <= -6.0)

    strict_all = bool(stage_b_closed and hard_mass_pass and flavor_derived_pass and m0_pass)
    extended_all = bool(stage_b_closed and m1_pass and m1_bic_justified)

    if strict_all:
        readiness = "TOE_STAGE_C_SINGLE_KERNEL_STRICT_CLOSED"
        verdict = "FINAL_SINGLE_KERNEL_TOE_GATE_PASS_STRICT"
        required_next = "FREEZE_STRICT_MODEL_AND_RUN_INDEPENDENT_CONFIRMATORY_EXECUTION"
    elif extended_all:
        readiness = "TOE_STAGE_C_SINGLE_KERNEL_CLOSED_WITH_ONE_SHARED_EXTENSION"
        verdict = "FINAL_SINGLE_KERNEL_TOE_GATE_PASS_EXTENDED"
        required_next = "PUBLISH_EXTENSION_WITH_COMPLEXITY_AUDIT_AND_RUN_CONFIRMATORY"
    elif stage_b_closed:
        readiness = "TOE_STAGE_C_BLOCKED"
        verdict = "FINAL_SINGLE_KERNEL_TOE_GATE_FAIL"
        required_next = "REPAIR_MASS_FLAVOR_LINK_UNDER_EXACT_SINGLE_KERNEL_CONSTRAINT"
    else:
        readiness = "TOE_STAGE_B_BLOCKED"
        verdict = "FINAL_SINGLE_KERNEL_TOE_GATE_BLOCKED_BY_STAGE_B"
        required_next = "RECOVER_STAGE_B_FIRST"

    blockers = compute_blockers_from_m0(d1941)

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "q1934_readiness": d1934.get("readiness"),
            "q1939_verdict": d1939.get("verdict"),
            "q1940_verdict": d1940.get("verdict"),
            "q1941_verdict": d1941.get("verdict"),
        },
        "flags": {
            "stage_b_closed": stage_b_closed,
            "hard_mass_pass": hard_mass_pass,
            "flavor_derived_pass": flavor_derived_pass,
            "m0_pass": m0_pass,
            "m1_pass": m1_pass,
            "m1_bic_justified": m1_bic_justified,
        },
        "blocker_ranking_m0": blockers,
        "readiness": readiness,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1942: FINAL SINGLE-KERNEL TOE GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Readiness: **{readiness}**",
        "",
        "## Flags",
        f"- stage_b_closed: {stage_b_closed}",
        f"- hard_mass_pass: {hard_mass_pass}",
        f"- flavor_derived_pass: {flavor_derived_pass}",
        f"- m0_pass: {m0_pass}",
        f"- m1_pass: {m1_pass}",
        f"- m1_bic_justified: {m1_bic_justified}",
        "",
        "## Top Blockers (M0 strict-derived)",
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

    print(f"[QW-1942] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1942] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1942] readiness={readiness} verdict={verdict}")


if __name__ == "__main__":
    main()

