#!/usr/bin/env python3
"""P1482 S4.32: strict-only next honest step toward F => L_SM + L_GR."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

SP1_LOCAL = GEN / "p1468_s418_qw2191_sp1_local_pilot_summary.json"
POLICY_V2 = GEN / "p1478_s428_qw2191_sp1_operating_policy_v2.json"
STRESS = GEN / "p1481_s431_qw2191_sp1_cross_scenario_policy_stress_summary.json"

SUMMARY = GEN / "p1482_s432_qw2191_strict_f_to_lsm_lgr_step_summary.json"
OBSTRUCTION = GEN / "p1482_s432_qw2191_strict_f_to_lsm_lgr_step_obstruction.json"


def main() -> None:
    sp1 = json.loads(SP1_LOCAL.read_text(encoding="utf-8"))
    policy = json.loads(POLICY_V2.read_text(encoding="utf-8"))
    stress = json.loads(STRESS.read_text(encoding="utf-8"))

    delta_local = float(sp1["arm_B_with_SP1_metric"]) - float(sp1["arm_A_no_selector_premise_metric"])
    stress_ok = stress["status"] == "PASS_SP1_CROSS_SCENARIO_LOCAL_ONLY"
    probe_shift = float(policy["worst_case_edge_shift"]) + float(policy["safety_margin"])
    in_policy_window = float(policy["safe_shift_min"]) <= probe_shift <= float(policy["safe_shift_max"])

    checks = {
        "strict_only": True,
        "no_legacy_bridge": True,
        "qw2191_closed": False,
        "sp1_local_delta_positive": delta_local > 0.0,
        "sp1_cross_scenario_pass": stress_ok,
        "recommended_shift_in_policy_window": in_policy_window,
    }

    readiness_checks = {k: v for k, v in checks.items() if k != "qw2191_closed"}
    all_green = all(readiness_checks.values())
    status = "PASS_STRICT_NEXT_STEP_READY" if all_green else "FAIL_STRICT_NEXT_STEP_BLOCKED"

    next_step = {
        "id": "S4.33",
        "name": "strict_sector_split_probe",
        "objective": "Construct explicit split candidate F => L_SM + L_GR with bounded mixed term under SP1 policy gate.",
        "acceptance": [
            "No legacy-bridge dependency",
            "No strict-core closure claim",
            "Explicit sector separation witness for SM-like and GR-like channels",
            "Selector obstruction QW-2191 stays explicit unless new selector premise is added",
        ],
    }

    summary = {
        "packet": "P1482",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "route": "strict_only_F_to_LSM_plus_LGR",
        "checks": checks,
        "observables": {
            "delta_local_sp1_minus_baseline": delta_local,
            "probe_shift": probe_shift,
            "safe_shift_min": policy["safe_shift_min"],
            "safe_shift_max": policy["safe_shift_max"],
        },
        "next_honest_step": next_step,
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if all_green:
        if OBSTRUCTION.exists():
            OBSTRUCTION.unlink()
    else:
        obstruction = {
            "packet": "P1482",
            "status": status,
            "rule": "strict next step readiness checks failed",
            "checks": checks,
        }
        OBSTRUCTION.write_text(json.dumps(obstruction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print(f"[P1482] status={status}")


if __name__ == "__main__":
    main()
