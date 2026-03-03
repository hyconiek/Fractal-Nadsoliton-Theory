#!/usr/bin/env python3
"""
QW-1957: Strict Stage-C+Info gate v2 (after QW-1955/1956 branch).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1957_strict_stage_c_plus_info_gate_v2.json"
OUT_MD = ROOT / "RAPORT_QW1957_STRICT_STAGE_C_PLUS_INFO_GATE_V2.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def rel_miss_hi(value: float, limit: float) -> float:
    return float(max(0.0, value - limit) / max(limit, 1e-12))


def rel_miss_lo(value: float, limit: float) -> float:
    return float(max(0.0, limit - value) / max(limit, 1e-12))


def main() -> None:
    d1946 = load("report_qw1946_final_single_kernel_gate_v2.json")
    d1955 = load("report_qw1955_nogo_and_minimal_operator_repair.json")
    d1956 = load("report_qw1956_two_state_observer_with_repaired_operator.json")

    stage_b_closed = bool(d1946["flags"]["stage_b_closed"])
    triad_strict_pass = bool(
        d1946["flags"]["stage_b_closed"]
        and d1946["flags"]["mass_pass_q1945_primary"]
        and d1946["flags"]["strict_flavor_gw_pass_q1944_baseline"]
    )
    repair_pass = bool(d1955.get("repair_pass", False))
    two_state_repair_pass = bool(d1956.get("pass", False))

    strict_all = bool(triad_strict_pass and repair_pass and two_state_repair_pass)

    if strict_all:
        readiness = "TOE_STAGE_C_PLUS_INFO_CLOSED_STRICT_V2"
        verdict = "STRICT_STAGE_C_PLUS_INFO_GATE_V2_PASS"
        required_next = "LOCK_V2_STACK_AND_EXECUTE_EXTERNAL_CONFIRMATORY_PACKAGE"
    elif stage_b_closed:
        readiness = "TOE_STAGE_C_PLUS_INFO_BLOCKED_V2"
        verdict = "STRICT_STAGE_C_PLUS_INFO_GATE_V2_FAIL"
        required_next = "REWORK_MINIMAL_REPAIR_OPERATOR_AND_MASS_FLAVOR_LINK"
    else:
        readiness = "TOE_STAGE_B_BLOCKED"
        verdict = "STRICT_STAGE_C_PLUS_INFO_GATE_V2_BLOCKED_BY_STAGE_B"
        required_next = "RECOVER_STAGE_B_FIRST"

    blockers: List[Dict[str, float]] = []
    for row in d1946.get("blocker_ranking", []):
        blockers.append({"metric": f"triad::{row['metric']}", "relative_miss": float(row["relative_miss"])})

    t55 = d1955["thresholds"]
    m55 = d1955["metrics"]["dual"]
    blockers.extend(
        [
            {"metric": "info1955::dual_accuracy", "relative_miss": rel_miss_lo(float(m55["accuracy"]), float(t55["dual_accuracy_min"]))},
            {
                "metric": "info1955::acc_gain_vs_control",
                "relative_miss": rel_miss_lo(float(d1955["metrics"]["acc_gain_dual_vs_control"]), float(t55["acc_gain_vs_control_min"])),
            },
            {
                "metric": "info1955::info_gain_vs_control",
                "relative_miss": rel_miss_lo(float(d1955["metrics"]["info_gain_dual_vs_control"]), float(t55["info_gain_vs_control_min"])),
            },
            {
                "metric": "info1955::channel_complementarity",
                "relative_miss": rel_miss_lo(float(m55["channel_complementarity"]), float(t55["channel_complementarity_min"])),
            },
        ]
    )

    t56 = d1956["thresholds"]
    m56 = d1956["metrics"]["closed"]
    blockers.extend(
        [
            {"metric": "info1956::closed_accuracy", "relative_miss": rel_miss_lo(float(m56["accuracy"]), float(t56["closed_accuracy_min"]))},
            {
                "metric": "info1956::closed_acc_gain_vs_open",
                "relative_miss": rel_miss_lo(float(d1956["metrics"]["closed_acc_gain_vs_open"]), float(t56["closed_acc_gain_vs_open_min"])),
            },
            {
                "metric": "info1956::closed_info_gain_vs_control",
                "relative_miss": rel_miss_lo(float(d1956["metrics"]["closed_info_gain_vs_control"]), float(t56["closed_info_gain_vs_control_min"])),
            },
            {
                "metric": "info1956::channel_complementarity",
                "relative_miss": rel_miss_lo(float(m56["channel_complementarity"]), float(t56["channel_complementarity_min"])),
            },
        ]
    )
    blockers = sorted(blockers, key=lambda x: x["relative_miss"], reverse=True)

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "q1946_readiness": d1946["readiness"],
            "q1955_verdict": d1955["verdict"],
            "q1956_verdict": d1956["verdict"],
        },
        "flags": {
            "stage_b_closed": stage_b_closed,
            "triad_strict_pass_q1946": triad_strict_pass,
            "minimal_repair_pass_q1955": repair_pass,
            "two_state_repair_pass_q1956": two_state_repair_pass,
        },
        "blocker_ranking": blockers,
        "readiness": readiness,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1957: STRICT STAGE-C+INFO GATE V2",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Readiness: **{readiness}**",
        "",
        "## Flags",
        f"- stage_b_closed: {stage_b_closed}",
        f"- triad_strict_pass_q1946: {triad_strict_pass}",
        f"- minimal_repair_pass_q1955: {repair_pass}",
        f"- two_state_repair_pass_q1956: {two_state_repair_pass}",
        "",
        "## Top Blockers",
    ]
    for row in blockers[:8]:
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

    print(f"[QW-1957] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1957] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1957] readiness={readiness} verdict={verdict}")


if __name__ == "__main__":
    main()

