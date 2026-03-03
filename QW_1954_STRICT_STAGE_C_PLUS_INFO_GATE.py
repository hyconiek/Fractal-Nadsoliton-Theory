#!/usr/bin/env python3
"""
QW-1954: Strict Stage-C+Info gate.

Integrates:
- QW-1946 (mass+flavor+GW, single-kernel v2),
- QW-1952 (information-channel de-degeneracy operator),
- QW-1953 (two-state internal observer).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1954_strict_stage_c_plus_info_gate.json"
OUT_MD = ROOT / "RAPORT_QW1954_STRICT_STAGE_C_PLUS_INFO_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def rel_miss_hi(value: float, limit: float) -> float:
    return float(max(0.0, value - limit) / max(limit, 1e-12))


def rel_miss_lo(value: float, limit: float) -> float:
    return float(max(0.0, limit - value) / max(limit, 1e-12))


def main() -> None:
    d1946 = load("report_qw1946_final_single_kernel_gate_v2.json")
    d1952 = load("report_qw1952_information_channel_dedegeneracy_operator.json")
    d1953 = load("report_qw1953_two_state_internal_observer.json")

    stage_b_closed = bool(d1946["flags"]["stage_b_closed"])
    triad_strict_pass = bool(
        d1946["flags"]["stage_b_closed"]
        and d1946["flags"]["mass_pass_q1945_primary"]
        and d1946["flags"]["strict_flavor_gw_pass_q1944_baseline"]
    )
    dedeg_pass = bool(d1952.get("dedeg_pass", False))
    two_state_pass = bool(d1953.get("two_state_pass", False))

    strict_all = bool(triad_strict_pass and dedeg_pass and two_state_pass)

    if strict_all:
        readiness = "TOE_STAGE_C_PLUS_INFO_CLOSED_STRICT"
        verdict = "STRICT_STAGE_C_PLUS_INFO_GATE_PASS"
        required_next = "LOCK_FULL_SINGLE_KERNEL_STACK_AND_PREPARE_EXTERNAL_CONFIRMATORY_EXECUTION"
    elif stage_b_closed:
        readiness = "TOE_STAGE_C_PLUS_INFO_BLOCKED"
        verdict = "STRICT_STAGE_C_PLUS_INFO_GATE_FAIL"
        required_next = "REWORK_INFORMATION_CHANNEL_DEGENERACY_AND_MASS_LINK_UNDER_FROZEN_KERNEL"
    else:
        readiness = "TOE_STAGE_B_BLOCKED"
        verdict = "STRICT_STAGE_C_PLUS_INFO_GATE_BLOCKED_BY_STAGE_B"
        required_next = "RECOVER_STAGE_B_FIRST"

    blockers: List[Dict[str, float]] = []
    for row in d1946.get("blocker_ranking", []):
        blockers.append({"metric": f"triad::{row['metric']}", "relative_miss": float(row["relative_miss"])})

    t52 = d1952["thresholds"]
    m52 = d1952["metrics"]["dual"]
    blockers.extend(
        [
            {
                "metric": "info1952::dual_accuracy",
                "relative_miss": rel_miss_lo(float(m52["accuracy"]), float(t52["dual_accuracy_min"])),
            },
            {
                "metric": "info1952::acc_gain_vs_control",
                "relative_miss": rel_miss_lo(
                    float(d1952["metrics"]["acc_gain_dual_vs_control"]),
                    float(t52["acc_gain_vs_control_min"]),
                ),
            },
            {
                "metric": "info1952::info_gain_vs_control",
                "relative_miss": rel_miss_lo(
                    float(d1952["metrics"]["info_gain_dual_vs_control"]),
                    float(t52["info_gain_vs_control_min"]),
                ),
            },
            {
                "metric": "info1952::channel_complementarity",
                "relative_miss": rel_miss_lo(
                    float(m52["channel_complementarity"]),
                    float(t52["channel_complementarity_min"]),
                ),
            },
        ]
    )

    t53 = d1953["thresholds"]
    m53 = d1953["metrics"]["closed"]
    blockers.extend(
        [
            {
                "metric": "info1953::closed_accuracy",
                "relative_miss": rel_miss_lo(float(m53["accuracy"]), float(t53["closed_accuracy_min"])),
            },
            {
                "metric": "info1953::closed_acc_gain_vs_open",
                "relative_miss": rel_miss_lo(
                    float(d1953["metrics"]["closed_acc_gain_vs_open"]),
                    float(t53["closed_acc_gain_vs_open_min"]),
                ),
            },
            {
                "metric": "info1953::closed_info_gain_vs_control",
                "relative_miss": rel_miss_lo(
                    float(d1953["metrics"]["closed_info_gain_vs_control"]),
                    float(t53["closed_info_gain_vs_control_min"]),
                ),
            },
            {
                "metric": "info1953::channel_complementarity",
                "relative_miss": rel_miss_lo(
                    float(m53["channel_complementarity"]),
                    float(d1952["thresholds"]["channel_complementarity_min"]),
                ),
            },
        ]
    )
    blockers = sorted(blockers, key=lambda x: x["relative_miss"], reverse=True)

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "q1946_readiness": d1946["readiness"],
            "q1952_verdict": d1952["verdict"],
            "q1953_verdict": d1953["verdict"],
        },
        "flags": {
            "stage_b_closed": stage_b_closed,
            "triad_strict_pass_q1946": triad_strict_pass,
            "info_dedeg_pass_q1952": dedeg_pass,
            "info_two_state_pass_q1953": two_state_pass,
        },
        "blocker_ranking": blockers,
        "readiness": readiness,
        "verdict": verdict,
        "required_next_step": required_next,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1954: STRICT STAGE-C+INFO GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Readiness: **{readiness}**",
        "",
        "## Flags",
        f"- stage_b_closed: {stage_b_closed}",
        f"- triad_strict_pass_q1946: {triad_strict_pass}",
        f"- info_dedeg_pass_q1952: {dedeg_pass}",
        f"- info_two_state_pass_q1953: {two_state_pass}",
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

    print(f"[QW-1954] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1954] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1954] readiness={readiness} verdict={verdict}")


if __name__ == "__main__":
    main()

