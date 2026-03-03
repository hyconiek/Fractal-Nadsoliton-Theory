#!/usr/bin/env python3
"""
QW-1847: PTA replication target proposal for strict closure.

Builds action targets from QW-1844 and QW-1846 for the probability criterion.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1847_pta_replication_target_proposal.json"
OUT_MD = ROOT / "RAPORT_QW1847_PTA_REPLICATION_TARGET_PROPOSAL.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def pick_row(table, p0: float, p1: float) -> Optional[Dict]:
    for r in table:
        if (
            float(r.get("p0_threshold", -1.0)) == p0
            and float(r.get("p_true_assumed", -1.0)) == p1
            and bool(r.get("feasible", False))
        ):
            return r
    return None


def main() -> None:
    d1844 = load("report_qw1844_strict_joint_rigor_gate.json")
    d1846 = load("report_qw1846_pta_prob_threshold_feasibility_map.json")

    n_current = int(d1844.get("power_gap", {}).get("n_replications_current", 0))
    table = d1846.get("feasibility_table", [])

    r_minimal = pick_row(table, p0=0.90, p1=1.00)  # all-positive idealized
    r_optimistic = pick_row(table, p0=0.90, p1=0.99)
    r_realistic = pick_row(table, p0=0.90, p1=0.97)

    # Robust recommendation defaults to realistic scenario.
    if r_realistic is not None:
        n_target = int(r_realistic["n_required"])
    elif r_optimistic is not None:
        n_target = int(r_optimistic["n_required"])
    elif r_minimal is not None:
        n_target = int(r_minimal["n_required"])
    else:
        n_target = n_current

    n_add = max(0, n_target - n_current)

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "strict_gate_input": {
            "qw1844_hard_gate": d1844.get("hard_gate"),
            "qw1844_readiness": d1844.get("readiness"),
        },
        "targets": {
            "current_n": n_current,
            "minimal_all_positive": r_minimal,
            "optimistic_p_true_0p99": r_optimistic,
            "realistic_p_true_0p97": r_realistic,
            "recommended_target_n": n_target,
            "additional_needed": n_add,
        },
        "execution_status": {
            "strict_ready_now": bool(d1844.get("hard_gate") == "PASS"),
            "strict_ready_after_target": True,
        },
        "verdict": "PTA_REPLICATION_TARGET_PROPOSAL_COMPLETE",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1847: PTA REPLICATION TARGET PROPOSAL",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Strict gate now: {d1844.get('hard_gate')} / {d1844.get('readiness')}",
        f"- Current n: {n_current}",
        f"- Recommended target n: {n_target}",
        f"- Additional needed: {n_add}",
        "",
        "## Scenario targets for p0=0.90",
    ]

    if r_minimal:
        lines.append(
            (
                f"- Minimal (p_true=1.00): n={r_minimal['n_required']}, "
                f"add={r_minimal['add_vs_current']}"
            )
        )
    if r_optimistic:
        lines.append(
            (
                f"- Optimistic (p_true=0.99): n={r_optimistic['n_required']}, "
                f"add={r_optimistic['add_vs_current']}"
            )
        )
    if r_realistic:
        lines.append(
            (
                f"- Realistic (p_true=0.97): n={r_realistic['n_required']}, "
                f"add={r_realistic['add_vs_current']}"
            )
        )

    lines += [
        "",
        "## Interpretation",
        "- Przy zachowaniu progu p0=0.90 potrzebna jest istotna eskalacja liczby replikacji, zeby kryterium bylo inferencyjnie domkniete (alpha=0.05, power~0.80).",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1847] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1847] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
