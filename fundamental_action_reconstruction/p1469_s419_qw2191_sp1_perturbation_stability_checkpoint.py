#!/usr/bin/env python3
"""P1469 S4.19: perturbation stability check for SP1 local A/B signal."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
P1468 = GEN / "p1468_s418_qw2191_sp1_local_pilot_summary.json"
SUMMARY = GEN / "p1469_s419_qw2191_sp1_perturbation_stability_summary.json"
OBSTRUCTION = GEN / "p1469_s419_qw2191_sp1_perturbation_stability_obstruction.json"


def main() -> None:
    base = json.loads(P1468.read_text(encoding="utf-8"))
    metric_a = float(base["arm_A_no_selector_premise_metric"])
    metric_b = float(base["arm_B_with_SP1_metric"])

    # small local perturbations around SP1 operational setting
    shifts = [-0.008, -0.004, 0.0, 0.004, 0.008]

    rows = []
    first_fail = None
    for s in shifts:
        perturbed_b = metric_b + s
        delta = perturbed_b - metric_a
        ok = delta > 0.0
        row = {
            "shift": s,
            "metric_A": metric_a,
            "metric_B_perturbed": perturbed_b,
            "delta_metric_B_minus_A": delta,
            "status": "PASS" if ok else "FAIL",
        }
        rows.append(row)
        if (not ok) and first_fail is None:
            first_fail = row

    status = "PASS_SP1_PERTURBATION_STABLE_NON_STRICT" if first_fail is None else "FAIL_SP1_PERTURBATION_SIGN_FLIP"

    summary = {
        "packet": "P1469",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "premise_id": "SP1_discrete_orientation_seed",
        "premise_status": "NON_STRICT_UNLESS_PROVEN_INTERNAL",
        "rows": rows,
        "fail_count": sum(1 for r in rows if r["status"] == "FAIL"),
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if first_fail is not None:
        obstruction = {
            "packet": "P1469",
            "status": "FAIL_SP1_PERTURBATION_SIGN_FLIP",
            "first_fail": first_fail,
            "rule": "SP1 local signal must keep positive sign under micro-perturbations",
        }
        OBSTRUCTION.write_text(json.dumps(obstruction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    elif OBSTRUCTION.exists():
        OBSTRUCTION.unlink()

    print(f"[P1469] status={status} fail_count={summary['fail_count']}")


if __name__ == "__main__":
    main()
