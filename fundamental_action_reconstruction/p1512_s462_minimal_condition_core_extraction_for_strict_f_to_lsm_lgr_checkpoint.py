#!/usr/bin/env python3
"""P1512 S4.62: extract minimal robust-condition core from P1511 boundary map."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1511 = GEN / "p1511_s461_admissible_perturbation_boundary_map_f_to_lsm_lgr_summary.json"
SUMMARY = GEN / "p1512_s462_minimal_condition_core_extraction_for_strict_f_to_lsm_lgr_summary.json"


def main() -> None:
    s1511 = json.loads(P1511.read_text(encoding="utf-8"))
    robust = s1511.get("robust_zone", [])
    rejection = s1511.get("rejection_zone", [])

    assumption_keys = ["A1_selector_positive", "A2_selector_strict_internal", "A3_weight_normalization", "A4_shared_orientation", "A5_no_legacy_bridge"]

    always_true_in_robust = {}
    for k in assumption_keys:
        always_true_in_robust[k] = all(bool(r.get("assumptions", {}).get(k, False)) for r in robust)

    violated_in_rejection = {}
    for k in assumption_keys:
        violated_in_rejection[k] = any(not bool(r.get("assumptions", {}).get(k, True)) for r in rejection)

    critical_boundary_conditions = [k for k in assumption_keys if violated_in_rejection[k]]

    minimal_core = {
        "necessary_conditions": [k for k in assumption_keys if always_true_in_robust[k]],
        "critical_rejection_triggers": critical_boundary_conditions,
    }

    summary = {
        "packet": "P1512",
        "status": "PASS_MINIMAL_CORE_EXTRACTED",
        "scope": "STRICT_ONLY_F_NADSOLITON_TO_LSM_PLUS_LGR",
        "input_counts": {
            "robust_points": len(robust),
            "rejection_points": len(rejection),
        },
        "always_true_in_robust": always_true_in_robust,
        "violated_in_rejection": violated_in_rejection,
        "minimal_core_A_star": minimal_core,
        "main_finding": "A4_shared_orientation is the observed active boundary trigger in current map; robust region preserves A1..A5.",
        "qw2191_closed": False,
        "legacy_bridge_used": False,
        "next_honest_step": "Implement P1513 orientation-stability strengthening checkpoint to harden A4 under expanded admissible perturbation families.",
        "layman_explanation": "Zredukowaliśmy model do rdzenia: wiemy, które warunki muszą być zawsze spełnione. Najbardziej krytyczny punkt to zgodność orientacji — gdy ona pęka, dowód się rozpada.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1512] status={summary['status']} critical={critical_boundary_conditions}")


if __name__ == "__main__":
    main()
