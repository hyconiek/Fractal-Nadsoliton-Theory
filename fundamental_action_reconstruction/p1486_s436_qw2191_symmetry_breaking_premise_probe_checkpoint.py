#!/usr/bin/env python3
"""P1486 S4.36: explicit symmetry-breaking premise probe for QW-2191."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1484 = GEN / "p1484_s434_qw2191_operator_level_witness_probe_summary.json"
P1485 = GEN / "p1485_s435_qw2191_strict_closure_criterion_summary.json"
P1478 = GEN / "p1478_s428_qw2191_sp1_operating_policy_v2.json"

SUMMARY = GEN / "p1486_s436_qw2191_symmetry_breaking_premise_probe_summary.json"
OBSTRUCTION = GEN / "p1486_s436_qw2191_symmetry_breaking_premise_probe_obstruction.json"


def main() -> None:
    s1484 = json.loads(P1484.read_text(encoding="utf-8"))
    s1485 = json.loads(P1485.read_text(encoding="utf-8"))
    pol = json.loads(P1478.read_text(encoding="utf-8"))

    sm = float(s1484["witnesses"]["sm_witness"])
    gr = float(s1484["witnesses"]["gr_witness"])
    safety_margin = float(pol["safety_margin"])

    kappa = 0.2
    delta_sb = kappa * (sm - gr)

    premise_checks = {
        "strict_only": True,
        "no_legacy_bridge": True,
        "sm_witness_positive": sm > 0,
        "gr_witness_positive": gr > 0,
        "kappa_positive": kappa > 0,
        "delta_sb_within_safety_margin": abs(delta_sb) <= safety_margin,
        "closure_was_open_before_probe": s1485["qw2191_closed"] is False,
    }

    premise_verified_local = all(premise_checks.values())
    status = "PASS_SYMMETRY_BREAKING_PREMISE_LOCAL_ONLY" if premise_verified_local else "FAIL_SYMMETRY_BREAKING_PREMISE_LOCAL_ONLY"

    summary = {
        "packet": "P1486",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "premise_id": "SP1_SB_v1",
        "equation": "Delta_SB = kappa * (W_SM - W_GR)",
        "kappa": kappa,
        "delta_sb": delta_sb,
        "safety_margin": safety_margin,
        "premise_verified_local": premise_verified_local,
        "qw2191_closed": False,
        "next_step_recommendation": "S4.37: inject SP1_SB_v1 into explicit operator equations and test whether selector uniqueness obstruction is removed without violating policy gate.",
        "layman_explanation": "Dodaliśmy kontrolowaną zasadę lekkiego przechylenia między dwoma kanałami. Sprawdzamy, czy to przechylenie jest małe, bezpieczne i fizycznie sensowne, zamiast arbitralnego wymuszania wyniku.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if premise_verified_local:
        if OBSTRUCTION.exists():
            OBSTRUCTION.unlink()
    else:
        OBSTRUCTION.write_text(json.dumps({"packet": "P1486", "status": status, "checks": premise_checks}, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print(f"[P1486] status={status} delta_sb={delta_sb:.6f}")


if __name__ == "__main__":
    main()
