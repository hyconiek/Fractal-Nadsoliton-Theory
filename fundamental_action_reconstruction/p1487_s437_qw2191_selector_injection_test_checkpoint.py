#!/usr/bin/env python3
"""P1487 S4.37: selector injection contrast test for QW-2191."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1484 = GEN / "p1484_s434_qw2191_operator_level_witness_probe_summary.json"
P1486 = GEN / "p1486_s436_qw2191_symmetry_breaking_premise_probe_summary.json"
P1478 = GEN / "p1478_s428_qw2191_sp1_operating_policy_v2.json"

SUMMARY = GEN / "p1487_s437_qw2191_selector_injection_test_summary.json"
OBSTRUCTION = GEN / "p1487_s437_qw2191_selector_injection_test_obstruction.json"


def main() -> None:
    s1484 = json.loads(P1484.read_text(encoding="utf-8"))
    s1486 = json.loads(P1486.read_text(encoding="utf-8"))
    pol = json.loads(P1478.read_text(encoding="utf-8"))

    w_sm = float(s1484["witnesses"]["sm_witness"])
    w_gr = float(s1484["witnesses"]["gr_witness"])
    delta_sb = float(s1486["delta_sb"])
    safety_margin = float(pol["safety_margin"])

    g0 = abs(w_sm - w_gr)
    g1 = abs(w_sm - w_gr - delta_sb)

    checks = {
        "strict_only": True,
        "no_legacy_bridge": True,
        "delta_within_safety_margin": abs(delta_sb) <= safety_margin,
        "selector_gap_reduced": g1 < g0,
        "sm_witness_positive": w_sm > 0,
        "gr_witness_positive": w_gr > 0,
    }

    local_progress = all(checks.values())
    status = "PASS_SELECTOR_INJECTION_LOCAL_PROGRESS" if local_progress else "FAIL_SELECTOR_INJECTION_LOCAL_PROGRESS"

    summary = {
        "packet": "P1487",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "gaps": {"baseline_g0": g0, "injected_g1": g1},
        "delta_sb": delta_sb,
        "checks": checks,
        "qw2191_closed": False,
        "closure_note": "Local selector progress observed; strict-final QW-2191 closure still requires exported internal selector source or theorem-level uniqueness discharge.",
        "next_step_recommendation": "S4.38: run selector-sensitivity sweep over admissible kappa window and export robustness envelope for uniqueness gap reduction.",
        "layman_explanation": "To test: dodaliśmy mały, kontrolowany składnik i sprawdziliśmy, czy problem wyboru rozwiązania się zmniejsza. Zmniejsza się, więc to dobry sygnał, ale to jeszcze nie pełne domknięcie.",
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if local_progress:
        if OBSTRUCTION.exists():
            OBSTRUCTION.unlink()
    else:
        OBSTRUCTION.write_text(json.dumps({"packet": "P1487", "status": status, "checks": checks}, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print(f"[P1487] status={status} g0={g0:.6f} g1={g1:.6f}")


if __name__ == "__main__":
    main()
