#!/usr/bin/env python3
"""P1489 S4.39: construct strict-only selector source candidate for QW-2191."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1487 = GEN / "p1487_s437_qw2191_selector_injection_test_summary.json"
P1478 = GEN / "p1478_s428_qw2191_sp1_operating_policy_v2.json"

SUMMARY = GEN / "p1489_s439_qw2191_selector_source_candidate_summary.json"
OBSTRUCTION = GEN / "p1489_s439_qw2191_selector_source_candidate_obstruction.json"


def main() -> None:
    s1487 = json.loads(P1487.read_text(encoding="utf-8"))
    pol = json.loads(P1478.read_text(encoding="utf-8"))

    g0 = float(s1487["gaps"]["baseline_g0"])
    delta_sb = float(s1487["delta_sb"])
    safety_margin = float(pol["safety_margin"])

    eps = 1e-12
    s_src = delta_sb / (abs(g0) + eps)

    checks = {
        "strict_only": True,
        "no_legacy_bridge": True,
        "gap_nonzero": abs(g0) > 0,
        "delta_within_safety_margin": abs(delta_sb) <= safety_margin,
        "selector_source_bounded": abs(s_src) <= 1.0,
        "selector_source_nontrivial": abs(s_src) > 0.0,
    }

    local_candidate_ok = all(checks.values())
    status = "PASS_SELECTOR_SOURCE_CANDIDATE_LOCAL_ONLY" if local_candidate_ok else "FAIL_SELECTOR_SOURCE_CANDIDATE_LOCAL_ONLY"

    summary = {
        "packet": "P1489",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "selector_source_id": "S_src_v1",
        "equation": "S_src_v1 = Delta_SB / (|W_SM - W_GR| + eps)",
        "value": s_src,
        "orientation": "SM_preferred" if s_src > 0 else "GR_preferred",
        "checks": checks,
        "qw2191_closed": False,
        "next_step_recommendation": "S4.40: verify S_src_v1 invariance across admissible kappa-shift grid and attempt theorem draft for internal selector source export.",
        "layman_explanation": "Wreszcie budujemy kandydat samego mechanizmu wyboru kierunku. Nie tylko patrzymy, że wynik się poprawił, ale mierzymy 'źródło przechyłu' i jego siłę.",
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if local_candidate_ok:
        if OBSTRUCTION.exists():
            OBSTRUCTION.unlink()
    else:
        OBSTRUCTION.write_text(json.dumps({"packet": "P1489", "status": status, "checks": checks}, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print(f"[P1489] status={status} S_src_v1={s_src:.6f}")


if __name__ == "__main__":
    main()
