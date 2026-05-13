#!/usr/bin/env python3
"""P1485 S4.35: strict closure criterion for QW-2191 (no legacy bridge)."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1484 = GEN / "p1484_s434_qw2191_operator_level_witness_probe_summary.json"
P1467 = GEN / "p1467_s417_qw2191_selector_premise_registry_summary.json"

SUMMARY = GEN / "p1485_s435_qw2191_strict_closure_criterion_summary.json"
OBSTRUCTION = GEN / "p1485_s435_qw2191_strict_closure_criterion_obstruction.json"


def main() -> None:
    s1484 = json.loads(P1484.read_text(encoding="utf-8"))
    s1467 = json.loads(P1467.read_text(encoding="utf-8"))

    checks_1484 = s1484["checks"]
    witnesses_ok = bool(checks_1484["sm_witness_positive"] and checks_1484["gr_witness_positive"] and checks_1484["mix_below_sm"] and checks_1484["mix_below_gr"])

    # registry currently tracks premise as non-strict unless internally proven
    premise_status = str(s1467.get("premise_status", "UNKNOWN"))
    strict_internal_selector_source = False
    explicit_symmetry_breaking_premise_verified = premise_status == "STRICT_INTERNAL_PROVEN"

    closure_criterion = {
        "needs_selector_source_or_verified_symmetry_breaking": True,
        "strict_internal_selector_source": strict_internal_selector_source,
        "explicit_symmetry_breaking_premise_verified": explicit_symmetry_breaking_premise_verified,
    }

    qw2191_closed = bool(
        witnesses_ok and (
            closure_criterion["strict_internal_selector_source"]
            or closure_criterion["explicit_symmetry_breaking_premise_verified"]
        )
    )

    status = "PASS_QW2191_STRICT_CLOSED" if qw2191_closed else "OPEN_QW2191_STRICT_NOT_CLOSED"

    summary = {
        "packet": "P1485",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "route": "strict_only_F_to_LSM_plus_LGR",
        "witnesses_ready": witnesses_ok,
        "closure_criterion": closure_criterion,
        "qw2191_closed": qw2191_closed,
        "next_step_recommendation": "S4.36: construct and test one physically explicit symmetry-breaking selector premise directly in operator equations for SM-like and GR-like channels.",
        "layman_explanation": "Mamy już dobre sygnały z dwóch sektorów, ale to jeszcze nie domyka problemu wyboru. Potrzebny jest dodatkowy, fizycznie konkretny mechanizm wyboru jednej gałęzi rozwiązania.",
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if qw2191_closed:
        if OBSTRUCTION.exists():
            OBSTRUCTION.unlink()
    else:
        obstruction = {
            "packet": "P1485",
            "status": status,
            "rule": "strict closure criterion not satisfied",
            "missing": [
                "strict_internal_selector_source OR explicit_symmetry_breaking_premise_verified"
            ],
        }
        OBSTRUCTION.write_text(json.dumps(obstruction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print(f"[P1485] status={status} witnesses_ready={witnesses_ok}")


if __name__ == "__main__":
    main()
