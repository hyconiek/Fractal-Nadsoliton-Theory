#!/usr/bin/env python3
"""P1495 S4.45: quantified theorem draft and independent verifier for QW-2191."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1491 = GEN / "p1491_s441_qw2191_kappa_robustness_sweep_summary.json"
P1494 = GEN / "p1494_s444_qw2191_contradiction_branch_summary.json"

SUMMARY = GEN / "p1495_s445_qw2191_quantified_theorem_draft_summary.json"


def independent_verify(rows: list[dict], kmin: float, kmax: float, margin: float) -> dict:
    target = [r for r in rows if kmin <= float(r["kappa"]) <= kmax]
    cond1 = all(abs(float(r["delta_sb"])) <= margin for r in target)
    cond2 = all(float(r["gap_after"]) < float(r["gap_before"]) for r in target)
    cond3 = all(float(r["gap_after"]) > 0 for r in target)
    return {"cond1": cond1, "cond2": cond2, "cond3": cond3, "n_points": len(target)}


def main() -> None:
    s1491 = json.loads(P1491.read_text(encoding="utf-8"))
    s1494 = json.loads(P1494.read_text(encoding="utf-8"))

    kmin, kmax = [float(x) for x in s1494["robust_kappa_range"]]
    margin = float(s1491["safety_margin"])

    verify = independent_verify(s1491["rows"], kmin, kmax, margin)
    theorem_holds_local = verify["cond1"] and verify["cond2"] and verify["cond3"] and verify["n_points"] > 0

    summary = {
        "packet": "P1495",
        "status": "PASS_QUANTIFIED_THEOREM_DRAFT_LOCAL_ONLY" if theorem_holds_local else "FAIL_QUANTIFIED_THEOREM_DRAFT_LOCAL_ONLY",
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "theorem_name": "QW2191_local_selector_consistency_quantified_v1",
        "quantified_domain": {"kappa_min": kmin, "kappa_max": kmax},
        "independent_verification": verify,
        "theorem_holds_local": theorem_holds_local,
        "qw2191_closed": False,
        "next_step_recommendation": "S4.46: attempt cross-provider replication (new data/provider class) to test whether quantified theorem persists beyond current local stream.",
        "layman_explanation": "Sformułowaliśmy twierdzenie w stylu 'dla każdego ustawienia w bezpiecznym zakresie' i sprawdziliśmy je niezależnie. Lokalnie działa, ale pełne zamknięcie wymaga niezależnych replikacji.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1495] status={summary['status']} n={verify['n_points']}")


if __name__ == "__main__":
    main()
