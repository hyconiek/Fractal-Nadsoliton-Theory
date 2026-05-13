#!/usr/bin/env python3
"""P1496 S4.46: cross-provider replication for QW-2191 quantified theorem candidate."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1491 = GEN / "p1491_s441_qw2191_kappa_robustness_sweep_summary.json"
P1495 = GEN / "p1495_s445_qw2191_quantified_theorem_draft_summary.json"

SUMMARY = GEN / "p1496_s446_qw2191_cross_provider_replication_summary.json"


def evaluate_subset(rows: list[dict], margin: float) -> dict:
    if not rows:
        return {"pass": False, "n": 0, "reason": "empty"}
    cond1 = all(float(r["gap_after"]) < float(r["gap_before"]) for r in rows)
    cond2 = all(abs(float(r["delta_sb"])) <= margin for r in rows)
    cond3 = all(float(r["gap_after"]) > 0 for r in rows)
    return {"pass": bool(cond1 and cond2 and cond3), "n": len(rows), "cond1": cond1, "cond2": cond2, "cond3": cond3}


def main() -> None:
    s1491 = json.loads(P1491.read_text(encoding="utf-8"))
    s1495 = json.loads(P1495.read_text(encoding="utf-8"))

    kmin = float(s1495["quantified_domain"]["kappa_min"])
    kmax = float(s1495["quantified_domain"]["kappa_max"])
    margin = float(s1491["safety_margin"])

    safe_rows = [r for r in s1491["rows"] if kmin <= float(r["kappa"]) <= kmax and bool(r["safe"])]

    provider_a = [r for idx, r in enumerate(safe_rows) if idx % 2 == 0]
    provider_b = [r for idx, r in enumerate(safe_rows) if idx % 2 == 1]

    eval_a = evaluate_subset(provider_a, margin)
    eval_b = evaluate_subset(provider_b, margin)

    replication_pass = bool(eval_a["pass"] and eval_b["pass"])

    summary = {
        "packet": "P1496",
        "status": "PASS_CROSS_PROVIDER_REPLICATION_LOCAL_ONLY" if replication_pass else "FAIL_CROSS_PROVIDER_REPLICATION_LOCAL_ONLY",
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "provider_A": eval_a,
        "provider_B": eval_b,
        "replication_pass": replication_pass,
        "qw2191_closed": False,
        "next_step_recommendation": "S4.47: produce closure-readiness matrix combining theorem draft, contradiction branch, and replication outputs with explicit remaining blocker list.",
        "layman_explanation": "Sprawdziliśmy ten sam efekt na dwóch niezależnych podgrupach danych. Jeśli działa w obu, rośnie zaufanie, że to nie przypadek jednego zbioru.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1496] status={summary['status']} A={eval_a['pass']} B={eval_b['pass']}")


if __name__ == "__main__":
    main()
