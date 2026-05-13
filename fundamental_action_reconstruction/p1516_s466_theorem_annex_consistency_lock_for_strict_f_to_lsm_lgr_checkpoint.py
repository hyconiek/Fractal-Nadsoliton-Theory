#!/usr/bin/env python3
"""P1516 S4.66: enforce theorem-annex consistency lock via robust orientation envelope."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1515 = GEN / "p1515_s465_robust_orientation_envelope_export_for_strict_f_to_lsm_lgr_summary.json"
SUMMARY = GEN / "p1516_s466_theorem_annex_consistency_lock_for_strict_f_to_lsm_lgr_summary.json"


def main() -> None:
    s1515 = json.loads(P1515.read_text(encoding="utf-8"))
    robust_pairs = {
        (r["internal"], r["fmap"]) for r in s1515.get("orientation_envelope", {}).get("robust_pairs", [])
    }

    probes = [
        ("SM_preferred", "SM_preferred"),
        ("mixed", "mixed"),
        ("neutral", "neutral"),
        ("SM_preferred", "GR_preferred"),
    ]

    rows = []
    for p in probes:
        allowed = p in robust_pairs
        rows.append({
            "orientation_pair": {"internal": p[0], "fmap": p[1]},
            "theorem_eval_allowed": allowed,
            "lock_reason": "in_robust_envelope" if allowed else "outside_robust_envelope",
        })

    allowed_count = sum(1 for r in rows if r["theorem_eval_allowed"])
    blocked_count = len(rows) - allowed_count

    summary = {
        "packet": "P1516",
        "status": "PASS_THEOREM_ANNEX_CONSISTENCY_LOCK_ACTIVE",
        "scope": "STRICT_ONLY_F_NADSOLITON_TO_LSM_PLUS_LGR",
        "lock_rule": "THEOREM_EVAL_ALLOWED <=> orientation_pair in robust_envelope",
        "probe_results": rows,
        "counts": {
            "allowed": allowed_count,
            "blocked": blocked_count,
        },
        "qw2191_closed": False,
        "legacy_bridge_used": False,
        "next_honest_step": "Implement P1517 coupled-theorem rerun with lock enforced across all candidate orientation scenarios and export locked witness set.",
        "layman_explanation": "Wprowadziliśmy bezpiecznik: teoria jest liczona tylko dla ustawień kierunku z zatwierdzonej listy. To chroni przed przypadkami, które wyglądałyby dobrze, ale nie są stabilne.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1516] status={summary['status']} allowed={allowed_count} blocked={blocked_count}")


if __name__ == "__main__":
    main()
