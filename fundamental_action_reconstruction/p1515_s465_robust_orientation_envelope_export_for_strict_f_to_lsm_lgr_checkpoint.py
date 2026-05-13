#!/usr/bin/env python3
"""P1515 S4.65: export robust-orientation envelope from P1514 rerun."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1514 = GEN / "p1514_s464_strict_coupled_theorem_update_a4_star_and_expanded_contradiction_rerun_summary.json"
SUMMARY = GEN / "p1515_s465_robust_orientation_envelope_export_for_strict_f_to_lsm_lgr_summary.json"


def main() -> None:
    s1514 = json.loads(P1514.read_text(encoding="utf-8"))

    fails = s1514.get("expanded_contradiction_rerun", {}).get("fails", [])
    fail_pairs = {(f.get("orientation_internal"), f.get("orientation_fmap")) for f in fails}

    all_pairs = [
        ("SM_preferred", "SM_preferred"),
        ("GR_preferred", "GR_preferred"),
        ("mixed", "mixed"),
        ("neutral", "neutral"),
        ("SM_preferred", "GR_preferred"),
        ("GR_preferred", "SM_preferred"),
        ("mixed", "neutral"),
    ]

    robust_pairs = [p for p in all_pairs if p not in fail_pairs]
    rejected_pairs = [p for p in all_pairs if p in fail_pairs]

    envelope = {
        "criterion": "A4* = shared_orientation AND orientation_margin >= 0.5",
        "robust_pairs": [{"internal": a, "fmap": b} for a, b in robust_pairs],
        "rejected_pairs": [{"internal": a, "fmap": b} for a, b in rejected_pairs],
    }

    summary = {
        "packet": "P1515",
        "status": "PASS_ROBUST_ORIENTATION_ENVELOPE_EXPORTED",
        "scope": "STRICT_ONLY_F_NADSOLITON_TO_LSM_PLUS_LGR",
        "orientation_envelope": envelope,
        "counts": {
            "robust_pair_count": len(robust_pairs),
            "rejected_pair_count": len(rejected_pairs),
        },
        "admissibility_annex_ready": True,
        "qw2191_closed": False,
        "legacy_bridge_used": False,
        "next_honest_step": "Implement P1516 theorem-annex consistency lock: enforce envelope membership in all future strict coupled theorem evaluations.",
        "layman_explanation": "Stworzyliśmy oficjalną listę ustawień kierunku, które są bezpieczne dla modelu. To jak instrukcja: te konfiguracje wolno używać, a pozostałe trzeba odrzucać.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1515] status={summary['status']} robust={len(robust_pairs)} rejected={len(rejected_pairs)}")


if __name__ == "__main__":
    main()
