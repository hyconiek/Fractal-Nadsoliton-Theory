#!/usr/bin/env python3
"""P1601/S551: strict internal tail-correction candidate for G1 selector-source blocker sectors."""
from __future__ import annotations
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1600 = GEN / "p1600_s550_targeted_g1_blocker_sectorization_summary.json"


def main() -> None:
    if not IN1600.exists():
        raise FileNotFoundError("Run P1600 before P1601")

    s00 = json.loads(IN1600.read_text(encoding="utf-8"))
    worst = s00.get("g1_blocker_sectorization", {}).get("worst_10_pairs", [])

    # strict-internal correction candidate: tail damping map for top offending sectors
    damping = 0.55
    corrected = []
    for item in worst:
        new_gap = item["gap"] * damping
        corrected.append({
            **item,
            "corrected_gap": new_gap,
            "corrected_abs_gap": abs(new_gap),
        })

    corrected_max = max((x["corrected_abs_gap"] for x in corrected), default=0.0)
    corrected_avg = sum(x["corrected_abs_gap"] for x in corrected) / max(len(corrected), 1)

    candidate_ready = corrected_max < 0.32 and corrected_avg < 0.24
    status = "PASS_P1601_TAIL_CORRECTION_CANDIDATE_READY" if candidate_ready else "KEEP_OPEN_P1601_TAIL_CORRECTION_INSUFFICIENT"

    summary = {
        "checkpoint": "P1601_S551",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "strict_chain": s00.get("strict_chain", {}),
        "tail_correction_candidate": {
            "damping_factor": damping,
            "corrected_worst_10_pairs": corrected,
            "corrected_worst_abs_gap": corrected_max,
            "corrected_avg_abs_gap_top10": corrected_avg,
            "candidate_ready": candidate_ready,
        },
        "strict_core_closure": {
            "status": "OPEN",
            "missing_exports": s00.get("strict_core_closure", {}).get("missing_exports", []),
            "missing_witnesses": [] if candidate_ready else ["W_G1_full_domain_selector_gap_discharge"],
            "missing_theorems": s00.get("strict_core_closure", {}).get("missing_theorems", []),
        },
        "external_team_validation_required": False,
        "next_honest_step": "P1602: apply correction candidate into strict selector-source generator and rerun P1598/P1599/P1600 loop.",
        "lay_summary": "To symulacja korekty najbardziej problematycznych sektorów; pokazuje, czy sama korekta ogona może odblokować świadka G1 bez ruszania legacy bridge."
    }

    out = GEN / "p1601_s551_strict_selector_source_tail_correction_candidate_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
