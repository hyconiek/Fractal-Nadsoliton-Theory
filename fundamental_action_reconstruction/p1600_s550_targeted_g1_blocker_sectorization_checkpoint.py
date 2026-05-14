#!/usr/bin/env python3
"""P1600/S550: targeted strict G1 blocker sectorization for witness discharge planning."""
from __future__ import annotations
import csv
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1581 = GEN / "p1581_s531_strict_selector_source_samples.csv"
IN1599 = GEN / "p1599_s549_replay_g3_with_refined_g1_witness_summary.json"


def main() -> None:
    if not IN1581.exists() or not IN1599.exists():
        raise FileNotFoundError("Run P1581 and P1599 before P1600")

    s99 = json.loads(IN1599.read_text(encoding="utf-8"))
    rows = []
    with IN1581.open("r", encoding="utf-8") as f:
        rd = csv.DictReader(f)
        for r in rd:
            rows.append({k: float(v) for k, v in r.items()})

    n = len(rows)
    sector = []
    for i in range(n // 2):
        j = n - 1 - i
        gap = rows[i]["selector_source"] + rows[j]["selector_source"]
        sector.append({"i": i, "j": j, "gap": gap, "abs_gap": abs(gap)})

    sector.sort(key=lambda x: x["abs_gap"], reverse=True)
    top = sector[:10]
    avg_top = sum(x["abs_gap"] for x in top) / max(len(top), 1)
    max_gap = top[0]["abs_gap"] if top else 0.0

    # targeted discharge candidate rule: can only declare candidate if extreme tail is reduced
    tail_candidate_ready = max_gap < 0.32 and avg_top < 0.24

    summary = {
        "checkpoint": "P1600_S550",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "PASS_P1600_TARGETED_G1_SECTORS_IDENTIFIED",
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "strict_chain": s99.get("strict_chain", {}),
        "g1_blocker_sectorization": {
            "worst_10_pairs": top,
            "worst_abs_gap": max_gap,
            "avg_abs_gap_top10": avg_top,
            "tail_candidate_ready": tail_candidate_ready,
        },
        "strict_core_closure": {
            "status": "OPEN",
            "missing_exports": s99.get("strict_core_closure", {}).get("missing_exports", []),
            "missing_witnesses": [] if tail_candidate_ready else ["W_G1_full_domain_selector_gap_discharge"],
            "missing_theorems": s99.get("strict_core_closure", {}).get("missing_theorems", []),
        },
        "external_team_validation_required": False,
        "next_honest_step": "P1601: apply strict internal selector-source correction on worst sectors and rerun P1598/P1599.",
        "lay_summary": "Wskazaliśmy dokładnie, które obszary domeny najbardziej psują świadka G1; to pozwala celować poprawki zamiast powtarzać cały cykl w ciemno."
    }

    out = GEN / "p1600_s550_targeted_g1_blocker_sectorization_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
