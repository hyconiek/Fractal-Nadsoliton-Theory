#!/usr/bin/env python3
"""P1598/S548: strict-only refinement of G1 full-domain selector witness for final closure path."""
from __future__ import annotations
import csv
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1581 = GEN / "p1581_s531_strict_selector_source_samples.csv"
IN1597 = GEN / "p1597_s547_final_g3_theorem_composition_object_summary.json"


def main() -> None:
    if not IN1581.exists() or not IN1597.exists():
        raise FileNotFoundError("Run P1581 and P1597 before P1598")

    s97 = json.loads(IN1597.read_text(encoding="utf-8"))
    vals = []
    with IN1581.open("r", encoding="utf-8") as f:
        rd = csv.DictReader(f)
        for r in rd:
            vals.append({k: float(v) for k, v in r.items()})

    # refined odd-symmetry witness under strict selector-source discipline
    n = len(vals)
    pair_gaps = []
    for i in range(n // 2):
        pair_gaps.append(vals[i]["selector_source"] + vals[n - 1 - i]["selector_source"])

    emax = max(abs(g) for g in pair_gaps) if pair_gaps else 0.0
    l2 = (sum(g * g for g in pair_gaps) / max(len(pair_gaps), 1)) ** 0.5
    mean = sum(pair_gaps) / max(len(pair_gaps), 1)

    # stricter theorem-facing thresholds
    witness_pass = (emax < 0.40) and (l2 < 0.20) and (abs(mean) < 0.03)
    status = "PASS_P1598_G1_WITNESS_REFINED" if witness_pass else "KEEP_OPEN_P1598_G1_WITNESS_STILL_OPEN"

    summary = {
        "checkpoint": "P1598_S548",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "strict_chain": s97.get("strict_chain", {}),
        "g1_refined_witness": {
            "full_domain_error_max": emax,
            "full_domain_error_l2": l2,
            "full_domain_error_mean": mean,
            "pass": witness_pass,
        },
        "strict_core_closure": {
            "status": "OPEN",
            "missing_exports": s97.get("strict_core_closure", {}).get("missing_exports", []),
            "missing_witnesses": [] if witness_pass else ["W_G1_full_domain_selector_gap_discharge"],
            "missing_theorems": s97.get("strict_core_closure", {}).get("missing_theorems", []),
        },
        "external_team_validation_required": False,
        "next_honest_step": "P1599: replay P1597 with refined G1 witness imported as discharged/non-discharged gate.",
        "lay_summary": "Ten krok wzmacnia świadka G1: sprawdza globalną symetrię źródła selektora surowszymi progami, aby uczciwie przygotować finalne domknięcie strict-core."
    }

    out = GEN / "p1598_s548_strict_g1_witness_refinement_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
