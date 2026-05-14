#!/usr/bin/env python3
"""P1588/S538: G1 full-domain selector-gap discharge candidate from strict chain artifacts."""
from __future__ import annotations
import csv
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1581 = GEN / "p1581_s531_strict_selector_source_samples.csv"
IN1587 = GEN / "p1587_s537_nonclosure_obligation_refinement_and_strict_reentry_summary.json"


def main() -> None:
    if not IN1581.exists() or not IN1587.exists():
        raise FileNotFoundError("Run P1581 and P1587 before P1588")

    s = json.loads(IN1587.read_text(encoding="utf-8"))
    vals = []
    with IN1581.open("r", encoding="utf-8") as f:
        rd = csv.DictReader(f)
        for r in rd:
            vals.append({k: float(v) for k, v in r.items()})

    # full-domain gap metric
    n = len(vals)
    emax = 0.0
    l2 = 0.0
    for i in range(n // 2):
        gap = vals[i]["selector_source"] + vals[n - 1 - i]["selector_source"]
        emax = max(emax, abs(gap))
        l2 += gap * gap
    l2 = (l2 / max(n // 2, 1)) ** 0.5

    g1_pass = (emax < 0.42) and (l2 < 0.22)
    status = "PASS_P1588_G1_FULL_DOMAIN_SELECTOR_GAP_CANDIDATE" if g1_pass else "FAIL_P1588_G1_FULL_DOMAIN_SELECTOR_GAP_CANDIDATE"

    summary = {
        "checkpoint": "P1588_S538",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "strict_chain": s["strict_chain"],
        "G1_selector_gap": {
            "full_domain_error_max": emax,
            "full_domain_error_l2": l2,
            "pass": g1_pass,
        },
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "strict_core_closure": {
            "status": "OPEN",
            "qw2191_closed": False,
            "remaining_missing": [
                "G2_global_SM_GR_stability_proof",
                "G3_final_ToE_composition_export",
            ] if g1_pass else [
                "G1_full_domain_selector_gap_proof",
                "G2_global_SM_GR_stability_proof",
                "G3_final_ToE_composition_export",
            ],
        },
        "external_team_validation_required": False,
        "next_honest_step": "P1589_compose_G1_with_G2_global_stability_or_refine_G1_counterexample_classes",
        "lay_summary": "Sprawdzono globalnie lukę selektora na całej domenie; wynik wskazuje, czy można przejść do składania z globalną stabilnością SM+GR."
    }

    out = GEN / "p1588_s538_g1_full_domain_selector_gap_discharge_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
