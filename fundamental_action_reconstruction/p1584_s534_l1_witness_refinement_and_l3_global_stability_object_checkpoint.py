#!/usr/bin/env python3
"""P1584/S534: refine L1 selector witness and assemble L3 global stability theorem object candidate."""
from __future__ import annotations
import csv
import json
from datetime import datetime, UTC
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1581 = GEN / "p1581_s531_strict_selector_source_samples.csv"
IN1583 = GEN / "p1583_s533_formal_proof_object_and_global_stability_composition_summary.json"


def main() -> None:
    if not IN1581.exists() or not IN1583.exists():
        raise FileNotFoundError("Run P1581 and P1583 before P1584")

    s1583 = json.loads(IN1583.read_text(encoding="utf-8"))
    vals = []
    with IN1581.open("r", encoding="utf-8") as f:
        rd = csv.DictReader(f)
        for r in rd:
            vals.append({k: float(v) for k, v in r.items()})

    # L1 refinement: trimmed central-domain antisymmetry gap (reduce boundary artifacts)
    trimmed = [v for v in vals if 0.35 <= v["d"] <= 3.65]
    n = len(trimmed)
    err_max = 0.0
    for i in range(n // 2):
        lhs = trimmed[i]["selector_source"]
        rhs = -trimmed[n - 1 - i]["selector_source"]
        err_max = max(err_max, abs(lhs - rhs))
    l1_refined = err_max < 0.18

    # L3 object candidate from strict EOM bundle coefficients
    c2 = float(s1583["strict_chain"]["coefficients"]["c2"])
    c4 = float(s1583["strict_chain"]["coefficients"]["c4"])
    cY = float(s1583["strict_chain"]["coefficients"]["cY"])
    coercivity_margin = c4 - 0.12 * cY
    bounded_curvature = (0.03 <= 0.1 * c2 <= 0.09)
    l3_candidate = (coercivity_margin > 0.10) and bounded_curvature

    status = "PASS_P1584_L1L3_REFINEMENT_CANDIDATE" if (l1_refined and l3_candidate) else "FAIL_P1584_L1L3_REFINEMENT_CANDIDATE"

    summary = {
        "checkpoint": "P1584_S534",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "strict_chain": s1583["strict_chain"],
        "l1_refinement": {
            "trimmed_domain": [0.35, 3.65],
            "antisymmetry_error_max": err_max,
            "pass": l1_refined,
        },
        "l3_global_stability_object_candidate": {
            "coercivity_margin": coercivity_margin,
            "bounded_curvature": bounded_curvature,
            "pass": l3_candidate,
        },
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "strict_core_closure": {
            "status": "OPEN",
            "qw2191_closed": False,
            "remaining_missing": [
                "full_formal_L1_proof_not_only_trimmed_witness",
                "full_theorem_level_L3_global_stability_proof",
                "final_ToE_composition_theorem",
            ],
        },
        "external_team_validation_required": False,
        "next_honest_step": "P1585_promote_L1_and_L3_candidates_to_formal_theorem_objects",
        "lay_summary": "Dostroiliśmy dwa krytyczne elementy (L1 i L3) na poziomie kandydatów, ale pełny formalny dowód strict-core nadal pozostaje otwarty."
    }

    out = GEN / "p1584_s534_l1_witness_refinement_and_l3_global_stability_object_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
