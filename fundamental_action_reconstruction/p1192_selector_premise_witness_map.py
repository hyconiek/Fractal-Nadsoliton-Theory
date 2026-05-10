#!/usr/bin/env python3
"""P1192 selector-premise -> witness-obligation map (strict-rigor)."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def main() -> None:
    p1190 = json.loads((GEN / "p1190_strict_promotion_nonclosure_guard_summary.json").read_text(encoding="utf-8"))

    guard_pass = bool(p1190.get("strict_promotion_nonclosure_guard_pass", False))

    obligations = [
        {
            "id": "W1_selector_uniqueness_bridge",
            "target": "QW-2191 strict-core uniqueness obstruction",
            "required_artifact": "formal selector source theorem OR explicit symmetry-breaking axiom with scoped non-strict tag",
            "pass_condition": "unique strict branch selected with reproducible witness",
            "status": "OPEN",
        },
        {
            "id": "W2_eom_role_assignment",
            "target": "explicit target eom-role witness family from P46 frontier",
            "required_artifact": "symbolic equality + numeric support across certified region",
            "pass_condition": "all declared role links validated",
            "status": "OPEN",
        },
        {
            "id": "W3_pairwise_mass_matching",
            "target": "pairwise matching set (e.g. m2_psi7=m2_psi10, m2_psi2=m2_psi5, m2_psi8=m2_psi11)",
            "required_artifact": "zero-defect certificates for each pair",
            "pass_condition": "all pair constraints pass without manual override",
            "status": "OPEN",
        },
        {
            "id": "W4_defect_polynomial_zeroes",
            "target": "direct m2 psi4 target action coefficient defect polynomial",
            "required_artifact": "explicit zero witness on common support",
            "pass_condition": "symbolic/numeric defect zero within tolerance",
            "status": "OPEN",
        },
    ]

    out = {
        "packet": "P1192",
        "as_of": "2026-05-10",
        "entry_condition_from_p1190": guard_pass,
        "candidate_status": "operationally_promoted_not_closed",
        "witness_obligations": obligations,
        "open_count": sum(1 for o in obligations if o["status"] == "OPEN"),
        "strict_closure_claim_allowed": False,
        "note": "Promotion is operational only until all witness obligations are discharged.",
    }

    out_path = GEN / "p1192_selector_premise_witness_map_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1192] open_witnesses={out['open_count']} strict_closure_claim_allowed={out['strict_closure_claim_allowed']} wrote {out_path}")


if __name__ == "__main__":
    main()
