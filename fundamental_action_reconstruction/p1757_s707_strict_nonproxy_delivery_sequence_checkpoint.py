#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1754 = GEN / "p1754_s704_strict_nonproxy_minimal_delivery_manifest_checkpoint.json"
IN1756 = GEN / "p1756_s706_strict_nonproxy_manifest_consistency_audit_checkpoint.json"
OUT = GEN / "p1757_s707_strict_nonproxy_delivery_sequence_checkpoint.json"


def main() -> None:
    p1754 = json.loads(IN1754.read_text(encoding="utf-8"))
    p1756 = json.loads(IN1756.read_text(encoding="utf-8"))

    manifest = p1754.get("minimal_nonproxy_delivery_manifest", {})
    consistency = p1756.get("consistency_result", {}).get("status", "UNKNOWN")

    sequence = [
        {
            "step": "D1_export_E_A_mu_nonproxy",
            "depends_on": [],
            "ready": manifest.get("M1_explicit_covariant_E_A_mu_expression_nonproxy", {}).get("status") == "MISSING",
        },
        {
            "step": "D2_export_E_H_nonproxy",
            "depends_on": [],
            "ready": manifest.get("M2_explicit_covariant_E_H_expression_nonproxy", {}).get("status") == "MISSING",
        },
        {
            "step": "D3_export_shared_background_family_contract",
            "depends_on": ["D1_export_E_A_mu_nonproxy", "D2_export_E_H_nonproxy"],
            "ready": manifest.get("M3_shared_background_family_contract", {}).get("status") == "MISSING",
        },
        {
            "step": "D4_export_boundary_term_control_clause",
            "depends_on": ["D3_export_shared_background_family_contract"],
            "ready": manifest.get("M4_boundary_term_control_clause", {}).get("status") == "MISSING",
        },
        {
            "step": "D5_finalize_boundary_control_contract",
            "depends_on": ["D4_export_boundary_term_control_clause"],
            "ready": manifest.get("M5_boundary_control_contract", {}).get("status") == "MISSING",
        },
        {
            "step": "D6_run_nonproxy_H1_4D",
            "depends_on": [
                "D1_export_E_A_mu_nonproxy",
                "D2_export_E_H_nonproxy",
                "D3_export_shared_background_family_contract",
                "D4_export_boundary_term_control_clause",
                "D5_finalize_boundary_control_contract",
            ],
            "ready": False,
        },
    ]

    payload = {
        "checkpoint": "P1757_S707",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "chain": "K_strict -> coefficients -> full non-skeleton L_total -> EOM -> nonproxy delivery sequence",
        "input_anchors": ["p1754_s704", "p1756_s706"],
        "consistency_precondition": consistency,
        "delivery_sequence": sequence,
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Zrealizować sekwencję D1..D5 bez skipów, następnie uruchomić D6 (nonproxy H1 4D) i dopiero po wyniku przejść do metric residual gate.",
        "lay_summary": "To jest plan wykonawczy krok po kroku. Najpierw trzeba dowieźć brakujące elementy techniczne, a dopiero potem uruchomić kluczowy test 4D.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
