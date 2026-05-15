#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1757 = GEN / "p1757_s707_strict_nonproxy_delivery_sequence_checkpoint.json"
IN1730 = GEN / "p1730_s680_strict_full_chain_physics_dossier_and_first_h1_witness_candidate_checkpoint.json"
OUT = GEN / "p1758_s708_strict_nonproxy_e_a_mu_template_delivery_checkpoint.json"


def main() -> None:
    p1757 = json.loads(IN1757.read_text(encoding="utf-8"))
    p1730 = json.loads(IN1730.read_text(encoding="utf-8"))

    full_L = p1730.get("full_lagrangian_density_nonskeleton_instantiated", {})

    e_a_mu_template = {
        "formal_name": "E_A_mu_nonproxy_template_v1",
        "definition": "E_A^mu := 1/sqrt(-g) * delta S / delta A_mu",
        "sector_anchor": full_L.get("L_gauge", ""),
        "mix_anchor": full_L.get("L_mix", ""),
        "expected_structure": [
            "nabla_nu( Z_a F_a^{nu mu} )",
            "J_matter^mu",
            "curvature-mix corrections from chi_RG and CT terms",
        ],
        "index_sign_lock_required": True,
        "background_family_contract_required": True,
        "status": "TEMPLATE_EXPORTED_NOT_YET_FULL_COMPONENTWISE",
    }

    payload = {
        "checkpoint": "P1758_S708",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "chain": "K_strict -> coefficients -> full non-skeleton L_total -> EOM -> D1 delivery (E_A^mu template)",
        "input_anchors": ["p1757_s707", "p1730_s680"],
        "d1_delivery": e_a_mu_template,
        "execution_sequence_status": {
            "D1_export_E_A_mu_nonproxy": "PARTIAL_TEMPLATE_DELIVERED",
            "D2_export_E_H_nonproxy": "PENDING",
            "D3_export_shared_background_family_contract": "PENDING",
            "D4_export_boundary_term_control_clause": "PENDING",
            "D5_finalize_boundary_control_contract": "PENDING",
            "D6_run_nonproxy_H1_4D": "BLOCKED",
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Dowieźć D2 (jawny E_H nonproxy) w tym samym standardzie i następnie zamknąć D3 kontraktem wspólnej rodziny teł.",
        "lay_summary": "Dostarczono pierwszy element techniczny do testu 4D: szablon równania dla pola cechowania. To jeszcze nie końcowy wzór obliczeniowy, ale konkretny krok naprzód.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
