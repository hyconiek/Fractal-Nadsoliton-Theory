#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1758 = GEN / "p1758_s708_strict_nonproxy_e_a_mu_template_delivery_checkpoint.json"
IN1730 = GEN / "p1730_s680_strict_full_chain_physics_dossier_and_first_h1_witness_candidate_checkpoint.json"
OUT = GEN / "p1759_s709_strict_nonproxy_e_h_template_delivery_checkpoint.json"


def main() -> None:
    p1758 = json.loads(IN1758.read_text(encoding="utf-8"))
    p1730 = json.loads(IN1730.read_text(encoding="utf-8"))

    full_L = p1730.get("full_lagrangian_density_nonskeleton_instantiated", {})

    e_h_template = {
        "formal_name": "E_H_nonproxy_template_v1",
        "definition": "E_H := 1/sqrt(-g) * delta S / delta H^dagger",
        "sector_anchor": full_L.get("L_higgs", ""),
        "yukawa_anchor": full_L.get("L_yukawa", ""),
        "mix_anchor": full_L.get("L_mix", ""),
        "expected_structure": [
            "D_mu D^mu H",
            "mu_H^2 H + 2 lambda_H (H^dagger H) H",
            "xi_HR R H",
            "yukawa source terms",
            "phi-mix and curvature-mix corrections",
        ],
        "index_sign_lock_required": True,
        "background_family_contract_required": True,
        "status": "TEMPLATE_EXPORTED_NOT_YET_FULL_COMPONENTWISE",
    }

    payload = {
        "checkpoint": "P1759_S709",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "chain": "K_strict -> coefficients -> full non-skeleton L_total -> EOM -> D2 delivery (E_H template)",
        "input_anchors": ["p1758_s708", "p1730_s680"],
        "d2_delivery": e_h_template,
        "execution_sequence_status": {
            "D1_export_E_A_mu_nonproxy": p1758.get("execution_sequence_status", {}).get("D1_export_E_A_mu_nonproxy", "PARTIAL_TEMPLATE_DELIVERED"),
            "D2_export_E_H_nonproxy": "PARTIAL_TEMPLATE_DELIVERED",
            "D3_export_shared_background_family_contract": "PENDING",
            "D4_export_boundary_term_control_clause": "PENDING",
            "D5_finalize_boundary_control_contract": "PENDING",
            "D6_run_nonproxy_H1_4D": "BLOCKED",
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Dostarczyć D3: wspólny kontrakt rodziny teł, który łączy E_A^mu i E_H w tym samym setupie 4D.",
        "lay_summary": "Dostarczono drugi kluczowy szablon równania (dla Higgsa). Teraz oba pola pary H1 mają bazowy template, ale trzeba jeszcze ustalić wspólne tło obliczeń.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
