#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1730 = GEN / "p1730_s680_strict_full_chain_physics_dossier_and_first_h1_witness_candidate_checkpoint.json"
IN1763 = GEN / "p1763_s713_strict_nonproxy_h1_4d_first_execution_attempt_checkpoint.json"
OUT = GEN / "p1764_s714_strict_covariant_nonproxy_ea_eh_explicit_export_checkpoint.json"


def main() -> None:
    p1730 = json.loads(IN1730.read_text(encoding="utf-8"))
    p1763 = json.loads(IN1763.read_text(encoding="utf-8"))
    L = p1730.get("full_lagrangian_density_nonskeleton_instantiated", {})

    export = {
        "E_A_mu_nonproxy_explicit_v2": {
            "definition": "E_A^mu := (1/sqrt(-g)) δS/δA_mu",
            "covariant_expression": [
                "nabla_nu( Z_a F_a^{nu mu} )",
                "+ 2*chi_RG*nabla_nu( R F_a^{nu mu} )",
                "+ J_matter,a^mu",
                "+ J_CT,a^mu",
                "= 0",
            ],
            "derived_from": ["L_gauge", "L_mix", "L_fermion", "L_higgs", "L_yukawa"],
            "sector_anchors": {
                "L_gauge": L.get("L_gauge", ""),
                "L_mix": L.get("L_mix", ""),
            },
            "nonproxy_scope": "covariant_4D_operator_level",
            "componentwise_state": "OPEN_NEEDS_INDEX_LEVEL_EXPANSION",
        },
        "E_H_nonproxy_explicit_v2": {
            "definition": "E_H := (1/sqrt(-g)) δS/δH^dagger",
            "covariant_expression": [
                "D_mu D^mu H",
                "+ mu_H^2 H + 2*lambda_H*(H^dagger H)H",
                "+ xi_HR * R * H",
                "+ lambda_phiH*phi^2*H",
                "+ d(CT_1loop)/dH^dagger",
                "- Y_u^dagger Q u - Y_d^dagger Q d - Y_e^dagger L e",
                "= 0",
            ],
            "derived_from": ["L_higgs", "L_yukawa", "L_scalar_phi", "L_mix"],
            "sector_anchors": {
                "L_higgs": L.get("L_higgs", ""),
                "L_yukawa": L.get("L_yukawa", ""),
                "L_scalar_phi": L.get("L_scalar_phi", ""),
            },
            "nonproxy_scope": "covariant_4D_operator_level",
            "componentwise_state": "OPEN_NEEDS_INDEX_LEVEL_EXPANSION",
        },
    }

    payload = {
        "checkpoint": "P1764_S714",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "chain": "K_strict -> coefficients -> full non-skeleton L_total -> explicit covariant E_A^mu/E_H nonproxy exports",
        "input_anchors": ["p1730_s680", "p1763_s713"],
        "blocking_context_from_p1763": p1763.get("execution_result", {}),
        "d1_d2_upgrade": export,
        "readiness_update": {
            "H1_4D_weak_form_readiness": "UPGRADED_TO_EXPLICIT_COVARIANT_OPERATOR_LEVEL",
            "H1_4D_strict_local_readiness": "OPEN_COMPONENTWISE_REQUIRED",
            "Bianchi_Ward_consistency_gate": "OPEN",
            "BRST_gate": "OPEN",
            "Cutkosky_gate": "OPEN",
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Wyeksportować indeksowo-komponentową postać E_A^mu i E_H na tej samej rodzinie teł i uruchomić H1: δE_A^mu/δH - δE_H/δA_mu.",
        "lay_summary": "Zrobiono ważny krok: równania cechowania i Higgsa są już jawnie zapisane w postaci kowariantnej (bez proxy), ale nadal trzeba je rozpisać komponentowo, by test H1 był końcowo rozstrzygający.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
