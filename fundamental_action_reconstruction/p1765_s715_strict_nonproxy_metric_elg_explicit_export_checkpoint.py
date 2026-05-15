#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1730 = GEN / "p1730_s680_strict_full_chain_physics_dossier_and_first_h1_witness_candidate_checkpoint.json"
IN1764 = GEN / "p1764_s714_strict_covariant_nonproxy_ea_eh_explicit_export_checkpoint.json"
OUT = GEN / "p1765_s715_strict_nonproxy_metric_elg_explicit_export_checkpoint.json"


def main() -> None:
    p1730 = json.loads(IN1730.read_text(encoding="utf-8"))
    p1764 = json.loads(IN1764.read_text(encoding="utf-8"))
    L = p1730.get("full_lagrangian_density_nonskeleton_instantiated", {})

    payload = {
        "checkpoint": "P1765_S715",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "chain": "K_strict -> coefficients -> full non-skeleton L_total -> explicit nonproxy metric EL_g export",
        "input_anchors": ["p1730_s680", "p1764_s714"],
        "explicit_nonproxy_EL_g_munu_v1": {
            "definition": "EL_g^{mu nu} := (2/sqrt(-g)) * delta S / delta g_{mu nu}",
            "covariant_structure": [
                "(M_P^2/2) G^{mu nu} + Lambda g^{mu nu}",
                "+ H_R2^{mu nu} + H_Ricci2^{mu nu} + H_Riemann2^{mu nu}",
                "- T_gauge^{mu nu} - T_fermion^{mu nu} - T_Higgs^{mu nu} - T_phi^{mu nu}",
                "- T_mix_RG^{mu nu} - T_CT^{mu nu}",
                "= 0",
            ],
            "sector_anchors": {
                "L_gravity": L.get("L_gravity", ""),
                "L_gauge": L.get("L_gauge", ""),
                "L_fermion": L.get("L_fermion", ""),
                "L_higgs": L.get("L_higgs", ""),
                "L_scalar_phi": L.get("L_scalar_phi", ""),
                "L_mix": L.get("L_mix", ""),
            },
            "residual_test_target": "EL_g_minus_E_munu_componentwise_on_B1_B2_B3_C1_C2",
            "classification": {
                "scope": "NONPROXY_COVARIANT_OPERATOR_LEVEL",
                "strict_local": "OPEN_COMPONENTWISE_REQUIRED",
                "global_theorem": "OPEN",
            },
        },
        "consistency_gates_update": {
            "boundary_term_control": "REUSED_FROM_P1762_FINALIZED_CONTRACT_TEMPLATE",
            "H1_4D_weak_form_readiness": p1764.get("readiness_update", {}).get("H1_4D_weak_form_readiness", "OPEN"),
            "Bianchi_Ward_consistency_gate": "OPEN_REQUIRES_COMPONENTWISE_DIVERGENCE_CHECK",
            "BRST_gate": "OPEN_REQUIRES_NILPOTENCY_WITNESS",
            "Cutkosky_gate": "OPEN_REQUIRES_UNITARITY_CUT_WITNESS",
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Na tej samej rodzinie teł wykonać komponentowy residual EL_g-E_munu (B1/B2/B3/C1/C2) oraz raport Bianchi/Ward divergence trace bez claimu PASS bez jawnego zera.",
        "lay_summary": "Mamy już jawny wzór metryczny na poziomie równań kowariantnych. Kolejny krok to policzyć go składnik po składniku i sprawdzić, czy wszystko się naprawdę zgadza.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
