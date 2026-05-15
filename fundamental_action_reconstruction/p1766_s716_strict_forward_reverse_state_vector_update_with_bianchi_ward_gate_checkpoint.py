#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1753 = GEN / "p1753_s703_strict_full_chain_forward_reverse_state_vector_checkpoint.json"
IN1764 = GEN / "p1764_s714_strict_covariant_nonproxy_ea_eh_explicit_export_checkpoint.json"
IN1765 = GEN / "p1765_s715_strict_nonproxy_metric_elg_explicit_export_checkpoint.json"
OUT = GEN / "p1766_s716_strict_forward_reverse_state_vector_update_with_bianchi_ward_gate_checkpoint.json"


def main() -> None:
    p1753 = json.loads(IN1753.read_text(encoding="utf-8"))
    p1764 = json.loads(IN1764.read_text(encoding="utf-8"))
    p1765 = json.loads(IN1765.read_text(encoding="utf-8"))

    prior = p1753.get("forward_reverse_state_vector", {})
    updated = {
        "forward_chain": {
            "K_strict": "LOCKED_INPUT",
            "coefficient_map": "EXPORTED",
            "full_nonskeleton_L_total": "EXPORTED",
            "covariant_E_A_mu_nonproxy": "EXPLICIT_OPERATOR_EXPORTED_P1764",
            "covariant_E_H_nonproxy": "EXPLICIT_OPERATOR_EXPORTED_P1764",
            "covariant_EL_g_nonproxy": "EXPLICIT_OPERATOR_EXPORTED_P1765",
        },
        "reverse_chain": {
            "H1_cross_variation": "OPEN_COMPONENTWISE_REQUIRED",
            "metric_residual_EL_g_minus_E_munu": "OPEN_COMPONENTWISE_REQUIRED",
            "Bianchi_Ward_divergence_gate": "OPEN_REQUIRES_COMPONENTWISE_DIVERGENCE_TRACE",
            "global_helmholtz_integrability": "OPEN",
        },
        "qg_theorem_gates": {
            "renormalization_counterterm_closure": "OPEN",
            "BRST_nilpotency": "OPEN",
            "Cutkosky_unitarity": "OPEN",
            "background_independence": "OPEN",
        },
    }

    payload = {
        "checkpoint": "P1766_S716",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "input_anchors": ["p1753_s703", "p1764_s714", "p1765_s715"],
        "prior_state_vector": prior,
        "updated_state_vector": updated,
        "gate_contract_update": {
            "H1_4D_weak_form_readiness": p1764.get("readiness_update", {}).get("H1_4D_weak_form_readiness", "OPEN"),
            "metric_ELg_export_status": p1765.get("explicit_nonproxy_EL_g_munu_v1", {}).get("classification", {}).get("scope", "OPEN"),
            "Bianchi_Ward_check_contract": {
                "target": "nabla_mu(E_total^{mu nu}) == 0 on declared background family",
                "input_required": [
                    "componentwise_EL_g_minus_E_munu_B1_B2_B3_C1_C2",
                    "shared_background_family_contract_v1",
                    "index_and_sign_lock",
                ],
                "admissible_outputs": ["PASS_ZERO", "OBSTRUCTION_WITH_DIVERGENCE_TRACE"],
            },
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Wykonać komponentowy EL_g-E_munu na B1/B2/B3/C1/C2 i natychmiast uruchomić divergence gate Bianchi/Ward z pełnym trace.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
