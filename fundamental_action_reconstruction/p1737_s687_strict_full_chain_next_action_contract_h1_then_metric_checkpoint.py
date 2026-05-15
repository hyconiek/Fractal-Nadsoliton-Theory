#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1736 = GEN / "p1736_s686_strict_kernel_to_eom_full_chain_and_reverse_qg_blocker_register_checkpoint.json"
OUT = GEN / "p1737_s687_strict_full_chain_next_action_contract_h1_then_metric_checkpoint.json"


def main() -> None:
    p1736 = json.loads(IN1736.read_text(encoding="utf-8"))
    blockers = p1736.get("qg_blocker_register", {})

    b1 = blockers.get("B1_H1_nonproxy_exports_missing", {})
    b2 = blockers.get("B2_metric_nonproxy_exports_missing", {})

    action_contract = {
        "step_1_h1_nonproxy_export": {
            "must_export": b1.get("missing", []),
            "must_compute": "deltaE_A_mu_over_deltaH_minus_deltaE_H_over_deltaA_mu",
            "decision": "PASS_ZERO_or_OBSTRUCTION_ONLY",
        },
        "step_2_metric_nonproxy_export": {
            "must_export": b2.get("missing", []),
            "must_compute": "EL_g_minus_E_munu_componentwise_on_B1_B2_B3_C1_C2",
            "decision": "PASS_ZERO_or_OBSTRUCTION_ONLY",
        },
        "step_3_qg_theorem_gate_updates": {
            "renormalization": "update_only_after_step_2",
            "unitarity": "update_only_after_step_2",
            "background_independence": "update_only_after_step_2",
        },
    }

    payload = {
        "checkpoint": "P1737_S687",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "chain": "K_strict -> coefficients -> full nonskeleton L_total -> EOM -> blocker register -> next-action contract",
        "full_lagrangian_density_nonskeleton_instantiated": p1736.get(
            "full_lagrangian_density_nonskeleton_instantiated", {}
        ),
        "next_action_contract": action_contract,
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Wykonać Step_1 (H1) zgodnie z contractem, natychmiast opublikować PASS_ZERO/OBSTRUCTION, potem wykonać Step_2 (EL_g-E_munu) i dopiero wtedy aktualizować bramki QG.",
        "lay_summary": "Kontrakt kolejnych działań mówi dokładnie: najpierw test H1, potem test grawitacyjny, a dopiero na końcu aktualizacja dużych twierdzeń kwantowej grawitacji.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
