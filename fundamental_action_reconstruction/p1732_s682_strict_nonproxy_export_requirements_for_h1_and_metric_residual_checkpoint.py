#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1731 = GEN / "p1731_s681_strict_h1_gauge_scalar_nonproxy_execution_attempt_checkpoint.json"
IN1727 = GEN / "p1727_s677_strict_componentwise_metric_residual_coefficient_vector_stub_checkpoint.json"
OUT = GEN / "p1732_s682_strict_nonproxy_export_requirements_for_h1_and_metric_residual_checkpoint.json"


def main() -> None:
    p1731 = json.loads(IN1731.read_text(encoding="utf-8"))
    p1727 = json.loads(IN1727.read_text(encoding="utf-8"))

    missing_h1 = p1731.get("missing_required_exports", [])
    metric_basis = p1727.get("residual_coefficient_vector_stub", {}).get("basis_order", [])

    export_contract = {
        "h1_gauge_scalar": {
            "required_nonproxy_exports": missing_h1,
            "target_expression": "deltaE_A_mu_over_deltaH_minus_deltaE_H_over_deltaA_mu",
            "decision_policy": "PASS_ZERO_or_OBSTRUCTION_ONLY",
        },
        "metric_residual": {
            "required_nonproxy_exports": [
                "explicit_nonproxy_EL_g_munu",
                "explicit_nonproxy_E_munu",
                "componentwise_curvature_variation_terms_R2_Ric2_Riem2",
                "shared_background_family_contract",
                "index_and_sign_convention_lock",
            ],
            "basis_order": metric_basis,
            "target_expression": "EL_g_minus_E_munu_componentwise_basis_projection",
            "decision_policy": "PASS_ZERO_or_OBSTRUCTION_ONLY",
        },
    }

    payload = {
        "checkpoint": "P1732_S682",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "chain": "K_strict -> coefficients -> full nonskeleton L_total -> EOM -> nonproxy export requirements contract (H1 + metric residual)",
        "input_obstruction_anchor": p1731.get("execution_result", {}),
        "nonproxy_export_requirements_contract": export_contract,
        "full_lagrangian_density_nonskeleton_instantiated": p1731.get(
            "input_full_lagrangian_density_nonskeleton", {}
        ),
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Wyeksportować minimalny pakiet nonproxy z contractu P1732 i natychmiast uruchomić dwa obliczenia: (1) H1 cross-variation, (2) EL_g-E_munu na bazie B1/B2/B3/C1/C2 z decyzją PASS_ZERO albo OBSTRUCTION.",
        "lay_summary": "Po pierwszej blokadzie H1 mamy teraz dokładną checklistę brakujących eksportów. To pozwala przejść od ogólnego planu do konkretnej pracy implementacyjnej bez powtarzania tych samych prób.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
