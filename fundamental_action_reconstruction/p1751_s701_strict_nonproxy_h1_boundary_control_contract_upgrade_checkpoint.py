#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1732 = GEN / "p1732_s682_strict_nonproxy_export_requirements_for_h1_and_metric_residual_checkpoint.json"
IN1750 = GEN / "p1750_s700_strict_reduced_h1_gauge_higgs_boundary_term_audit_checkpoint.json"
OUT = GEN / "p1751_s701_strict_nonproxy_h1_boundary_control_contract_upgrade_checkpoint.json"


def main() -> None:
    p1732 = json.loads(IN1732.read_text(encoding="utf-8"))
    p1750 = json.loads(IN1750.read_text(encoding="utf-8"))

    old_contract = p1732.get("nonproxy_export_requirements_contract", {})
    h1_old = old_contract.get("h1_gauge_scalar", {})

    boundary_audit = p1750.get("boundary_term_audit", {})
    boundary_exact = boundary_audit.get("exact_derivative_confirmed", False)

    h1_upgraded = {
        "required_nonproxy_exports": h1_old.get("required_nonproxy_exports", []),
        "target_expression": "deltaE_A_mu_over_deltaH_minus_deltaE_H_over_deltaA_mu",
        "decision_policy": "PASS_ZERO_or_OBSTRUCTION_ONLY",
        "boundary_control_contract": {
            "required": True,
            "reason": "P1750 boundary-sensitive obstruction requires explicit boundary clause handling in nonproxy 4D run.",
            "boundary_clause_template": "Integral boundary contribution must be explicitly exported and set by declared boundary family contract before PASS decisions.",
            "audit_anchor": "p1750_s700",
            "reduced_exact_derivative_confirmed": boundary_exact,
        },
        "verdict_schema": {
            "strict_local": ["PASS_ZERO", "OBSTRUCTION"],
            "weak_form_with_boundary": ["PASS_WEAK_FORM_WITH_BOUNDARY_CLAUSE", "OBSTRUCTION"],
            "promotion_rule": "weak-form pass alone cannot promote strict-core closure",
        },
    }

    upgraded_contract = dict(old_contract)
    upgraded_contract["h1_gauge_scalar"] = h1_upgraded

    payload = {
        "checkpoint": "P1751_S701",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "chain": "K_strict -> coefficients -> full nonskeleton L_total -> EOM -> H1 boundary-aware nonproxy contract upgrade",
        "input_contract_anchor": "p1732_s682",
        "input_boundary_audit_anchor": "p1750_s700",
        "nonproxy_export_requirements_contract_upgraded": upgraded_contract,
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Wykonać pierwszy 4D nonproxy run H1(A_mu,H) zgodnie z kontraktem boundary-aware i zwrócić dualny raport: strict_local + weak_form_with_boundary.",
        "lay_summary": "Po audycie brzegowym doprecyzowaliśmy zasady testu H1: zanim ogłosimy sukces, trzeba jawnie pokazać także kontrolę warunków brzegowych. To zwiększa rygor i usuwa ryzyko fałszywego passa.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
