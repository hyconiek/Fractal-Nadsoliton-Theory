#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1712 = GEN / "p1712_s662_strict_partial_sector_residual_zero_certificate_checkpoint.json"
IN1708 = GEN / "p1708_s658_strict_nonproxy_covariant_eom_first_explicit_formula_export_checkpoint.json"
OUT = GEN / "p1713_s663_strict_metric_residual_zero_attack_plan_checkpoint.json"


def main() -> None:
    p1712 = json.loads(IN1712.read_text(encoding="utf-8"))
    p1708 = json.loads(IN1708.read_text(encoding="utf-8"))

    metric_attack_plan = {
        "goal": "Produce first metric-sector EL residual-zero certificate in nonproxy-like convention.",
        "inputs": {
            "metric_template": p1708.get("first_explicit_nonproxy_formula_pack", {}).get("metric_eom_template", "MISSING"),
            "full_lagrangian_anchor": p1708.get("full_lagrangian_explicit_anchor", {}),
            "partial_sector_certificate": p1712.get("partial_sector_certificate", {}),
        },
        "execution_steps": [
            "freeze_metric_index_convention_and_signs",
            "instantiate_metric_variation_terms(delta_sqrtg, delta_R, delta_R2, delta_Ric2, delta_Riem2)",
            "assemble_Tmunu_from_matter_sectors_in_same_convention",
            "construct_E_munu_candidate",
            "compute_residual_ELg_minus_E_munu",
            "check_symbolic_zero_or_export_nonzero_obstruction",
        ],
        "success_condition": "EL_g(L_total) - E_munu == 0 in declared admissible class",
        "failure_condition": "export obstruction polynomial / unmatched tensor terms with provenance",
    }

    payload = {
        "checkpoint": "P1713_S663",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "chain": "K_strict -> coefficients -> full explicit L_total -> residual-zero partial certificate -> metric residual-zero attack plan",
        "current_partial_certificate": p1712.get("partial_sector_certificate", {}),
        "metric_residual_zero_attack_plan": metric_attack_plan,
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "open_obligations": [
            "metric_sector_residual_zero_nonproxy",
            "bianchi_ward_cross_consistency_certificate",
            "global_helmholtz_integrability_nonproxy",
            "brst_nilpotency_and_cohomology_nonproxy",
            "cutkosky_unitarity_full_sector",
            "counterterm_flow_renormalization_closure",
            "background_independence_family_theorem",
            "qw2191_selector_source_or_nonclosure_theorem",
        ],
        "next_honest_step": "Wykonać krok 1-4 planu i wyeksportować pierwszy jawny E_munu_candidate z rozpisanym T_munu oraz konwencją indeksową.",
        "lay_summary": "Mamy potwierdzone sektory cząstek, teraz atakujemy najtrudniejszy element: grawitację. Ten checkpoint ustawia konkretny plan, jak dojść do testu zero-reszty dla równań metrycznych.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
