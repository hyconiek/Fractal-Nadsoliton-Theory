#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1719 = GEN / "p1719_s669_strict_componentwise_sympy_runner_scaffold_checkpoint.json"
IN1721 = GEN / "p1721_s671_strict_componentwise_curvature_variation_pack_checkpoint.json"
OUT = GEN / "p1722_s672_strict_componentwise_runner_integration_step_checkpoint.json"


def main() -> None:
    p1719 = json.loads(IN1719.read_text(encoding="utf-8"))
    p1721 = json.loads(IN1721.read_text(encoding="utf-8"))

    integration = {
        "runner": p1719.get("runner_scaffold", {}).get("runner_name", "componentwise_sympy_metric_residual_runner_v1"),
        "variation_pack_loaded": True,
        "loaded_terms": list(p1721.get("componentwise_curvature_variation_pack", {}).keys()),
        "integration_stage": "SCHEMA_LEVEL_INTEGRATION_COMPLETE",
        "execution_stage": "RESIDUAL_COMPUTE_PENDING",
        "pending_blocks": [
            "componentwise_metric_ansatz_instantiation",
            "explicit_ELg_basis_assembly",
            "explicit_Emunu_basis_assembly",
            "residual_simplification_pass",
        ],
    }

    payload = {
        "checkpoint": "P1722_S672",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "chain": "K_strict -> coefficients -> full explicit L_total -> componentwise runner scaffold + curvature variation pack -> integration step",
        "runner_anchor": p1719.get("runner_scaffold", {}),
        "variation_pack_anchor": p1721.get("componentwise_curvature_variation_pack", {}),
        "integration_report": integration,
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "open_obligations": [
            "compute_first_metric_residual_expression",
            "classify_PASS_ZERO_or_OBSTRUCTION",
            "metric_sector_residual_zero_or_obstruction_certificate",
            "bianchi_ward_cross_consistency_certificate",
            "global_helmholtz_integrability_nonproxy",
            "brst_nilpotency_and_cohomology_nonproxy",
            "cutkosky_unitarity_full_sector",
            "counterterm_flow_renormalization_closure",
            "background_independence_family_theorem",
            "qw2191_selector_source_or_nonclosure_theorem",
        ],
        "next_honest_step": "Wykonać compute_first_metric_residual_expression i opublikować pierwszy jawny residual ELg-Emunu z klasyfikacją PASS_ZERO/OBSTRUCTION.",
        "lay_summary": "Zintegrowaliśmy brakujące wzory z narzędziem obliczeniowym. Następny krok to już właściwe policzenie reszty równania grawitacyjnego.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
