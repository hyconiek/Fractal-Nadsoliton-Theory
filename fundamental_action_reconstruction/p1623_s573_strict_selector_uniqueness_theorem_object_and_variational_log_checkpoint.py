#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1622 = GEN / "p1622_s572_full_strict_lagrangian_density_and_eom_summary.json"


def main() -> None:
    s22 = json.loads(IN1622.read_text(encoding="utf-8"))

    theorem_object = {
        "id": "T_qw2191_selector_uniqueness_strict",
        "premises": [
            "strict_only_pipeline",
            "noncyclic_provider_witness_required",
            "no_legacy_bridge_transfer",
        ],
        "claim": "selector_source is unique over full strict domain and induces deterministic closure branch",
        "proof_status": "OPEN",
        "blocking_obligations": [
            "construct_W_noncyclic_provider_for_selector_uniqueness",
            "export_E_selector_internal_source_full_domain",
        ],
    }

    variational_log_scaffold = {
        "id": "E_full_variational_proof_log_machine_checkable",
        "action": "S = ∫ d^4x sqrt(-g) L_total",
        "sectors": ["psi", "higgs", "gauge", "metric"],
        "checks": [
            "EL_psi_exact",
            "EL_higgs_exact",
            "EL_gauge_exact",
            "EL_metric_exact",
            "boundary_terms_cancellation",
        ],
        "proof_status": "OPEN",
    }

    summary = {
        "checkpoint": "P1623_S573",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "KEEP_OPEN_P1623_THEOREM_AND_VARIATIONAL_GAP",
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "input_chain": s22.get("full_lagrangian_density", {}),
        "input_eom": s22.get("euler_lagrange_eom", {}),
        "theorem_object": theorem_object,
        "variational_proof_log_scaffold": variational_log_scaffold,
        "strict_core_closure": {
            "status": "OPEN",
            "missing_exports": [
                "E_selector_internal_source_full_domain",
                "E_full_variational_proof_log_machine_checkable",
            ],
            "missing_witnesses": [
                "W_noncyclic_provider_for_selector_uniqueness",
            ],
            "missing_theorems": [
                "T_qw2191_selector_uniqueness_strict",
                "T_global_toe_closure_strict",
            ],
        },
        "strict_only": True,
        "legacy_bridge_used": False,
        "external_team_validation_required": False,
        "next_honest_step": "Zacząć konstrukcję W_noncyclic_provider_for_selector_uniqueness i podłączyć do formalnego dowodu T_qw2191.",
        "lay_summary": "Mamy już równania, ale potrzebujemy jeszcze formalnego dowodu, że teoria wybiera jedno, nieprzypadkowe rozwiązanie — to jest klucz do domknięcia.",
    }

    out = GEN / "p1623_s573_strict_selector_uniqueness_theorem_object_and_variational_log_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
