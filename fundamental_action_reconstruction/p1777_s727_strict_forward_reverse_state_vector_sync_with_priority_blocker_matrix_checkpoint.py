#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1766 = GEN / "p1766_s716_strict_forward_reverse_state_vector_update_with_bianchi_ward_gate_checkpoint.json"
IN1776 = GEN / "p1776_s726_strict_priority_gate_blocker_matrix_checkpoint.json"
OUT = GEN / "p1777_s727_strict_forward_reverse_state_vector_sync_with_priority_blocker_matrix_checkpoint.json"


def main() -> None:
    p1766 = json.loads(IN1766.read_text(encoding="utf-8"))
    p1776 = json.loads(IN1776.read_text(encoding="utf-8"))

    prev = p1766.get("updated_state_vector", {})
    blockers = p1776.get("priority_blocker_matrix", {})

    synced = {
        "forward_chain": {
            "K_strict": prev.get("forward_chain", {}).get("K_strict", "LOCKED_INPUT"),
            "coefficient_map": prev.get("forward_chain", {}).get("coefficient_map", "EXPORTED"),
            "full_nonskeleton_L_total": prev.get("forward_chain", {}).get("full_nonskeleton_L_total", "EXPORTED"),
            "covariant_E_A_mu_nonproxy": blockers.get("explicit_covariant_nonproxy_E_A_mu", {}).get("status", "OPEN"),
            "covariant_E_H_nonproxy": blockers.get("explicit_covariant_nonproxy_E_H", {}).get("status", "OPEN"),
            "covariant_EL_g_nonproxy": blockers.get("metric_EL_g_export", {}).get("status", "OPEN"),
            "boundary_term_control": blockers.get("boundary_term_control", {}).get("status", "OPEN"),
        },
        "reverse_chain": {
            "H1_cross_variation": blockers.get("H1_4D_weak_form_readiness", {}).get("status", "OPEN"),
            "metric_residual_EL_g_minus_E_munu": "OPEN_COMPONENTWISE_REQUIRED",
            "Bianchi_Ward_divergence_gate": blockers.get("Bianchi_Ward_consistency", {}).get("status", "OPEN"),
            "global_helmholtz_integrability": "OPEN",
        },
        "qg_theorem_gates": {
            "renormalization_counterterm_closure": "OPEN",
            "BRST_nilpotency": "BLOCKED_PENDING_BW_AND_W1",
            "Cutkosky_unitarity": "BLOCKED_PENDING_BW_AND_W1",
            "background_independence": "OPEN",
        },
    }

    payload = {
        "checkpoint": "P1777_S727",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "input_anchors": ["p1766_s716", "p1776_s726"],
        "previous_state_vector": prev,
        "priority_blocker_matrix_ref": blockers,
        "synced_state_vector": synced,
        "proof_update": {
            "proven": [
                "State-vector synchronization now reflects blocker-matrix statuses for E_A^mu/E_H/EL_g and BW/BRST/CUT dependencies.",
            ],
            "open": [
                "No new residual-zero theorem witness was produced in this synchronization step.",
            ],
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Zamknąć W1 do FULL_EXPORT i uruchomić komponentowe obliczenia H1 oraz EL_g-E_munu na tej samej rodzinie teł i konwencji indeksowej.",
        "lay_summary": "Ten krok porządkuje tablicę statusów: pokazuje dokładnie, co jest gotowe, a co nadal blokuje przejście do końcowych testów.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
