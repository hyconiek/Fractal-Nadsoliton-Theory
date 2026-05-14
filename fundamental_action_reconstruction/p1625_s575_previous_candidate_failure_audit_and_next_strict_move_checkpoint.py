#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

INPUTS = {
    "P1583": "p1583_s533_formal_proof_object_and_global_stability_composition_summary.json",
    "P1585": "p1585_s535_formal_l1_l3_theorem_objects_and_toe_composition_summary.json",
    "P1596": "p1596_s546_selector_uniqueness_bridge_theorem_object_summary.json",
    "P1605": "p1605_s555_np1_provider_instantiation_and_replay_summary.json",
    "P1623": "p1623_s573_strict_selector_uniqueness_theorem_object_and_variational_log_summary.json",
    "P1624": "p1624_s574_noncyclic_selector_witness_from_strict_kernel_summary.json",
}


def load(name: str) -> dict:
    return json.loads((GEN / INPUTS[name]).read_text(encoding="utf-8"))


def main() -> None:
    docs = {k: load(k) for k in INPUTS}

    prior_failures = {
        k: v.get("status", "UNKNOWN")
        for k, v in docs.items()
        if any(tok in v.get("status", "") for tok in ["FAIL", "KEEP_OPEN", "NOT_READY", "INCOMPLETE"])
    }

    candidate_verdict = {
        "previous_candidates_closed": False,
        "reason": "previous candidates remain FAIL/KEEP_OPEN/INCOMPLETE in strict selector-uniqueness path",
        "evidence": prior_failures,
    }

    s24 = docs["P1624"]
    chain_anchor = {
        "kernel": "K_strict(omega,phi,beta,eta)",
        "coefficients": s24["witness_object"]["coeff_values"],
        "lagrangian": "L_total = L_strict_scalar + L_SM + L_GR + L_mix",
        "eom": ["psi", "higgs", "gauge", "metric"],
    }

    summary = {
        "checkpoint": "P1625_S575",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "KEEP_OPEN_P1625_PREVIOUS_CANDIDATES_NOT_CLOSED",
        "candidate_audit": candidate_verdict,
        "strict_chain_anchor": chain_anchor,
        "strict_core_closure": {
            "status": "OPEN",
            "missing_exports": ["E_selector_internal_source_full_domain", "E_full_variational_proof_log_machine_checkable"],
            "missing_witnesses": ["W_noncyclic_provider_for_selector_uniqueness_PROVED"],
            "missing_theorems": ["T_qw2191_selector_uniqueness_strict", "T_global_toe_closure_strict"],
        },
        "strict_only": True,
        "legacy_bridge_used": False,
        "external_team_validation_required": False,
        "next_honest_step": "Sformalizować dowód, że selection_rule z P1624 daje jednoznaczny selector_source na full-domain i podpiąć do T_qw2191.",
        "lay_summary": "Sprawdziliśmy stare podejścia: nie domknęły teorii, więc nie powtarzamy ich, tylko budujemy brakujący formalny dowód na bazie obecnego pełnego łańcucha równań.",
    }

    out = GEN / "p1625_s575_previous_candidate_failure_audit_and_next_strict_move_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
