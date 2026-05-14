#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1632 = GEN / "p1632_s582_full_strict_lagrangian_and_closure_obligation_summary.json"
IN1637 = GEN / "p1637_s587_strict_identifiability_obstruction_from_kernel_to_coefficient_map_summary.json"
IN1638 = GEN / "p1638_s588_strict_selector_constraint_candidate_for_nullspace_removal_summary.json"


def main() -> None:
    s32 = json.loads(IN1632.read_text(encoding="utf-8"))
    s37 = json.loads(IN1637.read_text(encoding="utf-8"))
    s38 = json.loads(IN1638.read_text(encoding="utf-8"))

    summary = {
        "checkpoint": "P1639_S589",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "KEEP_OPEN_P1639_GLOBAL_THEOREM_REQUIRED",
        "route_target": s32["route_target"],
        "forward_chain": {
            "kernel": "K_strict(d)=cos(omega*d+phi)/(1+beta*d^eta)",
            "kernel_params": s32["kernel_to_coeff"],
            "coefficients": s32["coefficients"],
            "full_lagrangian_density": s32["full_lagrangian_density"],
            "action": s32["action"],
            "eom_bundle": s32["eom_bundle"],
        },
        "backward_chain": {
            "identifiability_obstruction": s37["identifiability_analysis"],
            "selector_constraint_candidate": s38["constraint_candidate"],
            "current_readout": "Local reverse consistency exists, but global uniqueness is not proven.",
        },
        "strict_core_closure": {
            "status": "OPEN",
            "missing_exports": ["E_selector_internal_source_full_domain", "E_full_variational_proof_log_machine_checkable"],
            "missing_witnesses": ["W_noncyclic_provider_for_selector_uniqueness_PROVED_GLOBAL"],
            "missing_theorems": ["T_qw2191_selector_uniqueness_strict_GLOBAL", "T_global_toe_closure_strict"],
        },
        "strict_only": True,
        "legacy_bridge_used": False,
        "next_honest_step": "Wyeksportować theorem-level global selector law eliminujący nullspace na pełnym coverze i domykający E_selector_internal_source_full_domain.",
        "lay_summary": "Mamy kompletny przepis od kernela do równań ruchu. Problemem nie jest już brak równań, lecz brak globalnego dowodu, że odwrotna droga wybiera jedno, unikalne rozwiązanie w całej przestrzeni.",
    }

    out = GEN / "p1639_s589_strict_full_chain_dossier_and_closure_blocker_map_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
