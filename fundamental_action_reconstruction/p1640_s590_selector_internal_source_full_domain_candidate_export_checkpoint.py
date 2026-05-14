#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1639 = GEN / "p1639_s589_strict_full_chain_dossier_and_closure_blocker_map_summary.json"
IN1631 = GEN / "p1631_s581_cover_wise_jacobian_invertibility_summary.json"
ATLAS = GEN / "selector_transition_global_c_v1_strict_v1.json"


def main() -> None:
    s39 = json.loads(IN1639.read_text(encoding="utf-8"))
    s31 = json.loads(IN1631.read_text(encoding="utf-8"))
    atlas = json.loads(ATLAS.read_text(encoding="utf-8"))

    candidate = {
        "export_name": "E_selector_internal_source_full_domain_candidate_v1",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "scope": "CANDIDATE_NON_THEOREM",
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": s39["route_target"],
        "kernel_to_eom_chain": s39["forward_chain"],
        "cover_local_invertibility": s31["cover_invertibility_report"],
        "atlas_transition_object": {
            "source": "generated/selector_transition_global_c_v1_strict_v1.json",
            "status": atlas.get("status", "unknown"),
            "hard_limits": atlas.get("hard_limits", []),
        },
        "closure_limits": {
            "not_theorem_level": True,
            "missing_global_operator_consistency_proof": True,
            "missing_qw2191_global_uniqueness_theorem": True,
        },
    }
    e_out = GEN / "e_selector_internal_source_full_domain_candidate_v1.json"
    e_out.write_text(json.dumps(candidate, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    summary = {
        "checkpoint": "P1640_S590",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "PASS_P1640_CANDIDATE_EXPORT_ASSEMBLED_KEEP_OPEN",
        "candidate_export": e_out.name,
        "strict_core_closure": {
            "status": "OPEN",
            "missing_exports": ["E_selector_internal_source_full_domain_THEOREM_LEVEL", "E_full_variational_proof_log_machine_checkable"],
            "missing_witnesses": ["W_noncyclic_provider_for_selector_uniqueness_PROVED_GLOBAL"],
            "missing_theorems": ["T_qw2191_selector_uniqueness_strict_GLOBAL", "T_global_toe_closure_strict"],
        },
        "next_honest_step": "Upgrade candidate export to theorem-level by proving global operator consistency and overlap composition on full strict domain.",
        "lay_summary": "Złożyliśmy brakujący pakiet w wersji kandydackiej. To jeszcze nie pełny dowód, ale konkretna baza do formalnego domknięcia.",
    }
    out = GEN / "p1640_s590_selector_internal_source_full_domain_candidate_export_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {e_out}\nWrote {out}")


if __name__ == "__main__":
    main()
