#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1640 = GEN / "p1640_s590_selector_internal_source_full_domain_candidate_export_summary.json"
IN1639 = GEN / "p1639_s589_strict_full_chain_dossier_and_closure_blocker_map_summary.json"


def main() -> None:
    s40 = json.loads(IN1640.read_text(encoding="utf-8"))
    s39 = json.loads(IN1639.read_text(encoding="utf-8"))

    req = [
        {
            "id": "R1",
            "name": "E_selector_internal_source_full_domain_THEOREM_LEVEL",
            "status": "OPEN",
            "evidence_now": "candidate export assembled (P1640)",
            "gap": "missing theorem-level global operator consistency proof",
        },
        {
            "id": "R2",
            "name": "W_noncyclic_provider_for_selector_uniqueness_PROVED_GLOBAL",
            "status": "OPEN",
            "evidence_now": "noncyclic cover + local invertibility reports",
            "gap": "missing global witness proof across full overlap composition",
        },
        {
            "id": "R3",
            "name": "T_qw2191_selector_uniqueness_strict_GLOBAL",
            "status": "OPEN",
            "evidence_now": "local proxies + obstruction map + constraint candidate",
            "gap": "no global uniqueness theorem object exported",
        },
        {
            "id": "R4",
            "name": "E_full_variational_proof_log_machine_checkable",
            "status": "OPEN",
            "evidence_now": "local CAS appendix exists",
            "gap": "full-chain machine-checkable global variational log missing",
        },
    ]

    summary = {
        "checkpoint": "P1641_S591",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "KEEP_OPEN_P1641_REQUIREMENT_MATRIX_EXPORT",
        "route_target": s39["route_target"],
        "forward_chain": s39["forward_chain"],
        "closure_requirement_matrix": req,
        "strict_core_closure": s39["strict_core_closure"],
        "strict_only": True,
        "legacy_bridge_used": False,
        "next_honest_step": "Sformalizować i wyeksportować R1 jako theorem-level object; potem użyć go do R2 i R3.",
        "lay_summary": "Mamy już prawie wszystkie klocki techniczne, ale brakuje czterech formalnych certyfikatów. Ta macierz mówi dokładnie, co trzeba jeszcze dowieść, krok po kroku.",
    }

    out = GEN / "p1641_s591_theorem_level_closure_requirement_matrix_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
