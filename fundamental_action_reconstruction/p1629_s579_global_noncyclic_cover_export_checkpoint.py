#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1628 = GEN / "p1628_s578_globalization_program_for_strict_selector_proof_summary.json"


def main() -> None:
    s28 = json.loads(IN1628.read_text(encoding="utf-8"))

    cover = {
        "id": "C_global_noncyclic_cover",
        "charts": [
            {"name": "U0", "domain_d": [0.0, 1.5], "anchor": "near-core"},
            {"name": "U1", "domain_d": [1.0, 3.5], "anchor": "mid-band"},
            {"name": "U2", "domain_d": [3.0, 6.0], "anchor": "far-tail"},
        ],
        "overlap_graph": [
            ["U0", "U1"],
            ["U1", "U2"],
        ],
        "noncyclic": True,
    }

    summary = {
        "checkpoint": "P1629_S579",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "PASS_P1629_GLOBAL_NONCYCLIC_COVER_EXPORTED",
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "full_chain_used": s28["full_chain_used"],
        "global_noncyclic_cover": cover,
        "strict_core_closure": {
            "status": "OPEN",
            "missing_exports": ["E_selector_internal_source_full_domain", "E_full_variational_proof_log_machine_checkable"],
            "missing_witnesses": ["W_noncyclic_provider_for_selector_uniqueness_PROVED_GLOBAL"],
            "missing_theorems": ["T_qw2191_selector_uniqueness_strict_GLOBAL", "T_global_toe_closure_strict"],
        },
        "strict_only": True,
        "legacy_bridge_used": False,
        "external_team_validation_required": False,
        "next_honest_step": "Udowodnić L_overlap_compatibility na (U0∩U1) i (U1∩U2), używając tego samego funkcjonału selektora M.",
        "lay_summary": "Podzieliliśmy całą domenę na obszary i pokazaliśmy jak się łączą; to mapa potrzebna do globalnego dowodu, że teoria działa wszędzie, a nie tylko lokalnie.",
    }

    out = GEN / "p1629_s579_global_noncyclic_cover_export_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
