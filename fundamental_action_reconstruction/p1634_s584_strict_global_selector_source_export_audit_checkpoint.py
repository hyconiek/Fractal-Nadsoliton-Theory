#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1633 = GEN / "p1633_s583_strict_overlap_compatibility_to_global_selector_uniqueness_summary.json"
ATLAS = GEN / "selector_transition_global_c_v1_strict_v1.json"


def main() -> None:
    s33 = json.loads(IN1633.read_text(encoding="utf-8"))
    atlas = json.loads(ATLAS.read_text(encoding="utf-8"))

    blockers = []
    limits = atlas.get("hard_limits", [])
    for lim in limits:
        if "does not" in lim.lower() or "not infer" in lim.lower():
            blockers.append(lim)

    has_explicit_overlap_domains = "overlap_domains" in atlas
    has_projector_cocycle = "cocycle_discipline" in atlas
    lane_scoped = any("lane-scoped" in s.lower() for s in limits)

    export_ready = has_explicit_overlap_domains and has_projector_cocycle and (not lane_scoped) and (len(blockers) == 0)

    summary = {
        "checkpoint": "P1634_S584",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "KEEP_OPEN_P1634_SELECTOR_SOURCE_GLOBAL_EXPORT_INCOMPLETE",
        "route_target": s33["route_target"],
        "strict_chain": s33["strict_chain"],
        "selector_source_export_audit": {
            "atlas_source": str(ATLAS.relative_to(ROOT)),
            "has_explicit_overlap_domains": has_explicit_overlap_domains,
            "has_projector_cocycle_discipline": has_projector_cocycle,
            "lane_scoped_limit_present": lane_scoped,
            "blocking_statements": blockers,
            "export_ready": export_ready,
        },
        "strict_core_closure": {
            "status": "OPEN",
            "missing_exports": ["E_selector_internal_source_full_domain", "E_full_variational_proof_log_machine_checkable"],
            "missing_witnesses": ["W_noncyclic_provider_for_selector_uniqueness_PROVED_GLOBAL"],
            "missing_theorems": ["T_qw2191_selector_uniqueness_strict_GLOBAL", "T_global_toe_closure_strict"],
        },
        "strict_only": True,
        "legacy_bridge_used": False,
        "next_honest_step": "Wyeksportować nielane-scoped global selector atlas na pełnym strict domain z dowodem operator-level transition consistency wymaganym do E_selector_internal_source_full_domain.",
        "lay_summary": "Sprawdziliśmy istniejące obiekty sklejeń. Są wartościowe, ale nadal dotyczą ograniczonego zakresu, więc nie zamykają jeszcze globalnego źródła selektora.",
    }

    out = GEN / "p1634_s584_strict_global_selector_source_export_audit_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
