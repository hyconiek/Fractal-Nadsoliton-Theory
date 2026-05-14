#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
FILES = [
    "p1636_s586_strict_bidirectional_kernel_coefficient_closure_audit_summary.json",
    "p1637_s587_strict_identifiability_obstruction_from_kernel_to_coefficient_map_summary.json",
    "p1638_s588_strict_selector_constraint_candidate_for_nullspace_removal_summary.json",
    "p1639_s589_strict_full_chain_dossier_and_closure_blocker_map_summary.json",
    "p1640_s590_selector_internal_source_full_domain_candidate_export_summary.json",
    "p1641_s591_theorem_level_closure_requirement_matrix_summary.json",
]


def corrected_status(status: str, closure_open: bool) -> str:
    if closure_open and status.startswith("PASS_") and "LOCAL" not in status and "CANDIDATE" not in status:
        return f"KEEP_OPEN_RECLASSIFIED_FROM_{status}"
    return status


def main() -> None:
    audits = []
    for fn in FILES:
        p = GEN / fn
        data = json.loads(p.read_text(encoding="utf-8"))
        status = data.get("status", "UNKNOWN")
        closure = data.get("strict_core_closure", {}).get("status", "OPEN")
        open_flag = closure == "OPEN"
        corr = corrected_status(status, open_flag)
        audits.append({
            "file": fn,
            "original_status": status,
            "closure_status": closure,
            "corrected_status": corr,
            "changed": corr != status,
        })

    changed_n = sum(1 for a in audits if a["changed"])
    summary = {
        "checkpoint": "P1642_S592",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "PASS_P1642_NO_FALSE_PASS_AUDIT_EXPORTED",
        "audited_files": audits,
        "changes_count": changed_n,
        "policy": "No final-closure-like PASS labels when strict_core_closure is OPEN; keep LOCAL/CANDIDATE passes explicitly scoped.",
        "strict_only": True,
        "legacy_bridge_used": False,
        "next_honest_step": "Zastosować corrected_status labels bezpośrednio w upstream checkpoint exports (P1640+), aby wszystkie raporty były semantycznie konserwatywne.",
        "lay_summary": "Sprawdziliśmy etykiety wyników, żeby żadna nie sugerowała zamknięcia teorii bez pełnego dowodu. To kontrola jakości przeciw fałszywym PASS.",
    }

    out = GEN / "p1642_s592_no_false_pass_status_consistency_audit_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
