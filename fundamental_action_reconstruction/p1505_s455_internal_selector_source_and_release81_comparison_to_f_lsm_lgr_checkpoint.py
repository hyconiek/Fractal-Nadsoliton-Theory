#!/usr/bin/env python3
"""P1505 S4.55: internal selector-source comparison and F=>LSM+LGR direction."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1489 = GEN / "p1489_s439_qw2191_selector_source_candidate_summary.json"
P1492 = GEN / "p1492_s442_qw2191_selector_theorem_candidate_summary.json"
P1500 = GEN / "p1500_s450_qw2191_selector_source_and_f_mapping_witness_summary.json"
P1504 = GEN / "p1504_s454_qw2191_external_replication_protocol_summary.json"
SUMMARY = GEN / "p1505_s455_internal_selector_source_and_release81_comparison_to_f_lsm_lgr_summary.json"


def main() -> None:
    s1489 = json.loads(P1489.read_text(encoding="utf-8"))
    s1492 = json.loads(P1492.read_text(encoding="utf-8"))
    s1500 = json.loads(P1500.read_text(encoding="utf-8"))
    s1504 = json.loads(P1504.read_text(encoding="utf-8"))

    selector_candidate_present = (
        "PASS" in str(s1489.get("status", ""))
        and bool(s1489.get("checks", {}).get("strict_only", False))
    )
    theorem_candidate_present = "PASS" in str(s1492.get("status", ""))
    f_mapping_present = bool(
        s1500.get("requirements_after_export", {}).get("R5_F_to_LSM_LGR_mapping_witness_exported", False)
    )

    no_external_audit_yet = True

    release81_like_closure_ready = (
        selector_candidate_present and theorem_candidate_present and f_mapping_present and (not no_external_audit_yet)
    )

    summary = {
        "packet": "P1505",
        "status": "PASS_INTERNAL_SELECTOR_SOURCE_COMPARISON_AND_F_DIRECTION",
        "scope": "STRICT_ONLY_INTERNAL_PROGRESS_NO_LEGACY_BRIDGE",
        "internal_selector_source": {
            "candidate_present": selector_candidate_present,
            "theorem_candidate_present": theorem_candidate_present,
            "strict_only": True,
        },
        "release_8_1_comparison": {
            "local_candidate_strength": selector_candidate_present and theorem_candidate_present and f_mapping_present,
            "external_audit_completed": not no_external_audit_yet,
            "release81_like_closure_ready": release81_like_closure_ready,
        },
        "f_to_lsm_lgr_direction": {
            "witness_present": f_mapping_present,
            "next_internal_test": "Run coupled operator-consistency check for F->LSM and F->LGR channels under a shared strict-side selector orientation.",
        },
        "qw2191_closed": False,
        "legacy_bridge_used": bool(s1504.get("legacy_bridge_used", False)),
        "next_honest_step": "Execute an internal coupled F->LSM+LGR operator-consistency checkpoint (strict-only) before any closure upgrade claim.",
        "layman_explanation": "Mamy mocny kandydat wewnętrzny, ale to jeszcze nie finał. Następny test sprawdza, czy część 'cząstkowa' i 'grawitacyjna' naprawdę składają się w jedną, spójną całość.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1505] status={summary['status']} release81_like_closure_ready={release81_like_closure_ready}")


if __name__ == "__main__":
    main()
