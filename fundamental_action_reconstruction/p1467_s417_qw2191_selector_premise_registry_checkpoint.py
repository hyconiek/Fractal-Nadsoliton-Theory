#!/usr/bin/env python3
"""P1467 S4.17: register non-strict selector-premise candidates for QW-2191 follow-up."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
SUMMARY = GEN / "p1467_s417_qw2191_selector_premise_registry_summary.json"


def main() -> None:
    candidates = [
        {
            "id": "SP1_discrete_orientation_seed",
            "premise_type": "explicit symmetry-breaking seed",
            "local_testability": "HIGH",
            "strict_internal_source_available": False,
            "closure_status": "NON_STRICT_UNLESS_PROVEN_INTERNAL",
        },
        {
            "id": "SP2_entropy_weighted_selector",
            "premise_type": "entropy-weighted selector prior",
            "local_testability": "MEDIUM",
            "strict_internal_source_available": False,
            "closure_status": "NON_STRICT_UNLESS_PROVEN_INTERNAL",
        },
        {
            "id": "SP3_minimal_phase_anchor",
            "premise_type": "minimal phase-anchor selector",
            "local_testability": "HIGH",
            "strict_internal_source_available": False,
            "closure_status": "NON_STRICT_UNLESS_PROVEN_INTERNAL",
        },
    ]

    # Conservative professor-level choice: highest testability with explicit controllability.
    selected = "SP1_discrete_orientation_seed"

    summary = {
        "packet": "P1467",
        "status": "PASS_QW2191_PREMISE_REGISTRY_LOCAL_ONLY",
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "note": "No strict-core closure claim; registry only.",
        "candidates": candidates,
        "selected_next_local_test_candidate": selected,
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1467] status={summary['status']} selected={selected}")


if __name__ == "__main__":
    main()
