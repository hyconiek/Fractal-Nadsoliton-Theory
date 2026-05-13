#!/usr/bin/env python3
"""
P1412 SST-1 strict selector-source checkpoint.

Strict-only checkpoint toward F_nadsoliton => L_SM + L_GR under no-false-pass discipline.
"""

from __future__ import annotations

import json
from pathlib import Path


def main() -> None:
    base = Path(__file__).resolve().parent
    out_dir = base / "generated"
    out_dir.mkdir(parents=True, exist_ok=True)

    summary = {
        "checkpoint_id": "P1412-SST1",
        "as_of": "2026-05-13",
        "target": "F_nadsoliton => L_SM + L_GR",
        "mode": "strict_only_no_legacy_bridge",
        "strict_sources": [
            "K_strict_gate",
            "nadsoliton_ontology_lane",
            "alpha_geo_strict_derived_v1=4ln2",
            "strict_route_constraints_exported",
        ],
        "assumption_contract": {
            "new_axioms_allowed": False,
            "silent_legacy_role_transfer_allowed": False,
            "cyclic_l5_l12_regeneration_allowed": False,
            "requires_selector_source_inside_strict_core": True,
        },
        "selector_source_test": {
            "selector_map_v1_exported": False,
            "uniqueness_discharge": False,
            "transport_robustness_pass": False,
            "dual_replay_pass": False,
        },
        "verdict": "FAIL_STRICT_SELECTOR_SOURCE_MISSING",
        "fail_reasons": [
            "QW-2191 uniqueness obstruction still active on current strict-core state",
            "No exported strict-internal selector map satisfying no-new-axiom rule",
        ],
        "next_required_artifacts": [
            "strict_assumptions_used_v1.json",
            "selector_map_v1.json",
            "selector_obstruction_v1.json",
            "sst1_dual_replay_summary_v1.json",
        ],
        "status": "NO_FALSE_PASS",
    }

    out_path = out_dir / "p1412_sst1_strict_selector_source_summary.json"
    out_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"written": str(out_path), "verdict": summary["verdict"]}, indent=2))


if __name__ == "__main__":
    main()
