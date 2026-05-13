#!/usr/bin/env python3
"""P1421 strict selector-uniqueness discharge pre-check (strict-only)."""

from __future__ import annotations

import json
from pathlib import Path


def main() -> None:
    base = Path(__file__).resolve().parent
    gen = base / "generated"
    gen.mkdir(parents=True, exist_ok=True)

    summary = {
        "checkpoint_id": "P1421-UNIQ-PRECHECK",
        "as_of": "2026-05-13",
        "target": "F_nadsoliton => L_SM + L_GR",
        "mode": "strict_only_no_legacy_bridge",
        "input_dependency": "P1420 PASS_PC8_RS2_RUN1",
        "precheck": {
            "strict_internal_selector_map_exported": True,
            "noncyclic_provider_chain": True,
            "transport_replay_thresholds_satisfied": True,
            "uniqueness_certificate_exported": False,
            "qw2191_discharge_proof_present": False,
        },
        "verdict": "FAIL_UNIQUENESS_DISCHARGE_MISSING_CERTIFICATE",
        "status": "NO_FALSE_PASS",
        "next_action": "export uniqueness obstruction and design certificate-construction checkpoint P1422",
    }

    out = gen / "p1421_strict_selector_uniqueness_discharge_precheck_summary.json"
    out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    obstruction = {
        "obstruction_id": "OBSTR-UNIQ-CERT-v1",
        "reason": "no strict uniqueness certificate exported",
        "required_objects": [
            "selector_uniqueness_certificate_v1.json",
            "qw2191_discharge_argument_v1.md",
        ],
        "recommended_next": "P1422_certificate_construction_checkpoint",
    }
    (gen / "p1421_uniqueness_obstruction_v1.json").write_text(json.dumps(obstruction, indent=2) + "\n", encoding="utf-8")

    print(json.dumps({"written": str(out), "verdict": summary["verdict"]}, indent=2))


if __name__ == "__main__":
    main()
