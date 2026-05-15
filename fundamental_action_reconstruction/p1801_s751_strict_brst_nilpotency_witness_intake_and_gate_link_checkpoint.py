#!/usr/bin/env python3
from __future__ import annotations
import json
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT_POLICY = ROOT / "generated" / "p1801_s751_strict_brst_nilpotency_witness_intake_and_gate_link_checkpoint.json"
OUT_TEMPLATE = ROOT / "generated" / "p1801_s751_brst_nilpotency_witness_intake_template.json"

policy = {
    "packet_id": "P1801_S751",
    "as_of": str(date(2026, 5, 15)),
    "tg2_pass_requirements": [
        "G_BW_PASS_ZERO",
        "brst_charge_definition_present",
        "Q_squared_simplified_zero",
        "cohomology_subspace_check_present",
        "ghost_sector_consistency_present",
        "shared_freeze_id",
        "all_check_digests_present"
    ],
    "tg2_fail_status": "OPEN_LOCKED_BY_NILPOTENCY_WITNESS"
}

template = {
    "shared_freeze_id": "FREEZE_ID_REQUIRED",
    "dependencies": {"G_BW_status": "PASS_ZERO"},
    "brst_witness_pack": {
        "brst_charge_definition_id": "REQUIRED",
        "nilpotency_check": {
            "expression": "Q^2",
            "final_simplified": "0",
            "check_digest": "HASH_REQUIRED"
        },
        "cohomology_subspace_check_id": "REQUIRED",
        "ghost_sector_consistency_id": "REQUIRED",
        "check_digests": ["HASH1", "HASH2", "HASH3"]
    }
}

for p, d in ((OUT_POLICY, policy), (OUT_TEMPLATE, template)):
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", encoding="utf-8") as f:
        json.dump(d, f, indent=2, ensure_ascii=False)
        f.write("\n")
    print(p)
