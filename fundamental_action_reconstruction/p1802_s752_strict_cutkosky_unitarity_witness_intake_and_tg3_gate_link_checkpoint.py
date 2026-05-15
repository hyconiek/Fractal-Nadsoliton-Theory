#!/usr/bin/env python3
from __future__ import annotations
import json
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT_POLICY = ROOT / "generated" / "p1802_s752_strict_cutkosky_unitarity_witness_intake_and_tg3_gate_link_checkpoint.json"
OUT_TEMPLATE = ROOT / "generated" / "p1802_s752_cutkosky_unitarity_witness_intake_template.json"

policy = {
    "packet_id": "P1802_S752",
    "as_of": str(date(2026, 5, 15)),
    "tg3_pass_requirements": [
        "G_BW_PASS_ZERO",
        "TG2_BRST_PASS",
        "cut_discontinuity_witness_present",
        "optical_theorem_match_true",
        "ghost_pole_exclusion_true",
        "background_family_consistency_present",
        "shared_freeze_id",
        "all_check_digests_present"
    ],
    "tg3_fail_status": "OPEN_LOCKED_BY_UNITARITY_WITNESS"
}

template = {
    "shared_freeze_id": "FREEZE_ID_REQUIRED",
    "dependencies": {
        "G_BW_status": "PASS_ZERO",
        "TG2_BRST_GLOBAL_NILPOTENCY": "PASS"
    },
    "cutkosky_witness_pack": {
        "cut_discontinuity_witness_id": "REQUIRED",
        "optical_theorem_match": {"status": "TRUE", "check_digest": "HASH_REQUIRED"},
        "ghost_pole_exclusion": {"status": "TRUE", "check_digest": "HASH_REQUIRED"},
        "background_family_consistency_id": "REQUIRED",
        "check_digests": ["HASH1", "HASH2", "HASH3"]
    }
}

for p, d in ((OUT_POLICY, policy), (OUT_TEMPLATE, template)):
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", encoding="utf-8") as f:
        json.dump(d, f, indent=2, ensure_ascii=False)
        f.write("\n")
    print(p)
