#!/usr/bin/env python3
from __future__ import annotations
import json
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT_POLICY = ROOT / "generated" / "p1803_s753_strict_tg1_tg2_tg3_unified_global_closure_verdict_protocol_checkpoint.json"
OUT_TEMPLATE = ROOT / "generated" / "p1803_s753_global_closure_verdict_input_template.json"

policy = {
    "packet_id": "P1803_S753",
    "as_of": str(date(2026, 5, 15)),
    "global_pass_requirements": [
        "TG1_PASS",
        "TG2_PASS",
        "TG3_PASS",
        "Helmholtz_integrability_PASS",
        "Counterterm_renormalization_closure_PASS",
        "Background_independence_consistency_PASS"
    ],
    "fail_verdict": "OPEN_OBSTRUCTION_WITH_TRACE",
    "forbidden_claim_if_fail": "QG_SOLVED"
}

template = {
    "theorem_gates": {
        "TG1_BIANCHI_WARD_GLOBAL": "OPEN",
        "TG2_BRST_GLOBAL_NILPOTENCY": "OPEN",
        "TG3_CUTKOSKY_GLOBAL_UNITARITY": "OPEN"
    },
    "global_closures": {
        "global_helmholtz_integrability_status": "OPEN",
        "counterterm_renormalization_closure_status": "OPEN",
        "background_independence_consistency_status": "OPEN"
    },
    "missing_witnesses": [],
    "blocking_gate_chain": []
}

for p, d in ((OUT_POLICY, policy), (OUT_TEMPLATE, template)):
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", encoding="utf-8") as f:
        json.dump(d, f, indent=2, ensure_ascii=False)
        f.write("\n")
    print(p)
