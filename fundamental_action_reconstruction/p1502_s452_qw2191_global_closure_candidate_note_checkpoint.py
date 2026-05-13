#!/usr/bin/env python3
"""P1502 S4.52: publish global-closure candidate note with explicit caveats."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1500 = GEN / "p1500_s450_qw2191_selector_source_and_f_mapping_witness_summary.json"
P1501 = GEN / "p1501_s451_qw2191_adversarial_falsifier_sweep_summary.json"
P1499 = GEN / "p1499_s449_qw2191_global_closure_requirements_summary.json"

SUMMARY = GEN / "p1502_s452_qw2191_global_closure_candidate_note_summary.json"


def main() -> None:
    s1500 = json.loads(P1500.read_text(encoding="utf-8"))
    s1501 = json.loads(P1501.read_text(encoding="utf-8"))
    s1499 = json.loads(P1499.read_text(encoding="utf-8"))

    req = s1500["requirements_after_export"]
    req_all = all(bool(v) for v in req.values())
    no_falsifier = int(s1501["falsifier_count"]) == 0

    global_candidate_note_ready = req_all and no_falsifier

    summary = {
        "packet": "P1502",
        "status": "PASS_GLOBAL_CLOSURE_CANDIDATE_NOTE_LOCAL_ONLY" if global_candidate_note_ready else "FAIL_GLOBAL_CLOSURE_CANDIDATE_NOTE_LOCAL_ONLY",
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "closed_items": {
            "requirements_R1_to_R5": req,
            "adversarial_falsifier_absent": no_falsifier,
        },
        "not_claimed_yet": [
            "full_external_certification",
            "cross-lab independent replication completed",
            "global_toe_finalization",
        ],
        "global_closure_candidate_note": global_candidate_note_ready,
        "qw2191_closed": False,
        "next_step_recommendation": "S4.53: run independent external replication package and require signed reproducibility report before upgrading from candidate note to global closure release.",
        "layman_explanation": "To jest oficjalne 'prawie zamknięte': wszystko wewnętrznie działa i nic nie obaliło modelu w testach. Ale zanim ogłosimy pełne globalne domknięcie, ktoś niezależny musi to odtworzyć.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1502] status={summary['status']} candidate_note={global_candidate_note_ready}")


if __name__ == "__main__":
    main()
