#!/usr/bin/env python3
"""P1503 S4.53: strict-physical next honest step plan for QW-2191."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1500 = GEN / "p1500_s450_qw2191_selector_source_and_f_mapping_witness_summary.json"
P1501 = GEN / "p1501_s451_qw2191_adversarial_falsifier_sweep_summary.json"
P1502 = GEN / "p1502_s452_qw2191_global_closure_candidate_note_summary.json"
SUMMARY = GEN / "p1503_s453_qw2191_strict_physical_closure_next_honest_step_summary.json"


def main() -> None:
    s1500 = json.loads(P1500.read_text(encoding="utf-8"))
    s1501 = json.loads(P1501.read_text(encoding="utf-8"))
    s1502 = json.loads(P1502.read_text(encoding="utf-8"))

    req_ok = all(bool(v) for v in s1500["requirements_after_export"].values())
    no_falsifier_local = int(s1501["falsifier_count"]) == 0
    candidate_note = bool(s1502.get("global_closure_candidate_note", False))

    ready_for_external_replication = req_ok and no_falsifier_local and candidate_note

    summary = {
        "packet": "P1503",
        "status": "PASS_STRICT_PHYSICAL_NEXT_HONEST_STEP_DEFINED" if ready_for_external_replication else "FAIL_STRICT_PHYSICAL_NEXT_HONEST_STEP_DEFINED",
        "scope": "STRICT_ONLY_LOCAL_TO_EXTERNAL_REPLICATION_GATE",
        "inputs": {
            "p1500_requirements_after_export_all_true": req_ok,
            "p1501_local_falsifier_absent": no_falsifier_local,
            "p1502_candidate_note_published": candidate_note,
        },
        "qw2191_closed": False,
        "legacy_bridge_used": False,
        "next_honest_step": "S4.53 external independent replication with signed reproducibility report before any closure upgrade.",
        "strict_physical_basis": [
            "K_strict_gate operational strict kernel",
            "nadsoliton ontology as primordial information",
            "alpha_geo_strict_derived_v1 = 4 ln 2 (strict-side source-upgrade object)",
        ],
        "layman_explanation": "W skrócie: lokalnie wszystko wygląda dobrze, ale naukowo to wciąż za mało. Teraz ktoś niezależny musi powtórzyć wynik i podpisać raport, zanim wolno ogłosić pełne domknięcie.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1503] status={summary['status']} ready_for_external_replication={ready_for_external_replication}")


if __name__ == "__main__":
    main()
