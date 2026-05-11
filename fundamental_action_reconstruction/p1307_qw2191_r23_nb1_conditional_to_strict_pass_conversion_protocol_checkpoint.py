#!/usr/bin/env python3
"""P1307: R23 NB1 conditional-to-strict pass conversion protocol checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1306", type=Path, default=GEN / "p1306_qw2191_r22_nb1_nontransfer_lemma_proof_attempt_batch_2_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1307_qw2191_r23_nb1_conditional_to_strict_pass_conversion_protocol_summary.json")
    args = parser.parse_args()

    p1306 = _read(args.p1306)
    if p1306.get("next_priority") != "R23_NB1_CONDITIONAL_TO_STRICT_PASS_CONVERSION_PROTOCOL":
        raise SystemExit("P1307 requires next_priority=R23_NB1_CONDITIONAL_TO_STRICT_PASS_CONVERSION_PROTOCOL from P1306.")
    if p1306.get("batch_status") != "NEAR_COMPLETE":
        raise SystemExit("P1307 requires NEAR_COMPLETE from P1306.")

    protocol = {
        "target_lemma": "LNB1_2",
        "required_exotic_classes": [
            "noncompact_transport_candidates",
            "piecewise_selector_transport_candidates",
            "axiom-tagged hybrid transport candidates",
        ],
        "strict_pass_criteria": [
            "each_exotic_class_violates_at_least_one_strict_invariant",
            "no_class_preserves_selector_uniqueness_guard_and_obstruction_consistency_together",
        ],
        "failure_exit": "DOWNGRADE_TO_BOUNDED_NONTRANSFER_CLAIM",
        "status": "PROTOCOL_READY",
    }

    out = {
        "packet": "P1307",
        "as_of": "2026-05-11",
        "lane": "NB1_NONBRIDGE_TRACK",
        "input": {"p1306": str(args.p1306)},
        "r23_conversion_protocol": protocol,
        "next_priority": "R24_NB1_EXOTIC_CLASS_SWEEP_AND_FINAL_LEMMA_STATUS",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1307] wrote {args.out}; status={protocol['status']}")


if __name__ == "__main__":
    main()
