#!/usr/bin/env python3
"""P1309: R25 NB1 formal closure statement and export packet checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1308", type=Path, default=GEN / "p1308_qw2191_r24_nb1_exotic_class_sweep_and_final_lemma_status_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1309_qw2191_r25_nb1_formal_closure_statement_and_export_packet_summary.json")
    args = parser.parse_args()

    p1308 = _read(args.p1308)
    if p1308.get("next_priority") != "R25_NB1_FORMAL_CLOSURE_STATEMENT_AND_EXPORT_PACKET":
        raise SystemExit("P1309 requires next_priority=R25_NB1_FORMAL_CLOSURE_STATEMENT_AND_EXPORT_PACKET from P1308.")
    if p1308.get("lnb1_2_final_status") != "PASS_STRICT":
        raise SystemExit("P1309 requires lnb1_2_final_status=PASS_STRICT from P1308.")

    closure = {
        "theorem_id": "NB1_QW2191_CLOSURE_R25",
        "scope": "NB1_NONBRIDGE_TRACK_ONLY",
        "statement": "QW-2191 admits non-transfer closure under NB assumptions with strict-role invariants preserved.",
        "non_exportables": [
            "legacy_to_strict_ontological_identity_claim",
            "global_strict_closure_outside_nb_scope",
        ],
        "status": "FORMAL_CLOSURE_DECLARED",
    }

    out = {
        "packet": "P1309",
        "as_of": "2026-05-11",
        "lane": "NB1_NONBRIDGE_TRACK",
        "input": {"p1308": str(args.p1308)},
        "r25_closure_statement": closure,
        "export_packet": {
            "allowed_claims": [
                "nb1_nontransfer_closure_within_scope",
                "strict_operational_predictions_under_nb_constraints",
            ],
            "disallowed_claims": closure["non_exportables"],
        },
        "next_priority": "R26_NB1_POSTCLOSURE_INDEPENDENT_REPLAY_AUDIT",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1309] wrote {args.out}; status={closure['status']}")


if __name__ == "__main__":
    main()
