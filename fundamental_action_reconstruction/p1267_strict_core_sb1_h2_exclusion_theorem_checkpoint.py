#!/usr/bin/env python3
"""P1267: strict-core H2 exclusion theorem checkpoint.

Goal: formalize that the selector-breaking step can be maintained without
non-strict axiom imports, under explicit evidence constraints.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1266", type=Path, default=GEN / "p1266_strict_core_sb1_hypothesis_discharge_matrix_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1267_strict_core_sb1_h2_exclusion_theorem_summary.json")
    args = parser.parse_args()

    p1266 = _read(args.p1266)
    rows = p1266.get("sb1_hypothesis_matrix", [])

    h2_row = next((r for r in rows if r.get("hypothesis") == "H2"), None)
    if not h2_row:
        raise SystemExit("P1267 requires H2 row in P1266 matrix.")

    theorem = {
        "id": "SB1-H2",
        "statement": (
            "Selector-breaking step is strict-core closed: no non-strict axiom import is required "
            "for admissible strict-core execution contexts."
        ),
        "status": "PARTIAL_DISCHARGE",
    }

    evidence = {
        "negative_control_artifact": "generated/p1231_w1_nonstrict_axiom_scenario_summary.json",
        "strict_lane_gate_artifact": "generated/p1258_strict_only_operational_lane_commitment_summary.json",
        "matrix_anchor": "generated/p1266_strict_core_sb1_hypothesis_discharge_matrix_summary.json",
    }

    conditions_open = [
        "Need formal strict-core derivation lemma replacing scenario-level exclusion evidence.",
        "Need explicit coupling to QW-2191 interface under strict-only selector step.",
    ]

    out = {
        "packet": "P1267",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1266": str(args.p1266)},
        "theorem": theorem,
        "evidence": evidence,
        "conditions_open": conditions_open,
        "h2_status_update": "PARTIAL",
        "strict_kernel_closure_ready": False,
        "closure_note": "H2 exclusion theorem partially discharged; strict-kernel closure still forbidden.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1267] wrote {args.out}")


if __name__ == "__main__":
    main()
