#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_P1372 = GEN / "p1372_b1_assumption_register_summary.json"
OUT = GEN / "p1373_b1_theorem_or_obstruction_summary.json"


def main() -> None:
    p1372 = json.loads(IN_P1372.read_text(encoding="utf-8"))

    out = {
        "artifact": "p1373_b1_theorem_or_obstruction_summary",
        "as_of": "2026-05-12",
        "input_dependency": IN_P1372.name,
        "input_status": p1372.get("status"),
        "theorem_export": "NOT_YET",
        "obstruction_export": "YES",
        "b1_status": "OPEN_WITH_EXPLICIT_OBSTRUCTION",
        "missing_lemmas": [
            "L-B1-01: explicit nadsoliton->gauge operator/functional",
            "L-B1-02: Phi_B1 invariance under allowed scale/scheme transport",
            "L-B1-03: strict-core selector source or explicit symmetry-breaking satisfying QW-2191 discipline",
            "L-B1-04: term-by-term compatibility proof with no-silent-transfer table",
        ],
        "next_packet": "P1374_L_B1_01_OPERATOR_CONSTRUCTION_ATTEMPT",
    }

    GEN.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"[p1373] wrote {OUT}")


if __name__ == "__main__":
    main()
