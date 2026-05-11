#!/usr/bin/env python3
"""P1268: strict-core SB1 ↔ QW-2191 compatibility theorem checkpoint."""
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
    parser.add_argument("--p1267", type=Path, default=GEN / "p1267_strict_core_sb1_h2_exclusion_theorem_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1268_strict_core_sb1_qw2191_compatibility_theorem_summary.json")
    args = parser.parse_args()

    p1266 = _read(args.p1266)
    p1267 = _read(args.p1267)

    if p1266.get("lane") != "STRICT_CORE_ONLY" or p1267.get("lane") != "STRICT_CORE_ONLY":
        raise SystemExit("P1268 requires STRICT_CORE_ONLY inputs from P1266 and P1267.")

    # Conservative theorem staging: compatibility is not yet fully discharged.
    compatibility_status = "PARTIAL_COMPATIBLE"

    out = {
        "packet": "P1268",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1266": str(args.p1266), "p1267": str(args.p1267)},
        "theorem": {
            "id": "SB1-QW2191",
            "statement": "SB1 selector-breaking source is compatible with QW-2191 under strict-core constraints.",
            "status": compatibility_status,
        },
        "supporting_points": [
            "H2 moved from OPEN to PARTIAL (P1267).",
            "QW2191 compatibility obligation remains explicit in SB1 matrix (P1266).",
            "No legacy-role import admitted in compatibility staging.",
        ],
        "remaining_obligations": [
            "Export strict-core lemma proving non-degeneracy survives full QW-2191 obstruction interface.",
            "Promote PARTIAL_COMPATIBLE to DISCHARGED with theorem-grade derivation chain.",
        ],
        "strict_kernel_closure_ready": False,
        "closure_policy": "STRICT_KERNEL_CLOSURE_FORBIDDEN_UNTIL_SB1_QW2191_DISCHARGED_AND_B1_OR_NB1",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1268] wrote {args.out}; status={compatibility_status}")


if __name__ == "__main__":
    main()
