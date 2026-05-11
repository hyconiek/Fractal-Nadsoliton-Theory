#!/usr/bin/env python3
"""P1278: obstruction-resolution design checkpoint after failed QW-2191 closure attempt."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1277", type=Path, default=GEN / "p1277_qw2191_full_closure_theorem_attempt_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1278_qw2191_obstruction_resolution_design_summary.json")
    args = parser.parse_args()

    p1277 = _read(args.p1277)
    if p1277.get("qw2191_closure_status") != "NOT_CLOSED":
        raise SystemExit("P1278 is only valid after an unresolved QW-2191 closure attempt.")

    design = {
        "goal": "Resolve primary QW-2191 blocker via minimal lemma set",
        "lemmas": [
            {"id": "R1", "task": "Prove full-neighborhood completeness of obstruction domain declaration.", "priority": "P1"},
            {"id": "R2", "task": "Derive strict-core-only bound transport across neighborhood boundary.", "priority": "P1"},
            {"id": "R3", "task": "Show selector non-degeneracy stability under transported bound.", "priority": "P2"},
        ],
    }

    out = {
        "packet": "P1278",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1277": str(args.p1277)},
        "resolution_design": design,
        "execution_mode": "THEOREM_REPAIR_PROGRAM",
        "strict_kernel_closure_ready": False,
        "global_closure_status": "OPEN",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1278] wrote {args.out}; lemmas={len(design['lemmas'])}")


if __name__ == "__main__":
    main()
