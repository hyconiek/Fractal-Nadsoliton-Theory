"""P1335: internal strict-side source construction attempt for A_branch_v1."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1335_internal_source_construction_for_a_branch_v1_report_v1.json"

PHASE_GRID = [0.25 + 0.005 * i for i in range(61)]
AMP_GRID = [0.45 + 0.005 * j for j in range(51)]
DELTA = 0.03
LAMBDA = 0.02


def score(phase: float, amp: float) -> float:
    return (phase - 0.40) + 0.9 * (amp - 0.58) - 0.6 * (phase - 0.40) ** 2


def dphase(phase: float) -> float:
    return 1.0 - 1.2 * (phase - 0.40)


def sign(x: float) -> int:
    return 1 if x >= 0 else -1


def strict_internal_selector(phase: float, amp: float) -> int:
    # candidate internal source: near-boundary tie-break from intrinsic slope
    s = score(phase, amp)
    if abs(s) >= DELTA:
        return sign(s)
    return sign(s + LAMBDA * dphase(phase))


def axiom_selector(phase: float, amp: float) -> int:
    return sign(score(phase, amp))


def main() -> None:
    total = 0
    boundary = 0
    mismatches = 0

    for p in PHASE_GRID:
        for a in AMP_GRID:
            total += 1
            si = strict_internal_selector(p, a)
            sa = axiom_selector(p, a)
            if abs(score(p, a)) < DELTA:
                boundary += 1
            if si != sa:
                mismatches += 1

    compatible = mismatches == 0

    payload = {
        "packet_id": "P1335_INTERNAL_SOURCE_CONSTRUCTION_FOR_A_BRANCH_V1_REPORT_V1",
        "date_utc": "2026-05-12",
        "delta": DELTA,
        "lambda": LAMBDA,
        "total_points": total,
        "boundary_points": boundary,
        "mismatches_vs_axiom": mismatches,
        "candidate_internal_source_compatible_with_axiom": compatible,
        "strict_internal_source_exportable": compatible,
        "status": "INTERNAL_SOURCE_CANDIDATE_COMPATIBLE" if compatible else "INTERNAL_SOURCE_CANDIDATE_INCOMPATIBLE",
        "qw2191_strict_status": "NOT_CLOSED",
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
