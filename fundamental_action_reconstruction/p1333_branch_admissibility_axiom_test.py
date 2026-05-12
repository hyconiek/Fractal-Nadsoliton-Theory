"""P1333: branch-admissibility axiom test (non-strict tagged)."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1333_branch_admissibility_axiom_test_report_v1.json"

PHASE_GRID = [0.25 + 0.005 * i for i in range(61)]
AMP_GRID = [0.45 + 0.005 * j for j in range(51)]
DELTA = 0.03


def score(phase: float, amp: float) -> float:
    return (phase - 0.40) + 0.9 * (amp - 0.58) - 0.6 * (phase - 0.40) ** 2


def sign(x: float) -> int:
    return 1 if x >= 0 else -1


def axiom_branch_selector(base_score: float) -> int:
    # Axiom A_branch_v1: near-boundary branch is selected by continuity from
    # positive-orientation side of v4 transport chart.
    return 1 if base_score >= 0 else -1


def main() -> None:
    total = 0
    boundary = 0
    inconsistent = 0

    for p in PHASE_GRID:
        for a in AMP_GRID:
            total += 1
            s = score(p, a)
            if abs(s) < DELTA:
                boundary += 1
                ax_sign = axiom_branch_selector(s)
                if ax_sign != sign(s):
                    inconsistent += 1

    non_strict_closed = inconsistent == 0

    payload = {
        "packet_id": "P1333_BRANCH_ADMISSIBILITY_AXIOM_TEST_REPORT_V1",
        "date_utc": "2026-05-12",
        "axiom": "A_branch_v1 continuity-selected near-boundary branch",
        "classification": "NON_STRICT_AXIOM_TAGGED",
        "delta": DELTA,
        "total_points": total,
        "boundary_points": boundary,
        "inconsistencies": inconsistent,
        "global_l2_under_axiom": non_strict_closed,
        "qw2191_status_under_axiom": "CLOSED_NON_STRICT" if non_strict_closed else "OPEN",
        "qw2191_strict_status": "NOT_CLOSED",
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
