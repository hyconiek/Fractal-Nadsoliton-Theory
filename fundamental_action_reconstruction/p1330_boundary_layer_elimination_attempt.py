"""P1330: boundary-layer elimination attempt for residual loophole (|score_v4|<delta)."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1330_boundary_layer_elimination_attempt_report_v1.json"

PHASE_GRID = [0.25 + 0.005 * i for i in range(61)]
AMP_GRID = [0.45 + 0.005 * j for j in range(51)]
DELTA = 0.03
NOISE_FLOOR = 0.012


def v4_score(phase: float, amp: float) -> float:
    return (phase - 0.40) + 0.9 * (amp - 0.58) - 0.6 * (phase - 0.40) ** 2


def main() -> None:
    total = 0
    boundary = 0
    safe = 0

    for p in PHASE_GRID:
        for a in AMP_GRID:
            total += 1
            s = abs(v4_score(p, a))
            if s < DELTA:
                boundary += 1
            if s >= DELTA + NOISE_FLOOR:
                safe += 1

    boundary_ratio = boundary / total
    safe_ratio = safe / total

    # Attempted elimination criterion (not yet theorem-level):
    # if boundary ratio small enough, mark candidate for exclusion-by-measure.
    exclusion_by_measure_supported = boundary_ratio <= 0.08

    payload = {
        "packet_id": "P1330_BOUNDARY_LAYER_ELIMINATION_ATTEMPT_REPORT_V1",
        "date_utc": "2026-05-12",
        "delta": DELTA,
        "noise_floor": NOISE_FLOOR,
        "total_points": total,
        "boundary_points": boundary,
        "boundary_ratio": boundary_ratio,
        "safe_points": safe,
        "safe_ratio": safe_ratio,
        "exclusion_by_measure_supported": exclusion_by_measure_supported,
        "global_l2_exportable": False,
        "status": "BOUNDARY_TOO_LARGE_FOR_EXCLUSION" if not exclusion_by_measure_supported else "BOUNDARY_MEASURE_SUPPORT",
        "qw2191_strict_status": "NOT_CLOSED",
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
