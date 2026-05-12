"""P1332: curvature-aware boundary suppression attempt for residual loophole."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1332_curvature_aware_boundary_suppression_attempt_report_v1.json"

PHASE_GRID = [0.25 + 0.005 * i for i in range(61)]
AMP_GRID = [0.45 + 0.005 * j for j in range(51)]
DELTA = 0.03
H = 0.02


def score(phase: float, amp: float) -> float:
    return (phase - 0.40) + 0.9 * (amp - 0.58) - 0.6 * (phase - 0.40) ** 2


def curvature_correction() -> float:
    # |d2 score / d phase^2| = 1.2, zero in amp.
    return 0.5 * 1.2 * (H ** 2)


def main() -> None:
    total = 0
    admissible = 0
    boundary_in_adm = 0
    corr = curvature_correction()

    for p in PHASE_GRID:
        for a in AMP_GRID:
            total += 1
            s = abs(score(p, a)) - corr
            if s <= 0.0:
                continue
            admissible += 1
            if s < DELTA:
                boundary_in_adm += 1

    boundary_ratio = boundary_in_adm / admissible if admissible else 1.0
    suppression_supported = boundary_ratio <= 0.10

    payload = {
        "packet_id": "P1332_CURVATURE_AWARE_BOUNDARY_SUPPRESSION_ATTEMPT_REPORT_V1",
        "date_utc": "2026-05-12",
        "delta": DELTA,
        "h": H,
        "curvature_correction": corr,
        "total_points": total,
        "admissible_points": admissible,
        "boundary_points_in_admissible": boundary_in_adm,
        "boundary_ratio_in_admissible": boundary_ratio,
        "suppression_supported": suppression_supported,
        "status": "CURVATURE_SUPPORT_WEAK" if not suppression_supported else "CURVATURE_SUPPORT_STRONG",
        "global_l2_exportable": False,
        "qw2191_strict_status": "NOT_CLOSED",
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
