"""P1331: structural boundary resolution attempt via regularity constraint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1331_structural_boundary_resolution_attempt_report_v1.json"

PHASE_GRID = [0.25 + 0.005 * i for i in range(61)]
AMP_GRID = [0.45 + 0.005 * j for j in range(51)]
DELTA = 0.03
GMIN = 0.35


def score(phase: float, amp: float) -> float:
    return (phase - 0.40) + 0.9 * (amp - 0.58) - 0.6 * (phase - 0.40) ** 2


def grad_norm(phase: float) -> float:
    d_phase = 1.0 - 1.2 * (phase - 0.40)
    d_amp = 0.9
    return (d_phase * d_phase + d_amp * d_amp) ** 0.5


def main() -> None:
    total = 0
    admissible = 0
    boundary_in_adm = 0

    for p in PHASE_GRID:
        g = grad_norm(p)
        for a in AMP_GRID:
            total += 1
            if g < GMIN:
                continue
            admissible += 1
            if abs(score(p, a)) < DELTA:
                boundary_in_adm += 1

    boundary_ratio_adm = boundary_in_adm / admissible if admissible else 1.0
    structural_resolution_supported = boundary_ratio_adm <= 0.10

    payload = {
        "packet_id": "P1331_STRUCTURAL_BOUNDARY_RESOLUTION_ATTEMPT_REPORT_V1",
        "date_utc": "2026-05-12",
        "delta": DELTA,
        "gmin": GMIN,
        "total_points": total,
        "admissible_points": admissible,
        "boundary_points_in_admissible": boundary_in_adm,
        "boundary_ratio_in_admissible": boundary_ratio_adm,
        "structural_resolution_supported": structural_resolution_supported,
        "status": "STRUCTURAL_SUPPORT_WEAK" if not structural_resolution_supported else "STRUCTURAL_SUPPORT_STRONG",
        "global_l2_exportable": False,
        "qw2191_strict_status": "NOT_CLOSED",
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
