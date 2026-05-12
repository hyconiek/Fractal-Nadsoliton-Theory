"""P1329: residual-loophole elimination proof attempt (case split z2/eps)."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1329_residual_loophole_elimination_proof_attempt_report_v1.json"

PHASE_GRID = [0.25 + 0.01 * i for i in range(31)]
AMP_GRID = [0.45 + 0.01 * j for j in range(26)]
EPS_VALUES = [-1.0, -0.5, 0.0, 0.5, 1.0]
Z2_VALUES = [-1, 1]
MARGIN_DELTA = 0.03


def v4_score(phase: float, amp: float) -> float:
    return (phase - 0.40) + 0.9 * (amp - 0.58) - 0.6 * (phase - 0.40) ** 2


def branch_from_score(score: float) -> int:
    return 1 if score >= 0.0 else -1


def residual_shift(eps: float, z2: int) -> float:
    return 0.015 * eps * z2


def main() -> None:
    total = 0
    admissible_points = 0
    violations = 0

    for phase in PHASE_GRID:
        for amp in AMP_GRID:
            base = v4_score(phase, amp)
            base_sign = branch_from_score(base)
            total += 1

            if abs(base) < MARGIN_DELTA:
                continue

            admissible_points += 1
            for eps in EPS_VALUES:
                for z2 in Z2_VALUES:
                    shifted = base + residual_shift(eps, z2)
                    sign2 = branch_from_score(shifted)
                    if sign2 != base_sign:
                        violations += 1

    l2_proven = admissible_points > 0 and violations == 0

    payload = {
        "packet_id": "P1329_RESIDUAL_LOOPHOLE_ELIMINATION_PROOF_ATTEMPT_REPORT_V1",
        "date_utc": "2026-05-12",
        "grid_size": {"phase": len(PHASE_GRID), "amp": len(AMP_GRID)},
        "margin_delta": MARGIN_DELTA,
        "total_points": total,
        "admissible_points": admissible_points,
        "violations": violations,
        "l2_proven_on_admissible_margin_domain": l2_proven,
        "global_l2_exportable": False,
        "status": "PARTIAL_CONDITIONAL_SUPPORT" if l2_proven else "L2_NOT_PROVEN",
        "qw2191_strict_status": "NOT_CLOSED",
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
