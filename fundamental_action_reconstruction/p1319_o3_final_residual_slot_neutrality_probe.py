"""P1319: residual-slot neutrality/elimination probe for O3 final gate.

Goal: test whether open(Z2/eps) residual slot can flip direction while
remaining within strict admissibility windows.
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1319_o3_final_residual_slot_neutrality_probe_report_v1.json"

EPS_VALUES = [-1.0, -0.5, 0.0, 0.5, 1.0]
Z2_VALUES = [-1, 1]


def induced_branch_sign(base_sign: int, eps: float, z2: int) -> int:
    # conservative strict-side rule: residual slot can induce branch inversion
    # only when nonzero eps couples with negative z2 sector.
    if eps != 0.0 and z2 < 0:
        return -base_sign
    return base_sign


def main() -> None:
    base_sign = +1
    probes = []
    distinct_signs = set()

    for eps in EPS_VALUES:
        for z2 in Z2_VALUES:
            sign = induced_branch_sign(base_sign, eps, z2)
            distinct_signs.add(sign)
            probes.append(
                {
                    "eps": eps,
                    "z2": z2,
                    "branch_sign": sign,
                    "admissible": True,
                }
            )

    neutrality_holds = len(distinct_signs) == 1

    payload = {
        "packet_id": "P1319_O3_FINAL_RESIDUAL_SLOT_NEUTRALITY_PROBE_REPORT_V1",
        "date_utc": "2026-05-12",
        "strict_only": True,
        "model": "residual_slot_open_Z2_eps",
        "probes": probes,
        "distinct_branch_signs": sorted(distinct_signs),
        "neutrality_holds": neutrality_holds,
        "verdict": "NEUTRALITY_PROVEN" if neutrality_holds else "NEUTRALITY_NOT_PROVEN",
        "nonuniqueness_residual_status": "CLOSED" if neutrality_holds else "OPEN",
        "qw2191_strict_status": "CLOSED" if neutrality_holds else "NOT_CLOSED",
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
