"""P1320: premise-augmented selector-law probe for QW-2191.

This probe evaluates a strict-lane candidate with an explicit symmetry-breaking
premise P_sel_v1: admissible branch requires z2 = +1.

Important: this is a premise-augmented (non-strict) closure test.
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1320_premise_augmented_selector_law_probe_report_v1.json"

EPS_VALUES = [-1.0, -0.5, 0.0, 0.5, 1.0]
Z2_VALUES = [-1, 1]


def admissible_under_p_sel_v1(z2: int) -> bool:
    return z2 == 1


def induced_branch_sign(base_sign: int, eps: float, z2: int) -> int:
    if eps != 0.0 and z2 < 0:
        return -base_sign
    return base_sign


def main() -> None:
    base_sign = +1
    probes = []
    admissible_signs = set()

    for eps in EPS_VALUES:
        for z2 in Z2_VALUES:
            admissible = admissible_under_p_sel_v1(z2)
            sign = induced_branch_sign(base_sign, eps, z2)
            if admissible:
                admissible_signs.add(sign)
            probes.append({"eps": eps, "z2": z2, "admissible": admissible, "branch_sign": sign})

    uniqueness_under_premise = len(admissible_signs) == 1

    payload = {
        "packet_id": "P1320_PREMISE_AUGMENTED_SELECTOR_LAW_PROBE_REPORT_V1",
        "date_utc": "2026-05-12",
        "premise": "P_sel_v1: z2 must equal +1 for admissible selector branch",
        "strict_core": False,
        "classification": "NON_STRICT_PREMISE_AUGMENTED",
        "probes": probes,
        "admissible_branch_signs": sorted(admissible_signs),
        "uniqueness_under_premise": uniqueness_under_premise,
        "qw2191_status_under_premise": "CLOSED_NON_STRICT" if uniqueness_under_premise else "OPEN",
        "qw2191_strict_status": "NOT_CLOSED",
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
