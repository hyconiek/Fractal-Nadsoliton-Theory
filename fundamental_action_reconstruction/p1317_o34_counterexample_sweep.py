"""P1317: O3.4 counterexample sweep over strict candidate classes C1-C4.

This script executes a reproducible strict-only counterexample sweep protocol
for the O3 pipeline (QW-2191 lane). It does not claim closure.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1317_o34_counterexample_sweep_report_v1.json"


@dataclass(frozen=True)
class Candidate:
    code: str
    branch_sign: int
    phase_class: float
    amplitude_class: float
    residual_slot: str
    admissible: bool


CANDIDATES = [
    Candidate("C1", +1, 0.41, 0.62, "open", True),
    Candidate("C2", +1, 0.43, 0.60, "open", True),
    Candidate("C3", +1, 0.40, 0.58, "open(Z2/eps)", True),
    Candidate("C4", +1, 0.42, 0.61, "open", True),
]

PHASE_TAU = 0.08
AMPLITUDE_TAU = 0.08


def pair_verdict(left: Candidate, right: Candidate) -> tuple[str, list[str]]:
    reasons: list[str] = []

    if not left.admissible or not right.admissible:
        reasons.append("inadmissible branch present")
        return "INADMISSIBLE_BRANCH", reasons

    if left.branch_sign != right.branch_sign:
        reasons.append("branch_sign mismatch")
        return "COUNTEREXAMPLE_FOUND", reasons

    phase_gap = abs(left.phase_class - right.phase_class)
    amp_gap = abs(left.amplitude_class - right.amplitude_class)

    if phase_gap > PHASE_TAU:
        reasons.append(f"phase gap {phase_gap:.3f} > tau")
        return "COUNTEREXAMPLE_FOUND", reasons

    if amp_gap > AMPLITUDE_TAU:
        reasons.append(f"amplitude gap {amp_gap:.3f} > tau")
        return "COUNTEREXAMPLE_FOUND", reasons

    if "open" in left.residual_slot or "open" in right.residual_slot:
        reasons.append("residual selector slot remains open")
        return "RESIDUAL_AMBIGUITY", reasons

    reasons.append("directionally equivalent under current tolerances")
    return "EQUIVALENT_DIRECTION", reasons


def main() -> None:
    pair_results = []
    verdict_counts = {
        "EQUIVALENT_DIRECTION": 0,
        "INADMISSIBLE_BRANCH": 0,
        "RESIDUAL_AMBIGUITY": 0,
        "COUNTEREXAMPLE_FOUND": 0,
    }

    for left, right in combinations(CANDIDATES, 2):
        verdict, reasons = pair_verdict(left, right)
        verdict_counts[verdict] += 1
        pair_results.append(
            {
                "pair": [left.code, right.code],
                "verdict": verdict,
                "reasons": reasons,
            }
        )

    o34_pass = verdict_counts["COUNTEREXAMPLE_FOUND"] == 0
    nonuniqueness_closed = o34_pass and verdict_counts["RESIDUAL_AMBIGUITY"] == 0

    payload = {
        "packet_id": "P1317_O3_4_COUNTEREXAMPLE_SWEEP_REPORT_V1",
        "date_utc": "2026-05-12",
        "strict_only": True,
        "input_candidates": [c.__dict__ for c in CANDIDATES],
        "tolerances": {"phase_tau": PHASE_TAU, "amplitude_tau": AMPLITUDE_TAU},
        "pair_results": pair_results,
        "verdict_counts": verdict_counts,
        "o3_4_status": "PASS_NO_COUNTEREXAMPLE" if o34_pass else "FAIL_COUNTEREXAMPLE_PRESENT",
        "nonuniqueness_residual_status": "CLOSED" if nonuniqueness_closed else "OPEN",
        "qw2191_strict_status": "NOT_CLOSED",
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
