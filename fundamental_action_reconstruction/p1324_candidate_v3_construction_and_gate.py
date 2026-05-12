"""P1324: candidate-v3 construction with robustness + informativeness gate."""

from __future__ import annotations

import json
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1324_candidate_v3_construction_and_gate_report_v1.json"

N_SAMPLES = 400
SEED = 1324


def s_sel_strict_v3(phase: float, amp: float) -> int:
    # conservative robust version: non-negative affine floor on current window
    score = 0.20 + (phase - 0.25) + 0.40 * (amp - 0.45)
    return 1 if score >= 0.0 else -1


def generate_sample(rng: random.Random) -> tuple[float, float]:
    return rng.uniform(0.25, 0.55), rng.uniform(0.45, 0.70)


def main() -> None:
    rng = random.Random(SEED)
    negatives = []
    positives = 0

    for i in range(N_SAMPLES):
        phase, amp = generate_sample(rng)
        sign = s_sel_strict_v3(phase, amp)
        if sign < 0:
            negatives.append({"id": i, "phase": phase, "amp": amp})
        else:
            positives += 1

    robust = len(negatives) == 0
    sign_diversity = 2 if positives and negatives else 1
    informative = sign_diversity == 2

    admissible_as_selector = robust and informative

    payload = {
        "packet_id": "P1324_CANDIDATE_V3_CONSTRUCTION_AND_GATE_REPORT_V1",
        "date_utc": "2026-05-12",
        "model": "S_sel_strict_v3",
        "formula": "sign(0.20 + (phase-0.25) + 0.40*(amp-0.45))",
        "sample_count": N_SAMPLES,
        "seed": SEED,
        "positive_count": positives,
        "negative_count": len(negatives),
        "counterexamples": negatives[:20],
        "robust_on_window": robust,
        "sign_diversity": sign_diversity,
        "informative_selector": informative,
        "admissible_as_strict_selector_candidate": admissible_as_selector,
        "status": "ROBUST_BUT_DEGENERATE" if (robust and not informative) else (
            "ROBUST_AND_INFORMATIVE" if admissible_as_selector else "NONROBUST"
        ),
        "qw2191_strict_status": "NOT_CLOSED",
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
