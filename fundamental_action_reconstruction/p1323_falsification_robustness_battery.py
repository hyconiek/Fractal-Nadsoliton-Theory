"""P1323: falsification and robustness battery for S_sel_strict_v2."""

from __future__ import annotations

import json
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1323_falsification_robustness_battery_report_v1.json"

N_SAMPLES = 200
SEED = 1323


def s_sel_strict_v2(phase: float, amp: float) -> int:
    score = (phase - 0.30) + 0.50 * (amp - 0.50)
    return 1 if score >= 0.0 else -1


def generate_sample(rng: random.Random) -> tuple[float, float]:
    phase = rng.uniform(0.25, 0.55)
    amp = rng.uniform(0.45, 0.70)
    return phase, amp


def main() -> None:
    rng = random.Random(SEED)
    negatives = []
    positives = 0

    for i in range(N_SAMPLES):
        phase, amp = generate_sample(rng)
        sign = s_sel_strict_v2(phase, amp)
        if sign < 0:
            negatives.append({"id": i, "phase": phase, "amp": amp})
        else:
            positives += 1

    robust_positive = len(negatives) == 0

    payload = {
        "packet_id": "P1323_FALSIFICATION_ROBUSTNESS_BATTERY_REPORT_V1",
        "date_utc": "2026-05-12",
        "model": "S_sel_strict_v2",
        "sample_count": N_SAMPLES,
        "seed": SEED,
        "positive_count": positives,
        "negative_count": len(negatives),
        "counterexamples": negatives[:20],
        "robust_positive_on_test_window": robust_positive,
        "status": "ROBUST_ON_TEST_WINDOW" if robust_positive else "COUNTEREXAMPLES_FOUND",
        "strict_core_selector_source_exported": False,
        "qw2191_strict_status": "NOT_CLOSED",
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
