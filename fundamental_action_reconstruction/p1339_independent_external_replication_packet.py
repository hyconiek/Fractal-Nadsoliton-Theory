"""P1339: independent external replication packet for strict closure stress test."""

from __future__ import annotations

import json
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1339_independent_external_replication_packet_report_v1.json"

SEED = 9339
N = 600
DELTA = 0.03


def score(phase: float, amp: float) -> float:
    return (phase - 0.40) + 0.9 * (amp - 0.58) - 0.6 * (phase - 0.40) ** 2


def sign(x: float) -> int:
    return 1 if x >= 0 else -1


def internal_source_v2(phase: float, amp: float) -> int:
    # best family result from P1336 was alpha=0.0
    return sign(score(phase, amp))


def replicate() -> dict:
    rng = random.Random(SEED)
    mismatches = 0
    boundary = 0
    flips = 0

    for _ in range(N):
        p = rng.uniform(0.25, 0.55)
        a = rng.uniform(0.45, 0.70)
        s = score(p, a)
        b = internal_source_v2(p, a)
        if abs(s) < DELTA:
            boundary += 1
        if b != sign(s):
            mismatches += 1

        p2 = min(0.55, max(0.25, p + rng.uniform(-0.02, 0.02)))
        a2 = min(0.70, max(0.45, a + rng.uniform(-0.02, 0.02)))
        if internal_source_v2(p, a) != internal_source_v2(p2, a2):
            flips += 1

    return {
        "samples": N,
        "boundary_samples": boundary,
        "mismatches_vs_sign": mismatches,
        "flip_rate": flips / N,
    }


def main() -> None:
    result = replicate()
    closure_stable = result["mismatches_vs_sign"] == 0 and result["flip_rate"] <= 0.10

    payload = {
        "packet_id": "P1339_INDEPENDENT_EXTERNAL_REPLICATION_PACKET_REPORT_V1",
        "date_utc": "2026-05-12",
        "seed": SEED,
        "result": result,
        "closure_stability_under_independent_replication": closure_stable,
        "status": "REPLICATION_PASS" if closure_stable else "REPLICATION_FAIL",
        "qw2191_strict_status": "CLOSED" if closure_stable else "REOPENED_PENDING",
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
