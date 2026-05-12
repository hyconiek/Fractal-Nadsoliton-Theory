"""P1326: independent replay/adversarial and O3 reintegration for S_sel_strict_v4."""

from __future__ import annotations

import json
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1326_v4_replay_adversarial_reintegration_report_v1.json"

SEEDS = [1326, 2326, 3326]
N_SAMPLES = 400
PERTURB = 0.02


def s_sel_strict_v4(phase: float, amp: float) -> int:
    score = (phase - 0.40) + 0.9 * (amp - 0.58) - 0.6 * (phase - 0.40) ** 2
    return 1 if score >= 0.0 else -1


def sample(rng: random.Random) -> tuple[float, float]:
    return rng.uniform(0.25, 0.55), rng.uniform(0.45, 0.70)


def run_once(seed: int) -> dict:
    rng = random.Random(seed)
    pos = 0
    neg = 0
    flips = 0
    for _ in range(N_SAMPLES):
        p, a = sample(rng)
        s = s_sel_strict_v4(p, a)
        if s > 0:
            pos += 1
        else:
            neg += 1
        p2 = min(0.55, max(0.25, p + rng.uniform(-PERTURB, PERTURB)))
        a2 = min(0.70, max(0.45, a + rng.uniform(-PERTURB, PERTURB)))
        s2 = s_sel_strict_v4(p2, a2)
        if s != s2:
            flips += 1
    return {
        "seed": seed,
        "positive_count": pos,
        "negative_count": neg,
        "sign_diversity": 2 if pos and neg else 1,
        "flip_count": flips,
        "flip_rate": flips / N_SAMPLES,
    }


def adversarial_edge_attack() -> dict:
    # sweep near estimated boundary in phase/amp space
    hits = 0
    total = 0
    unstable = 0
    for phase in [0.34 + 0.005 * i for i in range(25)]:
        for amp in [0.53 + 0.005 * j for j in range(25)]:
            total += 1
            s = s_sel_strict_v4(phase, amp)
            s2 = s_sel_strict_v4(min(0.55, phase + PERTURB), max(0.45, amp - PERTURB))
            if s != s2:
                unstable += 1
            if s > 0:
                hits += 1
    return {
        "grid_points": total,
        "positive_hits": hits,
        "unstable_points": unstable,
        "unstable_rate": unstable / total,
    }


def main() -> None:
    runs = [run_once(s) for s in SEEDS]
    replay_pass = all(r["sign_diversity"] == 2 and r["flip_rate"] <= 0.08 for r in runs)

    adv = adversarial_edge_attack()
    adversarial_pass = adv["unstable_rate"] <= 0.20

    o3_reintegration_ready = replay_pass and adversarial_pass

    payload = {
        "packet_id": "P1326_V4_REPLAY_ADVERSARIAL_REINTEGRATION_REPORT_V1",
        "date_utc": "2026-05-12",
        "model": "S_sel_strict_v4",
        "runs": runs,
        "replay_pass": replay_pass,
        "adversarial": adv,
        "adversarial_pass": adversarial_pass,
        "o3_reintegration_ready": o3_reintegration_ready,
        "strict_core_selector_source_exported": False,
        "qw2191_strict_status": "NOT_CLOSED",
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
