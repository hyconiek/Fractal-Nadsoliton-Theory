"""P1325: balanced v4 selector probe with revised non-contradictory gate."""

from __future__ import annotations

import json
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1325_balanced_v4_gate_revision_and_probe_report_v1.json"

N_SAMPLES = 400
SEED = 1325
PERTURB = 0.02


def s_sel_strict_v4(phase: float, amp: float) -> int:
    # balanced candidate: nonlinear blend to avoid degeneracy
    score = (phase - 0.40) + 0.9 * (amp - 0.58) - 0.6 * (phase - 0.40) ** 2
    return 1 if score >= 0.0 else -1


def generate_sample(rng: random.Random) -> tuple[float, float]:
    return rng.uniform(0.25, 0.55), rng.uniform(0.45, 0.70)


def main() -> None:
    rng = random.Random(SEED)
    rows = []
    positives = 0
    negatives = 0
    flip_count = 0

    for i in range(N_SAMPLES):
        phase, amp = generate_sample(rng)
        sign = s_sel_strict_v4(phase, amp)
        if sign > 0:
            positives += 1
        else:
            negatives += 1

        p2 = min(0.55, max(0.25, phase + rng.uniform(-PERTURB, PERTURB)))
        a2 = min(0.70, max(0.45, amp + rng.uniform(-PERTURB, PERTURB)))
        sign2 = s_sel_strict_v4(p2, a2)
        if sign != sign2:
            flip_count += 1

        if i < 20:
            rows.append({"id": i, "phase": phase, "amp": amp, "sign": sign, "perturbed_sign": sign2})

    sign_diversity = 2 if positives and negatives else 1
    flip_rate = flip_count / N_SAMPLES

    # Revised gate (non-contradictory):
    # informative if both signs appear; robust if local flip rate below threshold.
    informative = sign_diversity == 2
    robust_local = flip_rate <= 0.08
    admissible = informative and robust_local

    payload = {
        "packet_id": "P1325_BALANCED_V4_GATE_REVISION_AND_PROBE_REPORT_V1",
        "date_utc": "2026-05-12",
        "model": "S_sel_strict_v4",
        "formula": "sign((phase-0.40)+0.9*(amp-0.58)-0.6*(phase-0.40)^2)",
        "sample_count": N_SAMPLES,
        "seed": SEED,
        "positive_count": positives,
        "negative_count": negatives,
        "sign_diversity": sign_diversity,
        "flip_count": flip_count,
        "flip_rate": flip_rate,
        "revised_gate": {
            "informative_requires_sign_diversity_2": True,
            "robust_local_requires_flip_rate_le_0_08": True,
        },
        "informative": informative,
        "robust_local": robust_local,
        "admissible_as_strict_selector_candidate": admissible,
        "status": "CANDIDATE_V4_PASSES_REVISED_GATE" if admissible else "CANDIDATE_V4_FAILS_REVISED_GATE",
        "preview_rows": rows,
        "qw2191_strict_status": "NOT_CLOSED",
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
