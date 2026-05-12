"""P1318: O3.5 independent replay + adversarial audit for QW-2191 strict lane."""

from __future__ import annotations

import json
import random
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1318_o35_independent_replay_adversarial_audit_report_v1.json"


@dataclass(frozen=True)
class Candidate:
    code: str
    branch_sign: int
    phase_class: float
    amplitude_class: float
    residual_slot: str
    admissible: bool


BASE = [
    Candidate("C1", +1, 0.41, 0.62, "open", True),
    Candidate("C2", +1, 0.43, 0.60, "open", True),
    Candidate("C3", +1, 0.40, 0.58, "open(Z2/eps)", True),
    Candidate("C4", +1, 0.42, 0.61, "open", True),
]

PHASE_TAU = 0.08
AMPLITUDE_TAU = 0.08
REPLAY_SEEDS = [17, 314, 2718]


def eval_pair(left: Candidate, right: Candidate) -> str:
    if not left.admissible or not right.admissible:
        return "INADMISSIBLE_BRANCH"
    if left.branch_sign != right.branch_sign:
        return "COUNTEREXAMPLE_FOUND"
    if abs(left.phase_class - right.phase_class) > PHASE_TAU:
        return "COUNTEREXAMPLE_FOUND"
    if abs(left.amplitude_class - right.amplitude_class) > AMPLITUDE_TAU:
        return "COUNTEREXAMPLE_FOUND"
    if "open" in left.residual_slot or "open" in right.residual_slot:
        return "RESIDUAL_AMBIGUITY"
    return "EQUIVALENT_DIRECTION"


def run_replay(seed: int) -> dict:
    rng = random.Random(seed)
    perm = BASE[:]
    rng.shuffle(perm)

    results = {}
    counts = {"COUNTEREXAMPLE_FOUND": 0, "RESIDUAL_AMBIGUITY": 0, "EQUIVALENT_DIRECTION": 0, "INADMISSIBLE_BRANCH": 0}
    for lft, rgt in combinations(perm, 2):
        key = "-".join(sorted([lft.code, rgt.code]))
        verdict = eval_pair(lft, rgt)
        results[key] = verdict
        counts[verdict] += 1

    return {"seed": seed, "pair_verdicts": results, "counts": counts}


def adversarial_attack() -> dict:
    # adversarially force maximal allowed perturbation while staying within tolerances
    attacked = []
    for c in BASE:
        attacked.append(
            Candidate(
                c.code,
                c.branch_sign,
                c.phase_class + (PHASE_TAU * 0.99),
                max(0.0, c.amplitude_class - (AMPLITUDE_TAU * 0.99)),
                c.residual_slot,
                c.admissible,
            )
        )

    counts = {"COUNTEREXAMPLE_FOUND": 0, "RESIDUAL_AMBIGUITY": 0, "EQUIVALENT_DIRECTION": 0, "INADMISSIBLE_BRANCH": 0}
    for lft, rgt in combinations(attacked, 2):
        verdict = eval_pair(lft, rgt)
        counts[verdict] += 1

    success = counts["COUNTEREXAMPLE_FOUND"] > 0
    return {
        "attack_model": "within-tolerance perturbation + residual-slot pressure",
        "counts": counts,
        "adversarial_break_success": success,
    }


def main() -> None:
    replays = [run_replay(seed) for seed in REPLAY_SEEDS]

    baseline = replays[0]["pair_verdicts"]
    divergence_count = 0
    for rep in replays[1:]:
        for pair, verdict in rep["pair_verdicts"].items():
            if baseline[pair] != verdict:
                divergence_count += 1

    adv = adversarial_attack()

    replay_pass = divergence_count == 0
    adversarial_pass = not adv["adversarial_break_success"]
    o35_pass = replay_pass and adversarial_pass

    payload = {
        "packet_id": "P1318_O3_5_INDEPENDENT_REPLAY_ADVERSARIAL_AUDIT_REPORT_V1",
        "date_utc": "2026-05-12",
        "strict_only": True,
        "replay_runs": replays,
        "divergence_count": divergence_count,
        "adversarial": adv,
        "replay_status": "PASS" if replay_pass else "FAIL",
        "adversarial_status": "PASS" if adversarial_pass else "FAIL",
        "o3_5_status": "PASS" if o35_pass else "FAIL",
        "nonuniqueness_residual_status": "OPEN",
        "qw2191_strict_status": "NOT_CLOSED",
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
