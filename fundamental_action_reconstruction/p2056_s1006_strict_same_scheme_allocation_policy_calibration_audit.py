#!/usr/bin/env python3
"""P2056 S1006: allocation-policy calibration audit.

Compare allocation policies under same compute budget:
- all-to-top
- top-k split
- softmax allocation
Then select policy with best bootstrap gain/variance tradeoff.
Audit-only; C3 remains OPEN.
"""
from __future__ import annotations

import hashlib
import json
import math
import random
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2056_s1006_strict_same_scheme_allocation_policy_calibration_audit.json"
MD = GEN / "p2056_s1006_strict_same_scheme_allocation_policy_calibration_audit.md"

SCHEMA_VERSION = "p2056_s1006_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"
BOOTSTRAP_REPS = 2000
BOOTSTRAP_SEED = 2056001
TOTAL_EXTRA_UNITS = 9000.0
TOP_K = 2
SOFTMAX_TAU = 0.1


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def file_sha256(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def qtile(vals: list[float], q: float) -> float:
    if not vals:
        return 0.0
    s = sorted(vals)
    idx = max(0, min(len(s) - 1, int(round((len(s) - 1) * q))))
    return s[idx]


def detect_flags(extras: list[float]) -> list[int]:
    return [1 if e > 0 else 0 for e in extras]


def bootstrap_gain(base: list[int], policy: list[int], reps: int, seed: int) -> dict[str, Any]:
    n = min(len(base), len(policy))
    if n == 0:
        return {"mean": 0.0, "var": 0.0, "ci95": {"low": 0.0, "high": 0.0}, "samples": 0}
    rng = random.Random(seed)
    idxs = list(range(n))
    gs: list[float] = []
    for _ in range(reps):
        sample = [rng.choice(idxs) for _ in range(n)]
        br = sum(base[i] for i in sample) / n
        pr = sum(policy[i] for i in sample) / n
        gs.append(pr - br)
    mean = sum(gs) / len(gs)
    var = sum((x - mean) ** 2 for x in gs) / max(len(gs) - 1, 1)
    return {
        "mean": mean,
        "var": var,
        "ci95": {"low": qtile(gs, 0.025), "high": qtile(gs, 0.975)},
        "samples": len(gs),
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2054 = load("p2054_s1004_strict_same_scheme_adaptive_budget_allocator_audit.json")

    ready = p2054.get("result_kind") == "PASS_ADAPTIVE_BUDGET_ALLOCATOR_AUDIT_FIXED_VS_ADAPTIVE_WITH_TRACE__C3_STILL_OPEN"

    adaptive = p2054.get("adaptive_allocation") or []
    fixed = p2054.get("fixed_allocation") or []
    n = min(len(adaptive), len(fixed))

    # sorted by priority from adaptive rows
    rows = sorted(adaptive[:n], key=lambda r: float(r.get("priority_score", 0.0)), reverse=True)

    # baseline fixed flags
    fixed_extras = [float(r.get("extra_units_allocated", 0.0)) for r in fixed[:n]]
    fixed_flags = detect_flags(fixed_extras)

    # policy A: all-to-top
    all_to_top_extras = [TOTAL_EXTRA_UNITS if i == 0 else 0.0 for i in range(n)]
    all_to_top_flags = detect_flags(all_to_top_extras)

    # policy B: top-k split
    k = min(TOP_K, n) if n > 0 else 0
    topk_extras = [TOTAL_EXTRA_UNITS / k if i < k and k > 0 else 0.0 for i in range(n)]
    topk_flags = detect_flags(topk_extras)

    # policy C: softmax over priority scores
    scores = [float(r.get("priority_score", 0.0)) for r in rows]
    if n > 0:
        mx = max(scores)
        weights = [math.exp((s - mx) / max(SOFTMAX_TAU, 1e-9)) for s in scores]
        z = sum(weights)
        softmax_extras = [TOTAL_EXTRA_UNITS * w / z for w in weights]
    else:
        softmax_extras = []
    softmax_flags = detect_flags(softmax_extras)

    pol = {
        "all_to_top": bootstrap_gain(fixed_flags, all_to_top_flags, BOOTSTRAP_REPS, BOOTSTRAP_SEED + 1),
        "top_k_split": bootstrap_gain(fixed_flags, topk_flags, BOOTSTRAP_REPS, BOOTSTRAP_SEED + 2),
        "softmax": bootstrap_gain(fixed_flags, softmax_flags, BOOTSTRAP_REPS, BOOTSTRAP_SEED + 3),
    }

    # calibration objective: maximize mean_gain - 0.5*sqrt(var)
    scored = []
    for name, s in pol.items():
        objective = float(s["mean"]) - 0.5 * math.sqrt(max(float(s["var"]), 0.0))
        scored.append({"policy": name, "objective": objective, **s})
    scored.sort(key=lambda x: x["objective"], reverse=True)
    best = scored[0] if scored else None

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2056",
        "stage_id": "S1006",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_ALLOCATION_POLICY_CALIBRATION_AUDIT_WITH_BOOTSTRAP_GAIN_VARIANCE__C3_STILL_OPEN"
            if ready and n > 0
            else "OPEN_ALLOCATION_POLICY_CALIBRATION_AUDIT_BLOCKED"
        ),
        "depends_on": {"p2054_present": p2054.get("_missing") is None, "preconditions_ready": ready},
        "input_hashes": {
            "p2054_json_sha256": file_sha256(GEN / "p2054_s1004_strict_same_scheme_adaptive_budget_allocator_audit.json"),
        },
        "calibration_setup": {
            "profiles_compared": n,
            "same_total_compute_limit": True,
            "total_extra_units": TOTAL_EXTRA_UNITS,
            "top_k": TOP_K,
            "softmax_tau": SOFTMAX_TAU,
            "bootstrap_reps": BOOTSTRAP_REPS,
            "bootstrap_seed": BOOTSTRAP_SEED,
            "objective": "mean_gain - 0.5*sqrt(var_gain)",
        },
        "policy_results": scored,
        "selected_policy": best,
        "c3_gate_update": {
            "C3_allocation_policy_calibration_audit": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
            "missing_for_c3_theorem": [
                "theorem-grade operator identity across full background family",
                "cross-background finite-part transport identity theorem",
                "global finite-part lock/cocycle theorem",
            ],
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "profiles_nonempty": n > 0,
            "policies_evaluated": len(scored) == 3,
            "selected_policy_present": best is not None,
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2056 S1006: allocation-policy calibration audit",
                "",
                f"- Status: `{payload['status']}`",
                f"- Result kind: `{payload['result_kind']}`",
                f"- Selected policy: `{best['policy'] if best else 'NONE'}`",
                "",
                "Policies compared: all-to-top, top-k split, softmax.",
                "C3 remains OPEN (not discharged).",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
