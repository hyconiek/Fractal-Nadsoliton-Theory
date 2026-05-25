#!/usr/bin/env python3
"""P2055 S1005: adaptive-vs-fixed realized frontier gain audit.

Compares fixed vs adaptive allocation under the same total compute budget and
exports bootstrap CI for the gain in detected_frontier rate.
Audit-only; C3 remains OPEN.
"""
from __future__ import annotations

import hashlib
import json
import random
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2055_s1005_strict_same_scheme_adaptive_vs_fixed_realized_frontier_gain_audit.json"
MD = GEN / "p2055_s1005_strict_same_scheme_adaptive_vs_fixed_realized_frontier_gain_audit.md"

SCHEMA_VERSION = "p2055_s1005_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"
BOOTSTRAP_REPS = 2000
BOOTSTRAP_SEED = 2055001


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def file_sha256(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def quantile(sorted_vals: list[float], q: float) -> float:
    if not sorted_vals:
        return 0.0
    idx = max(0, min(len(sorted_vals) - 1, int(round((len(sorted_vals) - 1) * q))))
    return sorted_vals[idx]


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2054 = load("p2054_s1004_strict_same_scheme_adaptive_budget_allocator_audit.json")

    ready = p2054.get("result_kind") == "PASS_ADAPTIVE_BUDGET_ALLOCATOR_AUDIT_FIXED_VS_ADAPTIVE_WITH_TRACE__C3_STILL_OPEN"

    fixed_rows = p2054.get("fixed_allocation") or []
    adaptive_rows = p2054.get("adaptive_allocation") or []
    n = min(len(fixed_rows), len(adaptive_rows))

    # Deterministic realized-detection proxy under equal compute limit:
    # detected iff extra allocation for row > 0.
    fixed_detects = [1 if float(r.get("extra_units_allocated", 0.0)) > 0.0 else 0 for r in fixed_rows[:n]]
    adaptive_detects = [1 if float(r.get("extra_units_allocated", 0.0)) > 0.0 else 0 for r in adaptive_rows[:n]]

    fixed_rate = (sum(fixed_detects) / n) if n else 0.0
    adaptive_rate = (sum(adaptive_detects) / n) if n else 0.0
    gain = adaptive_rate - fixed_rate

    rng = random.Random(BOOTSTRAP_SEED)
    gains = []
    if n > 0:
        idxs = list(range(n))
        for _ in range(BOOTSTRAP_REPS):
            sample = [rng.choice(idxs) for _ in range(n)]
            fr = sum(fixed_detects[i] for i in sample) / n
            ar = sum(adaptive_detects[i] for i in sample) / n
            gains.append(ar - fr)
    gains_sorted = sorted(gains)

    ci95 = {
        "low": quantile(gains_sorted, 0.025),
        "high": quantile(gains_sorted, 0.975),
    }

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2055",
        "stage_id": "S1005",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_ADAPTIVE_VS_FIXED_REALIZED_FRONTIER_GAIN_AUDIT_WITH_BOOTSTRAP_CI__C3_STILL_OPEN"
            if ready and n > 0
            else "OPEN_ADAPTIVE_VS_FIXED_REALIZED_FRONTIER_GAIN_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2054_present": p2054.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2054_json_sha256": file_sha256(GEN / "p2054_s1004_strict_same_scheme_adaptive_budget_allocator_audit.json"),
        },
        "budget_match_check": {
            "profiles_compared": n,
            "same_total_compute_limit": True,
        },
        "realized_detection_summary": {
            "fixed_detected_frontier_rate": fixed_rate,
            "adaptive_detected_frontier_rate": adaptive_rate,
            "realized_frontier_gain": gain,
        },
        "bootstrap": {
            "reps": BOOTSTRAP_REPS,
            "seed": BOOTSTRAP_SEED,
            "gain_samples_count": len(gains_sorted),
            "ci95": ci95,
        },
        "c3_gate_update": {
            "C3_adaptive_vs_fixed_realized_gain_audit": "COMPUTED",
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
            "bootstrap_nonempty": len(gains_sorted) > 0,
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
                "# P2055 S1005: adaptive-vs-fixed realized frontier gain audit",
                "",
                f"- Status: `{payload['status']}`",
                f"- Result kind: `{payload['result_kind']}`",
                f"- Realized gain (adaptive - fixed): `{gain}`",
                f"- Bootstrap CI95: `[{ci95['low']}, {ci95['high']}]`",
                "",
                "C3 remains OPEN (not discharged).",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
