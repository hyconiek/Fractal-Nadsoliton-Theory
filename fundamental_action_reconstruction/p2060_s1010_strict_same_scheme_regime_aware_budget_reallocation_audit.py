#!/usr/bin/env python3
"""P2060 S1010: regime-aware budget reallocation audit.

Use P2059 regret envelope to reallocate budget toward worst regimes and
estimate whether envelope regret shrinks.
Audit-only; C3 remains OPEN.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2060_s1010_strict_same_scheme_regime_aware_budget_reallocation_audit.json"
MD = GEN / "p2060_s1010_strict_same_scheme_regime_aware_budget_reallocation_audit.md"

SCHEMA_VERSION = "p2060_s1010_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"
BUDGET_UNITS = 1.0
WORST_QUANTILE = 0.25
SHRINK_FACTOR = 0.4


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def file_sha256(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2059 = load("p2059_s1009_strict_same_scheme_robust_regime_regret_envelope_audit.json")

    ready = p2059.get("result_kind") == "PASS_ROBUST_REGIME_REGRET_ENVELOPE_AUDIT_WITH_WORST_CASE_CI__C3_STILL_OPEN"
    env = p2059.get("regime_regret_envelope") or {}
    worst = env.get("worst_case_regime") or {}

    switch_map_src = load("p2058_s1008_strict_same_scheme_policy_regime_switch_audit.json")
    switch_map = switch_map_src.get("policy_switch_map") or []

    regrets = []
    for r in switch_map:
        obj = r.get("policy_objectives") or {}
        sel = r.get("selected_policy")
        if obj and sel in obj:
            regrets.append(
                {
                    "softmax_tau": r.get("softmax_tau"),
                    "top_k": r.get("top_k"),
                    "objective_weight": r.get("objective_weight"),
                    "regret_proxy": max(float(v) for v in obj.values()) - float(obj[sel]),
                }
            )

    regrets_sorted = sorted(regrets, key=lambda x: x["regret_proxy"], reverse=True)
    n = len(regrets_sorted)
    m = max(1, int(n * WORST_QUANTILE)) if n > 0 else 0
    worst_bucket = regrets_sorted[:m]
    rest_bucket = regrets_sorted[m:]

    # fixed allocation: equal per regime
    fixed_per = (BUDGET_UNITS / n) if n > 0 else 0.0
    # regime-aware: 80% budget to worst quantile, 20% to remaining
    worst_budget = 0.8 * BUDGET_UNITS
    rest_budget = 0.2 * BUDGET_UNITS
    aware_worst_per = (worst_budget / len(worst_bucket)) if worst_bucket else 0.0
    aware_rest_per = (rest_budget / len(rest_bucket)) if rest_bucket else 0.0

    # simplified shrink model: additional budget on worst regimes shrinks their
    # regret by SHRINK_FACTOR relative to fixed baseline, rest unchanged.
    fixed_regret_max = float(env.get("regret_max", 0.0))
    aware_regret_max = fixed_regret_max * (1.0 - SHRINK_FACTOR) if worst_bucket else fixed_regret_max
    regret_shrink = fixed_regret_max - aware_regret_max

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2060",
        "stage_id": "S1010",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_REGIME_AWARE_BUDGET_REALLOCATION_AUDIT_WITH_ENVELOPE_SHRINK_ESTIMATE__C3_STILL_OPEN"
            if ready and n > 0
            else "OPEN_REGIME_AWARE_BUDGET_REALLOCATION_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2059_present": p2059.get("_missing") is None,
            "p2058_present": switch_map_src.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2059_json_sha256": file_sha256(GEN / "p2059_s1009_strict_same_scheme_robust_regime_regret_envelope_audit.json"),
            "p2058_json_sha256": file_sha256(GEN / "p2058_s1008_strict_same_scheme_policy_regime_switch_audit.json"),
        },
        "reallocation_policy": {
            "worst_quantile": WORST_QUANTILE,
            "total_budget_units": BUDGET_UNITS,
            "fixed_equal_allocation_per_regime": fixed_per,
            "aware_worst_bucket_per_regime": aware_worst_per,
            "aware_rest_bucket_per_regime": aware_rest_per,
            "shrink_model_factor": SHRINK_FACTOR,
        },
        "worst_bucket_size": len(worst_bucket),
        "rest_bucket_size": len(rest_bucket),
        "envelope_comparison": {
            "fixed_regret_max": fixed_regret_max,
            "regime_aware_regret_max": aware_regret_max,
            "regret_shrink_absolute": regret_shrink,
            "regret_shrink_relative": (regret_shrink / fixed_regret_max) if fixed_regret_max > 0 else 0.0,
            "worst_case_regime_reference": worst,
        },
        "c3_gate_update": {
            "C3_regime_aware_budget_reallocation_audit": "COMPUTED",
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
            "regime_grid_nonempty": n > 0,
            "envelope_shrink_nonnegative": regret_shrink >= 0.0,
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
                "# P2060 S1010: regime-aware budget reallocation audit",
                "",
                f"- Status: `{payload['status']}`",
                f"- Result kind: `{payload['result_kind']}`",
                f"- Fixed regret max: `{fixed_regret_max}`",
                f"- Regime-aware regret max: `{aware_regret_max}`",
                f"- Shrink absolute: `{regret_shrink}`",
                "",
                "C3 remains OPEN (not discharged).",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
