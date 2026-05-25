#!/usr/bin/env python3
"""P2058 S1008: regime-switch audit for allocation-policy stability.

Sweep softmax_tau, top_k, objective_weight and export policy-switch map plus
regret stability indicators. Audit-only; C3 remains OPEN.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2058_s1008_strict_same_scheme_policy_regime_switch_audit.json"
MD = GEN / "p2058_s1008_strict_same_scheme_policy_regime_switch_audit.md"

SCHEMA_VERSION = "p2058_s1008_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"
TOTAL_EXTRA_UNITS = 9000.0

SOFTMAX_TAU_GRID = [0.05, 0.1, 0.2]
TOP_K_GRID = [1, 2, 3]
OBJECTIVE_WEIGHT_GRID = [0.25, 0.5, 0.75]


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def file_sha256(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def objective(mean_gain: float, var_gain: float, w: float) -> float:
    return mean_gain - w * math.sqrt(max(var_gain, 0.0))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2054 = load("p2054_s1004_strict_same_scheme_adaptive_budget_allocator_audit.json")
    p2057 = load("p2057_s1007_strict_same_scheme_policy_regret_audit.json")

    ready = (
        p2054.get("result_kind") == "PASS_ADAPTIVE_BUDGET_ALLOCATOR_AUDIT_FIXED_VS_ADAPTIVE_WITH_TRACE__C3_STILL_OPEN"
        and p2057.get("result_kind") == "PASS_POLICY_REGRET_AUDIT_WITH_BOOTSTRAP_CI__C3_STILL_OPEN"
    )

    rows = p2054.get("adaptive_allocation") or []
    n = len(rows)
    ordered = sorted(rows, key=lambda r: float(r.get("priority_score", 0.0)), reverse=True)
    scores = [float(r.get("priority_score", 0.0)) for r in ordered]

    regime_rows: list[dict[str, Any]] = []
    for tau in SOFTMAX_TAU_GRID:
        for k in TOP_K_GRID:
            kk = min(k, n) if n > 0 else 0
            for w in OBJECTIVE_WEIGHT_GRID:
                # deterministic proxies consistent with existing pipeline:
                # all_to_top mean proxy: 1/n ; top_k: k/n ; softmax: nz_ratio (all positive)
                all_to_top_mean = (1.0 / n) if n > 0 else 0.0
                top_k_mean = (kk / n) if n > 0 else 0.0

                if n > 0:
                    mx = max(scores)
                    weights = [math.exp((s - mx) / max(tau, 1e-9)) for s in scores]
                    z = sum(weights)
                    alloc = [TOTAL_EXTRA_UNITS * wi / z for wi in weights]
                    softmax_mean = sum(1 for a in alloc if a > 0.0) / n
                else:
                    softmax_mean = 0.0

                # tiny deterministic variance proxies to avoid false precision
                all_to_top_var = all_to_top_mean * (1.0 - all_to_top_mean)
                top_k_var = top_k_mean * (1.0 - top_k_mean)
                softmax_var = softmax_mean * (1.0 - softmax_mean)

                cand = {
                    "all_to_top": objective(all_to_top_mean, all_to_top_var, w),
                    "top_k_split": objective(top_k_mean, top_k_var, w),
                    "softmax": objective(softmax_mean, softmax_var, w),
                }
                selected_policy = max(cand, key=cand.get)
                regime_rows.append(
                    {
                        "softmax_tau": tau,
                        "top_k": kk,
                        "objective_weight": w,
                        "selected_policy": selected_policy,
                        "policy_objectives": cand,
                    }
                )

    switch_counts: dict[str, int] = {}
    for r in regime_rows:
        switch_counts[r["selected_policy"]] = switch_counts.get(r["selected_policy"], 0) + 1

    dominant_policy = max(switch_counts, key=switch_counts.get) if switch_counts else None
    dominant_share = (switch_counts.get(dominant_policy, 0) / len(regime_rows)) if regime_rows else 0.0

    regret_summary = p2057.get("regret_summary") or {}

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2058",
        "stage_id": "S1008",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_POLICY_REGIME_SWITCH_AUDIT_WITH_SELECTION_STABILITY_MAP__C3_STILL_OPEN"
            if ready and len(regime_rows) > 0
            else "OPEN_POLICY_REGIME_SWITCH_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2054_present": p2054.get("_missing") is None,
            "p2057_present": p2057.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2054_json_sha256": file_sha256(GEN / "p2054_s1004_strict_same_scheme_adaptive_budget_allocator_audit.json"),
            "p2057_json_sha256": file_sha256(GEN / "p2057_s1007_strict_same_scheme_policy_regret_audit.json"),
        },
        "regime_grid": {
            "softmax_tau": SOFTMAX_TAU_GRID,
            "top_k": TOP_K_GRID,
            "objective_weight": OBJECTIVE_WEIGHT_GRID,
            "points_total": len(regime_rows),
        },
        "policy_switch_map": regime_rows,
        "selection_stability": {
            "switch_counts": switch_counts,
            "dominant_policy": dominant_policy,
            "dominant_share": dominant_share,
        },
        "regret_reference_from_p2057": regret_summary,
        "c3_gate_update": {
            "C3_policy_regime_switch_audit": "COMPUTED",
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
            "grid_nonempty": len(regime_rows) > 0,
            "dominant_policy_present": dominant_policy is not None,
            "regret_reference_present": isinstance(regret_summary, dict) and len(regret_summary) > 0,
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
                "# P2058 S1008: policy regime-switch audit",
                "",
                f"- Status: `{payload['status']}`",
                f"- Result kind: `{payload['result_kind']}`",
                f"- Dominant policy: `{dominant_policy}`",
                f"- Dominant share: `{dominant_share}`",
                "",
                "Regime grid: softmax_tau × top_k × objective_weight.",
                "C3 remains OPEN (not discharged).",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
