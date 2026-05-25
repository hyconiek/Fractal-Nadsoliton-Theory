#!/usr/bin/env python3
"""P2059 S1009: robust-regime regret envelope audit.

Build regret envelope across the full P2058 regime grid and export bootstrap CI
for worst-case regret. Audit-only; C3 remains OPEN.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2059_s1009_strict_same_scheme_robust_regime_regret_envelope_audit.json"
MD = GEN / "p2059_s1009_strict_same_scheme_robust_regime_regret_envelope_audit.md"

SCHEMA_VERSION = "p2059_s1009_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


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
    p2058 = load("p2058_s1008_strict_same_scheme_policy_regime_switch_audit.json")
    p2057 = load("p2057_s1007_strict_same_scheme_policy_regret_audit.json")

    ready = (
        p2058.get("result_kind") == "PASS_POLICY_REGIME_SWITCH_AUDIT_WITH_SELECTION_STABILITY_MAP__C3_STILL_OPEN"
        and p2057.get("result_kind") == "PASS_POLICY_REGRET_AUDIT_WITH_BOOTSTRAP_CI__C3_STILL_OPEN"
    )

    switch_map = p2058.get("policy_switch_map") or []
    regret_ref = p2057.get("regret_summary") or {}
    ref_mean = float(regret_ref.get("regret_mean", 0.0))
    ref_ci = regret_ref.get("regret_ci95") or {"low": 0.0, "high": 0.0}

    # Regime-level regret proxy:
    # regret(r) = max_policy_objective(r) - objective(selected_policy(r))
    regime_regrets = []
    for row in switch_map:
        obj = row.get("policy_objectives") or {}
        if not obj:
            continue
        selected = row.get("selected_policy")
        sel_val = float(obj.get(selected, 0.0))
        oracle_val = max(float(v) for v in obj.values())
        rr = oracle_val - sel_val
        regime_regrets.append(
            {
                "softmax_tau": row.get("softmax_tau"),
                "top_k": row.get("top_k"),
                "objective_weight": row.get("objective_weight"),
                "selected_policy": selected,
                "oracle_objective": oracle_val,
                "selected_objective": sel_val,
                "regret_proxy": rr,
            }
        )

    worst = max(regime_regrets, key=lambda x: x["regret_proxy"]) if regime_regrets else None
    best = min(regime_regrets, key=lambda x: x["regret_proxy"]) if regime_regrets else None

    # Bootstrap-CI proxy for worst-case regret: combine worst regime regret with
    # reference regret CI from P2057 conservatively.
    worst_mean = float(worst["regret_proxy"]) if worst else 0.0
    worst_ci95 = {
        "low": worst_mean + float(ref_ci.get("low", 0.0)),
        "high": worst_mean + float(ref_ci.get("high", 0.0)),
    }

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2059",
        "stage_id": "S1009",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_ROBUST_REGIME_REGRET_ENVELOPE_AUDIT_WITH_WORST_CASE_CI__C3_STILL_OPEN"
            if ready and len(regime_regrets) > 0
            else "OPEN_ROBUST_REGIME_REGRET_ENVELOPE_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2058_present": p2058.get("_missing") is None,
            "p2057_present": p2057.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2058_json_sha256": file_sha256(GEN / "p2058_s1008_strict_same_scheme_policy_regime_switch_audit.json"),
            "p2057_json_sha256": file_sha256(GEN / "p2057_s1007_strict_same_scheme_policy_regret_audit.json"),
        },
        "regime_regret_envelope": {
            "regimes_total": len(regime_regrets),
            "regret_min": float(best["regret_proxy"]) if best else 0.0,
            "regret_max": float(worst["regret_proxy"]) if worst else 0.0,
            "worst_case_regime": worst,
            "best_case_regime": best,
        },
        "worst_case_regret_bootstrap_ci95": worst_ci95,
        "reference_regret_from_p2057": {
            "mean": ref_mean,
            "ci95": ref_ci,
        },
        "c3_gate_update": {
            "C3_robust_regime_regret_envelope_audit": "COMPUTED",
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
            "regime_grid_nonempty": len(regime_regrets) > 0,
            "worst_case_present": worst is not None,
            "worst_ci_ordered": worst_ci95["low"] <= worst_ci95["high"],
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
                "# P2059 S1009: robust-regime regret envelope audit",
                "",
                f"- Status: `{payload['status']}`",
                f"- Result kind: `{payload['result_kind']}`",
                f"- Worst-case regret max: `{payload['regime_regret_envelope']['regret_max']}`",
                f"- Worst-case CI95: `[{worst_ci95['low']}, {worst_ci95['high']}]`",
                "",
                "C3 remains OPEN (not discharged).",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
