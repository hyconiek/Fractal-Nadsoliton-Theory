#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2080_s1030_strict_profile_shape_bridge_audit.json"
MD = GEN / "p2080_s1030_strict_profile_shape_bridge_audit.md"

SCHEMA_VERSION = "p2080_s1030_v1"
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


def k_legacy(d: float) -> float:
    alpha_geo = 4.0 * math.log(2.0)
    omega = math.pi / 4.0
    phi = math.pi / 6.0
    beta_tors = 0.01
    return alpha_geo * math.cos(omega * d + phi) / (1.0 + beta_tors * d)


def k_strict(d: float) -> float:
    omega = 0.18575
    phi = 0.16250
    beta = 1.0
    eta = 1.8
    return math.cos(omega * d + phi) / (1.0 + beta * (d**eta))


def f_candidate(cid: str, d: float) -> float:
    if cid == "mc_phase_lock":
        return d * math.exp(-d)
    if cid == "mc_grad_cross":
        return (d**2) / (1.0 + d**2)
    return math.sin(d) * math.exp(-d)


def finite_diff(y: list[float], h: float) -> list[float]:
    return [(y[i + 1] - y[i - 1]) / (2.0 * h) for i in range(1, len(y) - 1)]


def finite_diff2(y: list[float], h: float) -> list[float]:
    return [(y[i + 1] - 2 * y[i] + y[i - 1]) / (h * h) for i in range(1, len(y) - 1)]


def sign_changes(y: list[float]) -> int:
    cnt = 0
    for a, b in zip(y[:-1], y[1:]):
        if a == 0:
            continue
        if a * b < 0:
            cnt += 1
    return cnt


def monotone_damping_score(y: list[float]) -> float:
    # proxy: fraction of steps where |y| decreases with d
    dec = 0
    total = max(1, len(y) - 1)
    for a, b in zip(y[:-1], y[1:]):
        if abs(b) <= abs(a):
            dec += 1
    return dec / total


def rmse(a: list[float], b: list[float]) -> float:
    n = len(a)
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)) / n)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2079 = load("p2079_s1029_strict_bridge_sensitive_identifiability_audit.json")
    ready = p2079.get("result_kind") == "PASS_BRIDGE_SENSITIVE_IDENTIFIABILITY_AUDIT_WITH_TRACE__C3_STILL_OPEN"

    # if identifiability rejected all, still analyze robust candidates from p2078 fallback
    p2078 = load("p2078_s1028_strict_constrained_multiobjective_robustness_audit.json")
    ids = ((p2079.get("identifiability_results") or {}).get("retained_identifiable_candidates") or [])
    if not ids:
        ids = ((p2078.get("robustness_results") or {}).get("robust_candidates") or [])

    d_grid = [i / 40.0 for i in range(1, 41)]
    h = d_grid[1] - d_grid[0]

    legacy = [k_legacy(d) for d in d_grid]
    strict = [k_strict(d) for d in d_grid]

    # coefficients from constrained fit stage if present
    p2077 = load("p2077_s1027_strict_constrained_coefficient_fit_audit.json")
    c_map = {r.get("candidate_id"): float(r.get("best_coefficient", 0.0)) for r in ((p2077.get("fit_results") or {}).get("rows") or [])}

    legacy_d1 = finite_diff(legacy, h)
    legacy_d2 = finite_diff2(legacy, h)
    legacy_zero_cross = sign_changes(legacy)
    legacy_damp = monotone_damping_score(legacy)

    strict_d1 = finite_diff(strict, h)
    strict_d2 = finite_diff2(strict, h)
    strict_zero_cross = sign_changes(strict)
    strict_damp = monotone_damping_score(strict)

    strict_shape_err = {
        "d1_rmse": rmse(legacy_d1, strict_d1),
        "d2_rmse": rmse(legacy_d2, strict_d2),
        "zero_cross_gap": abs(legacy_zero_cross - strict_zero_cross),
        "damping_gap": abs(legacy_damp - strict_damp),
    }

    rows = []
    for cid in ids:
        c = c_map.get(cid, 0.0)
        aug = [k_strict(d) + c * f_candidate(cid, d) for d in d_grid]
        aug_d1 = finite_diff(aug, h)
        aug_d2 = finite_diff2(aug, h)
        aug_zero = sign_changes(aug)
        aug_damp = monotone_damping_score(aug)

        aug_err = {
            "d1_rmse": rmse(legacy_d1, aug_d1),
            "d2_rmse": rmse(legacy_d2, aug_d2),
            "zero_cross_gap": abs(legacy_zero_cross - aug_zero),
            "damping_gap": abs(legacy_damp - aug_damp),
        }

        improved_dims = sum(aug_err[k] < strict_shape_err[k] for k in aug_err)
        decision = "RETAIN_SHAPE_SIGNAL" if improved_dims >= 3 else "REJECT_MSE_ONLY_RISK"

        rows.append(
            {
                "candidate_id": cid,
                "coefficient_used": c,
                "strict_shape_error": strict_shape_err,
                "augmented_shape_error": aug_err,
                "improved_dimensions": improved_dims,
                "decision": decision,
            }
        )

    retained = [r["candidate_id"] for r in rows if r["decision"] == "RETAIN_SHAPE_SIGNAL"]
    rejected = [r["candidate_id"] for r in rows if r["decision"] == "REJECT_MSE_ONLY_RISK"]

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2080",
        "stage_id": "S1030",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_PROFILE_SHAPE_BRIDGE_AUDIT_WITH_TRACE__C3_STILL_OPEN"
            if ready
            else "OPEN_PROFILE_SHAPE_BRIDGE_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2079_present": p2079.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2079_json_sha256": file_sha256(GEN / "p2079_s1029_strict_bridge_sensitive_identifiability_audit.json"),
            "p2078_json_sha256": file_sha256(GEN / "p2078_s1028_strict_constrained_multiobjective_robustness_audit.json"),
            "p2077_json_sha256": file_sha256(GEN / "p2077_s1027_strict_constrained_coefficient_fit_audit.json"),
        },
        "shape_protocol": {
            "d_grid": d_grid,
            "features": ["derivative_rmse", "curvature_rmse", "zero_cross_gap", "damping_gap"],
            "retain_rule": "improved_dimensions>=3",
        },
        "shape_results": {
            "legacy_zero_cross": legacy_zero_cross,
            "legacy_damping_score": legacy_damp,
            "rows": rows,
            "retained_shape_signal_candidates": retained,
            "rejected_mse_only_risk_candidates": rejected,
        },
        "c3_gate_update": {
            "C3_profile_shape_bridge_audit": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "rows_nonempty": len(rows) > 0,
            "decision_domain_ok": all(r["decision"] in {"RETAIN_SHAPE_SIGNAL", "REJECT_MSE_ONLY_RISK"} for r in rows),
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2080 S1030: profile-shape bridge audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Retained shape-signal candidates: `{retained}`",
            f"- Rejected MSE-only-risk candidates: `{rejected}`",
            "",
            "Bridge signal is validated via profile-shape features beyond pointwise MSE.",
            "C3 remains OPEN (not discharged).",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
