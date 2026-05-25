#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2081_s1031_strict_feature_weighted_candidate_redesign_audit.json"
MD = GEN / "p2081_s1031_strict_feature_weighted_candidate_redesign_audit.md"

SCHEMA_VERSION = "p2081_s1031_v1"
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


def finite_diff(y: list[float], h: float) -> list[float]:
    return [(y[i + 1] - y[i - 1]) / (2.0 * h) for i in range(1, len(y) - 1)]


def finite_diff2(y: list[float], h: float) -> list[float]:
    return [(y[i + 1] - 2 * y[i] + y[i - 1]) / (h * h) for i in range(1, len(y) - 1)]


def rmse(a: list[float], b: list[float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)) / len(a))


def sign_changes(y: list[float]) -> int:
    c = 0
    for a, b in zip(y[:-1], y[1:]):
        if a != 0 and a * b < 0:
            c += 1
    return c


def damp_score(y: list[float]) -> float:
    total = max(1, len(y) - 1)
    return sum(1 for a, b in zip(y[:-1], y[1:]) if abs(b) <= abs(a)) / total


def f_redesign(cid: str, d: float) -> float:
    # redesigned to target derivative/curvature mismatch explicitly
    if cid == "mc_phase_lock":
        return d * (1.0 - d) * math.exp(-d)
    if cid == "mc_grad_cross":
        return (d * (1.0 - 0.5 * d)) / (1.0 + d * d)
    return d * math.exp(-d)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2080 = load("p2080_s1030_strict_profile_shape_bridge_audit.json")
    ready = p2080.get("result_kind") == "PASS_PROFILE_SHAPE_BRIDGE_AUDIT_WITH_TRACE__C3_STILL_OPEN"

    # redesign candidates from shape-rejected pool
    ids = ((p2080.get("shape_results") or {}).get("rejected_mse_only_risk_candidates") or [])
    if not ids:
        ids = [r.get("candidate_id") for r in ((p2080.get("shape_results") or {}).get("rows") or [])]

    d_grid = [i / 40.0 for i in range(1, 41)]
    h = d_grid[1] - d_grid[0]
    legacy = [k_legacy(d) for d in d_grid]
    strict = [k_strict(d) for d in d_grid]

    legacy_d1 = finite_diff(legacy, h)
    legacy_d2 = finite_diff2(legacy, h)
    strict_d1 = finite_diff(strict, h)
    strict_d2 = finite_diff2(strict, h)

    strict_shape = {
        "d1_rmse": rmse(legacy_d1, strict_d1),
        "d2_rmse": rmse(legacy_d2, strict_d2),
        "zero_gap": abs(sign_changes(legacy) - sign_changes(strict)),
        "damp_gap": abs(damp_score(legacy) - damp_score(strict)),
    }

    c_grid = [x / 100.0 for x in range(-20, 21)]
    lam = 0.05
    rows = []

    for cid in ids:
        best = None
        for c in c_grid:
            aug = [k_strict(d) + c * f_redesign(cid, d) for d in d_grid]
            aug_d1 = finite_diff(aug, h)
            aug_d2 = finite_diff2(aug, h)
            aug_shape = {
                "d1_rmse": rmse(legacy_d1, aug_d1),
                "d2_rmse": rmse(legacy_d2, aug_d2),
                "zero_gap": abs(sign_changes(legacy) - sign_changes(aug)),
                "damp_gap": abs(damp_score(legacy) - damp_score(aug)),
            }
            # weighted objective prioritizing derivative/curvature
            obj = 2.0 * aug_shape["d1_rmse"] + 2.0 * aug_shape["d2_rmse"] + aug_shape["zero_gap"] + aug_shape["damp_gap"] + lam * (c**2)
            if best is None or obj < best["obj"]:
                best = {"c": c, "shape": aug_shape, "obj": obj}

        assert best is not None
        improved_dims = sum(best["shape"][k] < strict_shape[k] for k in strict_shape)
        decision = "RETAIN_REDESIGNED" if improved_dims >= 3 else "REJECT_REDESIGNED"
        rows.append(
            {
                "candidate_id": cid,
                "best_coefficient": best["c"],
                "strict_shape_error": strict_shape,
                "redesigned_shape_error": best["shape"],
                "weighted_objective": best["obj"],
                "improved_dimensions": improved_dims,
                "decision": decision,
            }
        )

    retained = [r["candidate_id"] for r in rows if r["decision"] == "RETAIN_REDESIGNED"]
    rejected = [r["candidate_id"] for r in rows if r["decision"] == "REJECT_REDESIGNED"]

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2081",
        "stage_id": "S1031",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_FEATURE_WEIGHTED_CANDIDATE_REDESIGN_AUDIT_WITH_TRACE__C3_STILL_OPEN" if ready else "OPEN_FEATURE_WEIGHTED_CANDIDATE_REDESIGN_AUDIT_BLOCKED",
        "depends_on": {"p2080_present": p2080.get("_missing") is None, "preconditions_ready": ready},
        "input_hashes": {
            "p2080_json_sha256": file_sha256(GEN / "p2080_s1030_strict_profile_shape_bridge_audit.json"),
        },
        "redesign_protocol": {
            "d_grid": d_grid,
            "coefficient_bound": 0.2,
            "lambda_l2": lam,
            "feature_weights": {"d1_rmse": 2.0, "d2_rmse": 2.0, "zero_gap": 1.0, "damp_gap": 1.0},
            "retain_rule": "improved_dimensions>=3",
        },
        "redesign_results": {
            "rows": rows,
            "retained_redesigned_candidates": retained,
            "rejected_redesigned_candidates": rejected,
        },
        "c3_gate_update": {
            "C3_feature_weighted_candidate_redesign_audit": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "rows_nonempty": len(rows) > 0,
            "coefficients_within_bound": all(abs(float(r["best_coefficient"])) <= 0.2 for r in rows),
            "decision_domain_ok": all(r["decision"] in {"RETAIN_REDESIGNED", "REJECT_REDESIGNED"} for r in rows),
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
            "# P2081 S1031: feature-weighted candidate redesign audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Retained redesigned: `{retained}`",
            f"- Rejected redesigned: `{rejected}`",
            "",
            "Candidates are redesigned to target derivative/curvature gaps, then re-fit under constraints.",
            "C3 remains OPEN (not discharged).",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
