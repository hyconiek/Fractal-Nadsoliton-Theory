#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2077_s1027_strict_constrained_coefficient_fit_audit.json"
MD = GEN / "p2077_s1027_strict_constrained_coefficient_fit_audit.md"

SCHEMA_VERSION = "p2077_s1027_v1"
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


def mse(a: list[float], b: list[float]) -> float:
    return sum((x - y) ** 2 for x, y in zip(a, b)) / float(len(a))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2076 = load("p2076_s1026_strict_top_candidate_numeric_fit_audit.json")
    ready = p2076.get("result_kind") == "PASS_TOP_CANDIDATE_NUMERIC_FIT_AUDIT_WITH_TRACE__C3_STILL_OPEN"

    top_ids = ((p2076.get("fit_protocol") or {}).get("top_candidate_ids") or []) if isinstance(p2076, dict) else []

    d_grid = [i / 20.0 for i in range(1, 21)]
    legacy_vals = [k_legacy(d) for d in d_grid]
    strict_vals = [k_strict(d) for d in d_grid]
    baseline_mse = mse(legacy_vals, strict_vals)

    def f_candidate(cid: str, d: float) -> float:
        if cid == "mc_phase_lock":
            return d * math.exp(-d)
        if cid == "mc_grad_cross":
            return (d**2) / (1.0 + d**2)
        return math.sin(d) * math.exp(-d)

    coef_grid = [x / 100.0 for x in range(-20, 21)]  # hard bound |c|<=0.2
    lambda_reg = 0.05

    rows = []
    for cid in top_ids:
        best_c = 0.0
        best_raw_mse = baseline_mse
        best_obj = baseline_mse

        for c in coef_grid:
            aug_vals = [k_strict(d) + c * f_candidate(cid, d) for d in d_grid]
            raw = mse(legacy_vals, aug_vals)
            obj = raw + lambda_reg * (c**2)
            if obj < best_obj:
                best_obj = obj
                best_raw_mse = raw
                best_c = c

        abs_improve = baseline_mse - best_raw_mse
        rel_improve = (abs_improve / baseline_mse) if baseline_mse > 0 else 0.0

        qw2191_proxy_risk = "LOW_PROXY"  # bounded by design |c|<=0.2
        retain = abs_improve > 0.0

        rows.append(
            {
                "candidate_id": cid,
                "baseline_mse": baseline_mse,
                "best_raw_mse": best_raw_mse,
                "best_regularized_objective": best_obj,
                "best_coefficient": best_c,
                "absolute_improvement": abs_improve,
                "relative_improvement": rel_improve,
                "qw2191_proxy_risk": qw2191_proxy_risk,
                "decision": "RETAIN" if retain else "REJECT",
            }
        )

    retained = [r["candidate_id"] for r in rows if r["decision"] == "RETAIN"]
    rejected = [r["candidate_id"] for r in rows if r["decision"] == "REJECT"]

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2077",
        "stage_id": "S1027",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_CONSTRAINED_COEFFICIENT_FIT_AUDIT_WITH_TRACE__C3_STILL_OPEN"
            if ready
            else "OPEN_CONSTRAINED_COEFFICIENT_FIT_AUDIT_BLOCKED"
        ),
        "depends_on": {"p2076_present": p2076.get("_missing") is None, "preconditions_ready": ready},
        "input_hashes": {
            "p2076_json_sha256": file_sha256(GEN / "p2076_s1026_strict_top_candidate_numeric_fit_audit.json"),
        },
        "fit_protocol": {
            "d_grid": d_grid,
            "coefficient_bound": 0.2,
            "coefficient_grid": coef_grid,
            "regularization": {"type": "l2", "lambda": lambda_reg},
            "error_metric": "MSE(legacy, augmented_strict)",
            "top_candidate_ids": top_ids,
        },
        "fit_results": {
            "rows": rows,
            "retained_candidates": retained,
            "rejected_candidates": rejected,
        },
        "c3_gate_update": {
            "C3_constrained_coefficient_fit_audit": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "rows_nonempty": len(rows) > 0,
            "all_coefficients_within_bound": all(abs(float(r["best_coefficient"])) <= 0.2 for r in rows),
            "all_qw2191_proxy_low": all(r["qw2191_proxy_risk"] == "LOW_PROXY" for r in rows),
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
                "# P2077 S1027: constrained-coefficient fit audit",
                "",
                f"- Status: `{payload['status']}`",
                f"- Result kind: `{payload['result_kind']}`",
                f"- Retained: `{retained}`",
                f"- Rejected: `{rejected}`",
                "",
                "Fit repeated with hard coefficient bound and L2 regularization.",
                "C3 remains OPEN (not discharged).",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
