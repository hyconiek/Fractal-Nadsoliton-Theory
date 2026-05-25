#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2076_s1026_strict_top_candidate_numeric_fit_audit.json"
MD = GEN / "p2076_s1026_strict_top_candidate_numeric_fit_audit.md"

SCHEMA_VERSION = "p2076_s1026_v1"
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


def mse(values_a: list[float], values_b: list[float]) -> float:
    n = len(values_a)
    return sum((a - b) ** 2 for a, b in zip(values_a, values_b)) / float(n)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2075 = load("p2075_s1025_strict_missing_characteristic_candidate_screening_audit.json")
    ready = p2075.get("result_kind") == "PASS_MISSING_CHARACTERISTIC_CANDIDATE_SCREENING_AUDIT_WITH_TRACE__C3_STILL_OPEN"

    rows = ((p2075.get("screening_table") or {}).get("rows") or []) if isinstance(p2075, dict) else []
    by_priority = sorted(rows, key=lambda r: int(r.get("priority_rank", 99)))
    top_ids = [r.get("candidate_id") for r in by_priority[:2]]

    d_grid = [i / 20.0 for i in range(1, 21)]
    legacy_vals = [k_legacy(d) for d in d_grid]
    strict_vals = [k_strict(d) for d in d_grid]
    baseline_err = mse(legacy_vals, strict_vals)

    # controlled one-parameter augmentations over strict profile
    # K_aug = K_strict + c * f_i(d) where f_i is candidate-dependent
    def f_candidate(cid: str, d: float) -> float:
        if cid == "mc_phase_lock":
            return d * math.exp(-d)
        if cid == "mc_grad_cross":
            return (d**2) / (1.0 + d**2)
        return math.sin(d) * math.exp(-d)

    coef_grid = [x / 100.0 for x in range(-50, 51)]
    fit_rows = []
    for cid in top_ids:
        best_c = 0.0
        best_err = baseline_err
        for c in coef_grid:
            aug_vals = [k_strict(d) + c * f_candidate(cid, d) for d in d_grid]
            err = mse(legacy_vals, aug_vals)
            if err < best_err:
                best_err = err
                best_c = c

        improvement = baseline_err - best_err
        rel_improvement = (improvement / baseline_err) if baseline_err > 0 else 0.0

        # conservative QW-2191 proxy risk: prefer small |c|
        qw2191_risk = "LOW_PROXY" if abs(best_c) <= 0.2 else "ELEVATED_PROXY"
        admit = (improvement > 0.0) and (qw2191_risk == "LOW_PROXY")

        fit_rows.append(
            {
                "candidate_id": cid,
                "baseline_mse": baseline_err,
                "best_mse": best_err,
                "best_coefficient": best_c,
                "absolute_improvement": improvement,
                "relative_improvement": rel_improvement,
                "qw2191_proxy_risk": qw2191_risk,
                "decision": "RETAIN" if admit else "REJECT",
            }
        )

    retained = [r for r in fit_rows if r["decision"] == "RETAIN"]

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2076",
        "stage_id": "S1026",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_TOP_CANDIDATE_NUMERIC_FIT_AUDIT_WITH_TRACE__C3_STILL_OPEN"
            if ready
            else "OPEN_TOP_CANDIDATE_NUMERIC_FIT_AUDIT_BLOCKED"
        ),
        "depends_on": {"p2075_present": p2075.get("_missing") is None, "preconditions_ready": ready},
        "input_hashes": {
            "p2075_json_sha256": file_sha256(GEN / "p2075_s1025_strict_missing_characteristic_candidate_screening_audit.json"),
        },
        "fit_protocol": {
            "legacy_kernel": "alpha_geo*cos(omega*d+phi)/(1+beta_tors*d)",
            "strict_kernel": "cos(omega*d+phi)/(1+beta*d^eta)",
            "d_grid": d_grid,
            "coefficient_grid": coef_grid,
            "error_metric": "MSE(legacy, augmented_strict)",
            "top_candidate_ids": top_ids,
        },
        "fit_results": {
            "rows": fit_rows,
            "retained_candidates": [r["candidate_id"] for r in retained],
            "rejected_candidates": [r["candidate_id"] for r in fit_rows if r["decision"] == "REJECT"],
        },
        "c3_gate_update": {
            "C3_top_candidate_numeric_fit_audit": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
            "missing_for_theorem_grade": [
                "proof-level legacy<->strict bridge/non-bridge theorem",
                "global selector-safe closure beyond QW-2191 proxy risk",
                "full tensorial EOM closure on declared background family",
            ],
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "fit_rows_nonempty": len(fit_rows) > 0,
            "improvement_computed": all("absolute_improvement" in r for r in fit_rows),
            "has_retain_or_reject_decisions": all(r.get("decision") in {"RETAIN", "REJECT"} for r in fit_rows),
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
            "# P2076 S1026: top-candidate numeric fit audit",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Top candidates fitted: `{top_ids}`",
            f"- Retained: `{payload['fit_results']['retained_candidates']}`",
            f"- Rejected: `{payload['fit_results']['rejected_candidates']}`",
            "",
            "Controlled numeric fit executed with conservative QW-2191 proxy screening.",
            "C3 remains OPEN (not discharged).",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
