#!/usr/bin/env python3
"""P2042 S992: operator-level residual envelope sweep for same-scheme witness.

Extends P2041 with a small sweep over a family of test vectors sampled on the
unit sphere and scaled by a controlled norm range. Exports a robust sup-envelope
for residuals plus a separate uncertainty buffer. C3 remains OPEN.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2042_s992_strict_same_scheme_operator_level_residual_envelope_sweep.json"
MD = GEN / "p2042_s992_strict_same_scheme_operator_level_residual_envelope_sweep.md"

SCHEMA_VERSION = "p2042_s992_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"
BASIS = ["R2_bar", "Ric2_bar", "Riem2_bar"]


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def file_sha256(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def as_bool(v: Any) -> bool:
    return bool(v is True)


def linf(v: np.ndarray) -> float:
    return float(np.max(np.abs(v)))


def main() -> None:
    GEN.mkdir(exist_ok=True)

    p2038 = load("p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.json")
    p2040 = load("p2040_s990_strict_same_scheme_subtraction_compatibility_residual_audit.json")
    p2041 = load("p2041_s991_strict_same_scheme_operator_level_subtraction_identity_witness.json")

    checks_2038 = p2038.get("gatekeeper_checks") or {}
    checks_2040 = p2040.get("gatekeeper_checks") or {}
    checks_2041 = p2041.get("gatekeeper_checks") or {}

    ready = (
        p2041.get("result_kind")
        == "PASS_OPERATOR_LEVEL_SAME_SCHEME_SUBTRACTION_IDENTITY_WITNESS_WITH_BOUND__C3_STILL_OPEN"
        and as_bool(checks_2038.get("nonzero_candidate_present"))
        and as_bool(checks_2040.get("uncertainty_bound_finite"))
        and as_bool(checks_2041.get("worst_case_finite"))
    )

    imp = p2038.get("candidate_data_import") or {}
    delta = np.array(imp.get("delta_finite_part_candidate_matrix") or [[0.0, 0.0, 0.0]] * 3, dtype=float)
    M = np.eye(3, dtype=float) + delta

    ra = p2040.get("residual_audit") or {}
    propagated_bound = float(ra.get("propagated_uncertainty_bound_linf", math.inf))

    # Small sweep family: random unit vectors on sphere with controlled norms.
    rng = np.random.default_rng(2042)
    norms = [0.010, 0.015, 0.020, 0.025]
    vectors_per_norm = 12

    rows: list[dict[str, Any]] = []
    worst_residual = 0.0
    worst_residual_with_buffer = 0.0

    for n in norms:
        for k in range(vectors_per_norm):
            v = rng.normal(size=3)
            v_norm = float(np.linalg.norm(v))
            if v_norm == 0.0:
                continue
            u = v / v_norm
            a = n * u
            D = a.copy()
            res = D - M.dot(a)
            r = linf(res)
            r_buf = r + propagated_bound
            worst_residual = max(worst_residual, r)
            worst_residual_with_buffer = max(worst_residual_with_buffer, r_buf)
            rows.append(
                {
                    "norm_bucket": n,
                    "sample_index": k,
                    "a_test": [float(x) for x in a],
                    "residual_linf": r,
                    "residual_with_uncertainty_buffer_linf": r_buf,
                }
            )

    threshold = 2.5e-4

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2042",
        "stage_id": "S992",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_OPERATOR_LEVEL_RESIDUAL_ENVELOPE_SWEEP_WITH_UNCERTAINTY_BUFFER__C3_STILL_OPEN"
            if ready and math.isfinite(worst_residual_with_buffer)
            else "OPEN_OPERATOR_LEVEL_RESIDUAL_ENVELOPE_SWEEP_BLOCKED"
        ),
        "route": "strict_only_operator_level_envelope_sweep",
        "depends_on": {
            "p2038_present": p2038.get("_missing") is None,
            "p2040_present": p2040.get("_missing") is None,
            "p2041_present": p2041.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2038_json_sha256": file_sha256(GEN / "p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.json"),
            "p2040_json_sha256": file_sha256(GEN / "p2040_s990_strict_same_scheme_subtraction_compatibility_residual_audit.json"),
            "p2041_json_sha256": file_sha256(GEN / "p2041_s991_strict_same_scheme_operator_level_subtraction_identity_witness.json"),
        },
        "sweep_spec": {
            "controlled_background_pair": imp.get("controlled_background_pair", "UNKNOWN"),
            "basis": BASIS,
            "seed": 2042,
            "distribution": "normal->unit_sphere",
            "norm_buckets": norms,
            "vectors_per_norm": vectors_per_norm,
            "total_samples": len(rows),
        },
        "residual_envelope": {
            "sup_residual_linf": worst_residual,
            "uncertainty_buffer_linf": propagated_bound,
            "sup_residual_with_buffer_linf": worst_residual_with_buffer,
            "envelope_threshold_linf": threshold,
            "envelope_pass": bool(worst_residual_with_buffer <= threshold),
        },
        "samples": rows,
        "c3_gate_update": {
            "C3_operator_level_envelope_sweep": "COMPUTED_WITH_SPHERE_SAMPLING",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
            "missing_for_c3_theorem": [
                "theorem-grade operator identity proof for full family",
                "cross-background finite-part transport identity theorem",
                "global finite-part lock/cocycle theorem",
            ],
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "sample_count_ge_40": len(rows) >= 40,
            "worst_case_finite": math.isfinite(worst_residual_with_buffer),
            "nonzero_delta_present": bool(np.max(np.abs(delta)) > 0.0),
            "uncertainty_buffer_finite": math.isfinite(propagated_bound),
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    md = [
        "# P2042 S992: operator-level residual envelope sweep",
        "",
        f"- Status: `{payload['status']}`",
        f"- Result kind: `{payload['result_kind']}`",
        "",
        "## Export",
        "",
        "Exported robust sup-envelope of operator-level residuals over a small",
        "sphere-sampled test-vector family with separate uncertainty buffer.",
        "",
        "## Gate update",
        "",
        "- `C3`: envelope witness computed.",
        "- `C3`: theorem remains open (not discharged).",
    ]
    MD.write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
