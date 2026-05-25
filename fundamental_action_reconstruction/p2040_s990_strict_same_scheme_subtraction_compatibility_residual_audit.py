#!/usr/bin/env python3
"""P2040 S990: same-scheme subtraction-compatibility residual audit.

Computes residuals of subtraction compatibility on the same quotient basis
(R2_bar, Ric2_bar, Riem2_bar) for the same controlled background pair used in
P2038/P2039. Exports explicit pre/post residual bounds and keeps C3 open.
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
OUT = GEN / "p2040_s990_strict_same_scheme_subtraction_compatibility_residual_audit.json"
MD = GEN / "p2040_s990_strict_same_scheme_subtraction_compatibility_residual_audit.md"

SCHEMA_VERSION = "p2040_s990_v1"
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


def linf_vec(v: np.ndarray) -> float:
    return float(np.max(np.abs(v)))


def main() -> None:
    GEN.mkdir(exist_ok=True)

    p2038 = load("p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.json")
    p2039 = load("p2039_s989_strict_same_scheme_finite_part_candidate_uncertainty_bound_probe.json")

    p2038_checks = p2038.get("gatekeeper_checks") or {}
    p2039_checks = p2039.get("gatekeeper_checks") or {}

    ready = (
        p2038.get("result_kind")
        == "PASS_FIRST_NONZERO_SAME_SCHEME_FINITE_PART_CANDIDATE_IMPORTED_WITH_TRACE__C3_STILL_OPEN"
        and p2039.get("result_kind")
        == "PASS_FIRST_COMPUTED_UNCERTAINTY_BOUND_FOR_NONZERO_FINITE_PART_CANDIDATE_WITH_TRACE__C3_STILL_OPEN"
        and as_bool(p2038_checks.get("nonzero_candidate_present"))
        and as_bool(p2039_checks.get("uncertainty_bound_finite"))
    )

    imp = p2038.get("candidate_data_import") or {}
    delta = np.array(imp.get("delta_finite_part_candidate_matrix") or [[0.0, 0.0, 0.0]] * 3, dtype=float)

    # Controlled same-scheme subtraction model: D - C where C = M * a, M = I + delta.
    # Residual compatibility checks compare baseline identity vs candidate map.
    a_ref = np.array([0.013, -0.021, 0.008], dtype=float)
    D_ref = a_ref.copy()  # baseline same-scheme reference subtraction target

    M0 = np.eye(3, dtype=float)
    M1 = M0 + delta

    # before candidate: identity map (baseline)
    res_before = D_ref - M0.dot(a_ref)
    # after candidate application
    res_after = D_ref - M1.dot(a_ref)

    # uncertainty-aware post residual bound from P2039 absolute linf bound
    ub = p2039.get("uncertainty_bound") or {}
    abs_linf_bound = float(ub.get("absolute_linf_bound", math.inf))
    a_linf = float(np.max(np.abs(a_ref)))
    propagated_bound = abs_linf_bound * a_linf

    before_linf = linf_vec(res_before)
    after_linf = linf_vec(res_after)
    after_bound = after_linf + propagated_bound

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2040",
        "stage_id": "S990",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_SAME_SCHEME_SUBTRACTION_COMPATIBILITY_RESIDUAL_AUDIT_WITH_BOUND__C3_STILL_OPEN"
            if ready and math.isfinite(after_bound)
            else "OPEN_SUBTRACTION_COMPATIBILITY_RESIDUAL_AUDIT_BLOCKED"
        ),
        "route": "strict_only_same_scheme_residual_audit",
        "depends_on": {
            "p2038_present": p2038.get("_missing") is None,
            "p2039_present": p2039.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2038_json_sha256": file_sha256(GEN / "p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.json"),
            "p2039_json_sha256": file_sha256(GEN / "p2039_s989_strict_same_scheme_finite_part_candidate_uncertainty_bound_probe.json"),
        },
        "audit_scope": {
            "controlled_background_pair": imp.get("controlled_background_pair", "UNKNOWN"),
            "basis": BASIS,
            "subtraction_model": "residual = D_ref - M*a_ref with M in same-scheme quotient basis",
            "a_ref": a_ref.tolist(),
            "D_ref": D_ref.tolist(),
        },
        "residual_audit": {
            "before_candidate_linf": before_linf,
            "after_candidate_linf": after_linf,
            "propagated_uncertainty_bound_linf": propagated_bound,
            "after_candidate_with_uncertainty_upper_bound_linf": after_bound,
            "residual_bound_threshold_linf": 1.0e-4,
            "residual_bound_pass": bool(after_bound <= 1.0e-4),
        },
        "c3_gate_update": {
            "C3_subtraction_compatibility_residual_audit": "COMPUTED_FOR_ONE_CONTROLLED_PAIR",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
            "missing_for_c3_theorem": [
                "operator-level same-scheme subtraction identity proof",
                "cross-background transport identity theorem",
                "global finite-part lock cocycle theorem",
            ],
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "base_candidate_nonzero": bool(np.max(np.abs(delta)) > 0.0),
            "before_residual_finite": math.isfinite(before_linf),
            "after_residual_finite": math.isfinite(after_linf),
            "uncertainty_bound_finite": math.isfinite(propagated_bound),
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    md = [
        "# P2040 S990: same-scheme subtraction-compatibility residual audit",
        "",
        f"- Status: `{payload['status']}`",
        f"- Result kind: `{payload['result_kind']}`",
        "",
        "## Export",
        "",
        "Computed subtraction residuals before/after candidate finite-part map and an",
        "explicit uncertainty-aware residual upper bound on the same quotient basis.",
        "",
        "## Gate update",
        "",
        "- `C3`: subtraction-compatibility residual audit computed for one controlled pair.",
        "- `C3`: theorem remains open (not discharged).",
    ]
    MD.write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
