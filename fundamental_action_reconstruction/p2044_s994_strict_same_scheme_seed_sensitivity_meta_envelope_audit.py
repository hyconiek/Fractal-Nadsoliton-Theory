#!/usr/bin/env python3
"""P2044 S994: stress-test seed sensitivity meta-envelope audit.

Repeats the P2043-style envelope stratification across several RNG seeds and
exports a meta-envelope (sup across seeds) plus between-seed confidence bounds.
This is still audit-level and does not discharge C3.
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
OUT = GEN / "p2044_s994_strict_same_scheme_seed_sensitivity_meta_envelope_audit.json"
MD = GEN / "p2044_s994_strict_same_scheme_seed_sensitivity_meta_envelope_audit.md"

SCHEMA_VERSION = "p2044_s994_v1"
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


def run_seed_sweep(seed: int, M: np.ndarray, propagated_bound: float) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    norms = [0.010, 0.015, 0.020, 0.025]
    vectors_per_norm = 12

    by_bucket: dict[float, list[float]] = {n: [] for n in norms}
    worst = 0.0
    rows = 0

    for n in norms:
        for _ in range(vectors_per_norm):
            v = rng.normal(size=3)
            nv = float(np.linalg.norm(v))
            if nv == 0.0:
                continue
            a = n * (v / nv)
            D = a.copy()
            r = linf(D - M.dot(a)) + propagated_bound
            by_bucket[n].append(r)
            worst = max(worst, r)
            rows += 1

    bucket_sup = {str(k): float(max(v)) for k, v in by_bucket.items()}
    return {
        "seed": seed,
        "total_samples": rows,
        "bucket_sup_residual_with_buffer_linf": bucket_sup,
        "global_sup_residual_with_buffer_linf": float(worst),
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)

    p2038 = load("p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.json")
    p2043 = load("p2043_s993_strict_same_scheme_stratified_residual_envelope_audit.json")

    checks_2038 = p2038.get("gatekeeper_checks") or {}
    checks_2043 = p2043.get("gatekeeper_checks") or {}

    ready = (
        p2043.get("result_kind")
        == "PASS_STRATIFIED_RESIDUAL_ENVELOPE_AUDIT_WITH_BOOTSTRAP_CI__C3_STILL_OPEN"
        and as_bool(checks_2038.get("nonzero_candidate_present"))
        and as_bool(checks_2043.get("bootstrap_ci_finite"))
    )

    imp = p2038.get("candidate_data_import") or {}
    delta = np.array(imp.get("delta_finite_part_candidate_matrix") or [[0.0, 0.0, 0.0]] * 3, dtype=float)
    M = np.eye(3, dtype=float) + delta

    propagated_bound = float(((load("p2040_s990_strict_same_scheme_subtraction_compatibility_residual_audit.json").get("residual_audit") or {}).get("propagated_uncertainty_bound_linf", math.inf)))

    seeds = [2037, 2042, 2043, 2051, 2063]
    per_seed = [run_seed_sweep(s, M, propagated_bound) for s in seeds]

    per_seed_sup = np.array([row["global_sup_residual_with_buffer_linf"] for row in per_seed], dtype=float)
    meta_sup = float(np.max(per_seed_sup)) if len(per_seed_sup) else math.inf

    # CI between seeds (empirical quantile CI on seed-level sup values)
    if len(per_seed_sup) > 0:
        ci_low = float(np.quantile(per_seed_sup, 0.025))
        ci_high = float(np.quantile(per_seed_sup, 0.975))
    else:
        ci_low = math.inf
        ci_high = math.inf

    bucket_keys = ["0.01", "0.015", "0.02", "0.025"]
    meta_bucket_sup = {
        k: float(max(row["bucket_sup_residual_with_buffer_linf"][k] for row in per_seed))
        for k in bucket_keys
    }

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2044",
        "stage_id": "S994",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_SEED_SENSITIVITY_META_ENVELOPE_AUDIT_WITH_INTERSEED_CI__C3_STILL_OPEN"
            if ready and math.isfinite(meta_sup)
            else "OPEN_SEED_SENSITIVITY_META_ENVELOPE_AUDIT_BLOCKED"
        ),
        "route": "strict_only_seed_sensitivity_meta_envelope",
        "depends_on": {
            "p2038_present": p2038.get("_missing") is None,
            "p2043_present": p2043.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2038_json_sha256": file_sha256(GEN / "p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.json"),
            "p2043_json_sha256": file_sha256(GEN / "p2043_s993_strict_same_scheme_stratified_residual_envelope_audit.json"),
        },
        "seed_sensitivity_spec": {
            "controlled_background_pair": imp.get("controlled_background_pair", "UNKNOWN"),
            "basis": BASIS,
            "seeds": seeds,
            "vectors_per_norm": 12,
            "norm_buckets": [0.010, 0.015, 0.020, 0.025],
        },
        "per_seed_results": per_seed,
        "meta_envelope": {
            "meta_sup_residual_with_buffer_linf": meta_sup,
            "meta_bucket_sup_residual_with_buffer_linf": meta_bucket_sup,
            "interseed_ci_level": 0.95,
            "interseed_ci_low": ci_low,
            "interseed_ci_high": ci_high,
            "interseed_ci_width": (ci_high - ci_low) if math.isfinite(ci_low) and math.isfinite(ci_high) else math.inf,
        },
        "c3_gate_update": {
            "C3_seed_sensitivity_meta_envelope_audit": "COMPUTED",
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
            "seed_count_ge_5": len(seeds) >= 5,
            "meta_sup_finite": math.isfinite(meta_sup),
            "interseed_ci_finite": math.isfinite(ci_low) and math.isfinite(ci_high),
            "nonzero_delta_present": bool(np.max(np.abs(delta)) > 0.0),
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    md = [
        "# P2044 S994: seed sensitivity meta-envelope audit",
        "",
        f"- Status: `{payload['status']}`",
        f"- Result kind: `{payload['result_kind']}`",
        "",
        "## Export",
        "",
        "Repeated stratified-envelope sweep across multiple RNG seeds and exported",
        "meta-envelope (sup over seeds) with inter-seed CI.",
        "",
        "## Gate update",
        "",
        "- `C3`: seed-sensitivity meta-envelope audit computed.",
        "- `C3`: theorem remains open (not discharged).",
    ]
    MD.write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
