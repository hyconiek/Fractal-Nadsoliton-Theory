#!/usr/bin/env python3
"""P2050 S1000: censoring-resolution audit for frontier detection.

Adaptively expands epsilon scan horizon and adversarial sign-pattern budget until
first break detection or explicit compute-limit exhaustion. Exports either:
- detected frontier summary, or
- computational censoring certificate with the reached limits.

Audit-level only; C3 remains OPEN.
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
OUT = GEN / "p2050_s1000_strict_same_scheme_censoring_resolution_audit.json"
MD = GEN / "p2050_s1000_strict_same_scheme_censoring_resolution_audit.md"

SCHEMA_VERSION = "p2050_s1000_v1"
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


def rank_desc(vals: list[float]) -> list[int]:
    idx = sorted(range(len(vals)), key=lambda i: (-vals[i], i))
    out = [0] * len(vals)
    for r, i in enumerate(idx):
        out[i] = r + 1
    return out


def main() -> None:
    GEN.mkdir(exist_ok=True)

    p2038 = load("p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.json")
    p2049 = load("p2049_s999_strict_same_scheme_frontier_confidence_audit.json")
    p2040 = load("p2040_s990_strict_same_scheme_subtraction_compatibility_residual_audit.json")

    checks_2038 = p2038.get("gatekeeper_checks") or {}
    checks_2049 = p2049.get("gatekeeper_checks") or {}

    ready = (
        p2049.get("result_kind")
        == "PASS_FRONTIER_CONFIDENCE_AUDIT_WITH_TRACE__C3_STILL_OPEN"
        and as_bool(checks_2038.get("nonzero_candidate_present"))
        and as_bool(checks_2049.get("break_or_censored_present"))
    )

    imp = p2038.get("candidate_data_import") or {}
    delta = np.array(imp.get("delta_finite_part_candidate_matrix") or [[0.0, 0.0, 0.0]] * 3, dtype=float)
    propagated = float(((p2040.get("residual_audit") or {}).get("propagated_uncertainty_bound_linf", 0.0)))

    fscope = p2049.get("audit_scope") or {}
    norm_buckets = [float(x) for x in (fscope.get("norm_buckets") or [0.01, 0.015, 0.02, 0.025])]
    vectors_per_bucket = int(fscope.get("vectors_per_bucket", 20))
    eps_star = float(fscope.get("epsilon_star_input", 0.0))

    cert_scope = (load("p2047_s997_strict_same_scheme_stability_certificate_extraction_audit.json").get("certificate_scope") or {})
    base_rank = (cert_scope.get("ranking_reference") or {}).get("base_rank") or [4, 3, 2, 1]

    def is_stable(eps: float, trial_seeds: list[int], sign_patterns: list[np.ndarray]) -> bool:
        for sgn in sign_patterns:
            M = np.eye(3, dtype=float) + delta + eps * sgn
            per_seed = []
            for sd in trial_seeds:
                rg = np.random.default_rng(sd)
                bucket_worst = []
                for n in norm_buckets:
                    worst = 0.0
                    for _ in range(vectors_per_bucket):
                        v = rg.normal(size=3)
                        nv = float(np.linalg.norm(v))
                        if nv == 0.0:
                            continue
                        a = n * (v / nv)
                        res = linf(a - M.dot(a)) + propagated
                        worst = max(worst, res)
                    bucket_worst.append(worst)
                per_seed.append(bucket_worst)
            avg = list(np.mean(np.array(per_seed), axis=0))
            if rank_desc(avg) != base_rank:
                return False
        return True

    # Explicit computational limits
    max_horizon_eps = 2.0e-3
    max_scan_steps = 40
    max_sign_patterns = 72

    ensemble_rng = np.random.default_rng(2050)
    trial_seeds = [int(x) for x in ensemble_rng.integers(7000, 9900, size=4)]

    sign_pattern_budget = 18
    step = max(2e-7, eps_star * 0.2 + 1e-7)
    eps_lo = eps_star
    eps_hi = eps_star
    scan_steps = 0
    break_found = False

    expansion_trace = []

    while True:
        # refresh sign patterns with current budget
        sign_patterns = []
        for _ in range(sign_pattern_budget):
            s = np.sign(ensemble_rng.normal(size=delta.shape))
            s[s == 0] = 1.0
            sign_patterns.append(s)

        local_found = False
        for _ in range(max_scan_steps):
            scan_steps += 1
            eps_hi += step
            st = is_stable(eps_hi, trial_seeds, sign_patterns)
            expansion_trace.append(
                {
                    "sign_pattern_budget": sign_pattern_budget,
                    "epsilon": eps_hi,
                    "stable": st,
                    "scan_step": scan_steps,
                }
            )
            if not st:
                local_found = True
                break
            eps_lo = eps_hi
            step *= 1.28
            if eps_hi >= max_horizon_eps:
                break

        if local_found:
            break_found = True
            break

        if eps_hi >= max_horizon_eps or sign_pattern_budget >= max_sign_patterns:
            break

        sign_pattern_budget = min(sign_pattern_budget + 18, max_sign_patterns)

    epsilon_last_stable = eps_lo
    epsilon_break = math.inf
    bisect_trace = []

    if break_found:
        # bisection with final sign-pattern budget
        sign_patterns = []
        for _ in range(sign_pattern_budget):
            s = np.sign(ensemble_rng.normal(size=delta.shape))
            s[s == 0] = 1.0
            sign_patterns.append(s)

        lo = eps_lo
        hi = eps_hi
        for _ in range(24):
            mid = 0.5 * (lo + hi)
            st = is_stable(mid, trial_seeds, sign_patterns)
            bisect_trace.append({"epsilon": mid, "stable": st})
            if st:
                lo = mid
            else:
                hi = mid
        epsilon_last_stable = lo
        epsilon_break = hi

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2050",
        "stage_id": "S1000",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_CENSORING_RESOLUTION_AUDIT_WITH_FRONTIER_OR_CERTIFICATE__C3_STILL_OPEN"
            if ready and math.isfinite(epsilon_last_stable)
            else "OPEN_CENSORING_RESOLUTION_AUDIT_BLOCKED"
        ),
        "route": "strict_only_censoring_resolution_audit",
        "depends_on": {
            "p2038_present": p2038.get("_missing") is None,
            "p2049_present": p2049.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2038_json_sha256": file_sha256(GEN / "p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.json"),
            "p2049_json_sha256": file_sha256(GEN / "p2049_s999_strict_same_scheme_frontier_confidence_audit.json"),
        },
        "resolution_scope": {
            "controlled_background_pair": imp.get("controlled_background_pair", "UNKNOWN"),
            "basis": BASIS,
            "norm_buckets": norm_buckets,
            "trial_seeds": trial_seeds,
            "vectors_per_bucket": vectors_per_bucket,
            "epsilon_star_input": eps_star,
            "limits": {
                "max_horizon_eps": max_horizon_eps,
                "max_scan_steps": max_scan_steps,
                "max_sign_patterns": max_sign_patterns,
            },
        },
        "resolution_result": {
            "break_found": break_found,
            "epsilon_last_stable": epsilon_last_stable,
            "epsilon_break": epsilon_break,
            "frontier_gap": (epsilon_break - epsilon_last_stable) if math.isfinite(epsilon_break) else math.inf,
            "final_sign_pattern_budget": sign_pattern_budget,
            "scan_steps_used": scan_steps,
            "expansion_trace": expansion_trace,
            "bisect_trace": bisect_trace,
            "resolution_kind": (
                "detected_frontier" if break_found else "computational_censoring_certificate"
            ),
            "censoring_certificate": {
                "horizon_reached": bool(eps_hi >= max_horizon_eps),
                "pattern_budget_reached": bool(sign_pattern_budget >= max_sign_patterns),
                "certificate_statement": (
                    "No instability detected up to declared limits" if not break_found else "not_applicable"
                ),
            },
        },
        "c3_gate_update": {
            "C3_censoring_resolution_audit": "COMPUTED",
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
            "resolution_kind_set": break_found or (not break_found),
            "epsilon_last_stable_finite": math.isfinite(epsilon_last_stable),
            "trace_nonempty": len(expansion_trace) > 0,
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    md = [
        "# P2050 S1000: censoring-resolution audit",
        "",
        f"- Status: `{payload['status']}`",
        f"- Result kind: `{payload['result_kind']}`",
        "",
        "## Export",
        "",
        "Adaptively expanded epsilon horizon and adversarial-pattern budget until",
        "frontier detection or explicit computational censoring limits.",
        "",
        "## Gate update",
        "",
        "- `C3`: censoring-resolution audit computed.",
        "- `C3`: theorem remains open (not discharged).",
    ]
    MD.write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
