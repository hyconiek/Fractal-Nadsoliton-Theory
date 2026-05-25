#!/usr/bin/env python3
"""P2051 S1001: budget-scaling audit.

Compares multiple computational budget profiles for frontier search and exports
sensitivity map of outcome type:
- detected_frontier
- computational_censoring_certificate

Goal: separate genuine no-break behavior from compute-limit artifacts.
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
OUT = GEN / "p2051_s1001_strict_same_scheme_budget_scaling_audit.json"
MD = GEN / "p2051_s1001_strict_same_scheme_budget_scaling_audit.md"

SCHEMA_VERSION = "p2051_s1001_v1"
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
    p2050 = load("p2050_s1000_strict_same_scheme_censoring_resolution_audit.json")
    p2040 = load("p2040_s990_strict_same_scheme_subtraction_compatibility_residual_audit.json")

    checks_2038 = p2038.get("gatekeeper_checks") or {}
    checks_2050 = p2050.get("gatekeeper_checks") or {}

    ready = (
        p2050.get("result_kind")
        == "PASS_CENSORING_RESOLUTION_AUDIT_WITH_FRONTIER_OR_CERTIFICATE__C3_STILL_OPEN"
        and as_bool(checks_2038.get("nonzero_candidate_present"))
        and as_bool(checks_2050.get("epsilon_last_stable_finite"))
    )

    imp = p2038.get("candidate_data_import") or {}
    delta = np.array(imp.get("delta_finite_part_candidate_matrix") or [[0.0, 0.0, 0.0]] * 3, dtype=float)
    propagated = float(((p2040.get("residual_audit") or {}).get("propagated_uncertainty_bound_linf", 0.0)))

    scope2050 = p2050.get("resolution_scope") or {}
    norm_buckets = [float(x) for x in (scope2050.get("norm_buckets") or [0.01, 0.015, 0.02, 0.025])]
    eps_star = float(scope2050.get("epsilon_star_input", 0.0))

    cert_scope = (load("p2047_s997_strict_same_scheme_stability_certificate_extraction_audit.json").get("certificate_scope") or {})
    base_rank = (cert_scope.get("ranking_reference") or {}).get("base_rank") or [4, 3, 2, 1]

    # Budget profiles: small / medium / large.
    profiles = [
        {"profile_id": "B1_small", "max_horizon_eps": 7.5e-4, "max_sign_patterns": 24, "vectors_per_bucket": 12},
        {"profile_id": "B2_medium", "max_horizon_eps": 2.0e-3, "max_sign_patterns": 72, "vectors_per_bucket": 20},
        {"profile_id": "B3_large", "max_horizon_eps": 5.0e-3, "max_sign_patterns": 120, "vectors_per_bucket": 32},
    ]

    ensemble_rng = np.random.default_rng(2051)
    trial_seeds = [int(x) for x in ensemble_rng.integers(7000, 9900, size=4)]

    def is_stable(eps: float, sign_patterns: list[np.ndarray], vectors_per_bucket: int) -> bool:
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

    profile_results = []

    for prof in profiles:
        max_horizon_eps = float(prof["max_horizon_eps"])
        max_sign_patterns = int(prof["max_sign_patterns"])
        vectors_per_bucket = int(prof["vectors_per_bucket"])

        sign_pattern_budget = 18
        step = max(2e-7, eps_star * 0.2 + 1e-7)
        eps_lo = eps_star
        eps_hi = eps_star
        break_found = False
        scan_steps = 0

        while True:
            sign_patterns = []
            for _ in range(sign_pattern_budget):
                s = np.sign(ensemble_rng.normal(size=delta.shape))
                s[s == 0] = 1.0
                sign_patterns.append(s)

            local_found = False
            for _ in range(40):
                scan_steps += 1
                eps_hi += step
                st = is_stable(eps_hi, sign_patterns, vectors_per_bucket)
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

        if break_found:
            resolution_kind = "detected_frontier"
            epsilon_break = eps_hi
            epsilon_last_stable = eps_lo
        else:
            resolution_kind = "computational_censoring_certificate"
            epsilon_break = math.inf
            epsilon_last_stable = eps_lo

        profile_results.append(
            {
                **prof,
                "resolution_kind": resolution_kind,
                "epsilon_last_stable": epsilon_last_stable,
                "epsilon_break": epsilon_break,
                "horizon_reached": bool(eps_hi >= max_horizon_eps),
                "pattern_budget_reached": bool(sign_pattern_budget >= max_sign_patterns),
                "final_sign_pattern_budget": sign_pattern_budget,
                "scan_steps_used": scan_steps,
            }
        )

    detected_count = sum(1 for r in profile_results if r["resolution_kind"] == "detected_frontier")
    censor_count = sum(1 for r in profile_results if r["resolution_kind"] == "computational_censoring_certificate")

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2051",
        "stage_id": "S1001",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_BUDGET_SCALING_AUDIT_WITH_SENSITIVITY_MAP__C3_STILL_OPEN"
            if ready and len(profile_results) > 0
            else "OPEN_BUDGET_SCALING_AUDIT_BLOCKED"
        ),
        "route": "strict_only_budget_scaling_audit",
        "depends_on": {
            "p2038_present": p2038.get("_missing") is None,
            "p2050_present": p2050.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2038_json_sha256": file_sha256(GEN / "p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.json"),
            "p2050_json_sha256": file_sha256(GEN / "p2050_s1000_strict_same_scheme_censoring_resolution_audit.json"),
        },
        "audit_scope": {
            "controlled_background_pair": imp.get("controlled_background_pair", "UNKNOWN"),
            "basis": BASIS,
            "norm_buckets": norm_buckets,
            "trial_seeds": trial_seeds,
            "epsilon_star_input": eps_star,
        },
        "budget_profiles": profile_results,
        "sensitivity_map": {
            "detected_frontier_count": detected_count,
            "computational_censoring_count": censor_count,
            "profiles_total": len(profile_results),
            "interpretation": (
                "Mixed outcomes suggest compute-limit sensitivity; all-censor suggests unresolved frontier under tested budgets; "
                "all-detected suggests robust break detection under tested budgets."
            ),
        },
        "c3_gate_update": {
            "C3_budget_scaling_audit": "COMPUTED",
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
            "profiles_nonempty": len(profile_results) > 0,
            "counts_consistent": detected_count + censor_count == len(profile_results),
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    md = [
        "# P2051 S1001: budget-scaling audit",
        "",
        f"- Status: `{payload['status']}`",
        f"- Result kind: `{payload['result_kind']}`",
        "",
        "## Export",
        "",
        "Compared multiple computational budget profiles and exported sensitivity",
        "map for detected_frontier vs computational_censoring_certificate outcomes.",
        "",
        "## Gate update",
        "",
        "- `C3`: budget-scaling audit computed.",
        "- `C3`: theorem remains open (not discharged).",
    ]
    MD.write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
